#!/usr/bin/env perl
#
# Quilter.pl - Alex Nord - 2016
#
# USAGE: $ perl Quilter.pl <protein.fa> <genome.fa> <index.gtf> <species> [OPT.s]
#
# ABOUT: This program identifies high-similarity alignments of protein
#        sequences to a genome using a gtf-formatted index file and
#        specification of the species that the genome corresponds to.
#
use warnings;
use strict;
use POSIX;
use Time::HiRes;
use File::Basename;
use lib dirname (__FILE__);
use Cwd;
use DiagonalSets;



sub MAX;
sub MIN;
sub TranslateCodon;
sub GetStopPoints;
sub BuildGeneIndex;
sub AddExon;
sub GetChromosomeLengths;
sub RunFastDiagonals;
sub FindFullAlignments;
sub ListCodonCenters;
sub ExonAssistedSPALN;
sub BLATAssistedSPALN;
sub ParseSPALNOutput;
sub AdjustSPALNOutput;
sub SelenocysteineCheck;
sub LooksRepetitive;
sub GenSortIndex;




####################
####            ####
####   SCRIPT   ####
####            ####
####################


# Bail if we have insufficient arguments
my $expectedArgs = 4;
if (@ARGV < $expectedArgs) {
    print "\n\n";
    print "  USAGE:  ./Quilter.pl  <isoforms>  <genome>  <gtf index>  <species>  [options]\n\n"; 
    print "  WHERE:  <isoforms>  is a FASTA-formatted file containing one or more\n";
    print "                      amino-acid encoded isoforms.\n";
    print "          <genome>    is a FASTA-formatted file containing a DNA-encoded\n";
    print "                      genome.\n";
    print "          <gtf index> is a gtf-formatted index corresponding to the\n";
    print "                      provided genome.\n";
    print "          <species>   is the (case-insensitive) species name that the program\n";
    print "                      will match in the gtf index.\n\n";
    print "  OPT.s:  -n <int>    Select number of CPU cores to use.\n";
    print "          -o <string> Write output to a file with a specific name.\n";
    print "          -v          Verbose output.\n";
    die "\n\n";
}


# Option variables.
my $debug   = 0; # Print out a TON of stuff -- best for small cases.
my $timing  = 0; # Collect timing data
my $overw   = 0; # Don't overwrite (default)
my $resdir  = 0; # Use a user provided results directory
my $CPUs    = 2; # Default: use 2 cores
my $verbose = 0; # Use verbose printing
my $spalner = 0; # Are we just running this as a means for getting SPALN output?
my $noFD    = 0; # WHAT?! (Only run SPALN, not FD)


# Stop points (per CPU, relative to protein DB)
my @StopPoints;

# Check for additional options
my $i = $expectedArgs;
while ($i < @ARGV) {
    my $opt = lc($ARGV[$i]);
    if ($opt =~ /\-debug$/) {
	$debug = 1;
    } elsif ($opt =~ /\-overwrite$/) {
	$overw = 1;
    } elsif ($opt =~ /\-o$/) {
	$i++;
	$resdir = $i; # We just want the argument index
	$overw  = 1;
    } elsif ($opt =~ /\-setcwd$/) {
	$i++;
	chdir($ARGV[$i]);
    } elsif ($opt =~ /\-v$/) {
	$verbose = 1;
    } elsif ($opt =~ /\-n$/) {
	$i++;
	$CPUs = int($ARGV[$i]);
	if ($CPUs <= 0) {
	    print "\n\tUnsupported number of CPUs requested ($CPUs)\n";
	    print "\tReverting to default (2)\n\n";
	    $CPUs = 2;
	}
    } elsif ($opt =~ /\-nplus$/) {
	$i++;
	$CPUs = int($ARGV[$i]);
	if ($CPUs <= 0) {
	    print "\n\tUnsupported number of CPUs requested ($CPUs)\n";
	    print "\tReverting to default (2)\n\n";
	    $CPUs = 2;
	} else {
	    for (my $j=0; $j<$CPUs; $j++) {
		$i++;
		push(@StopPoints,$ARGV[$i]);
		if ($j && $StopPoints[$j] <= $StopPoints[$j-1]) {
		    die "\n  -nplus stop point inconsistency detected (around $j: $StopPoints[$j-1] >= $StopPoints[$j])\n\n";
		}
	    }
	}
    } elsif ($opt =~ /\-spaln$/) {
	$spalner = 1;
    } elsif ($opt =~ /\-time$/) {
	$timing = 1;
    } elsif ($opt =~ /\-fast$/) {
	$noFD = 1;
    } else {
	print "\n\tUnrecognized option '$ARGV[$i]' ignored.\n\n";
    }
    $i++;
}


# Do we need to come up with our stop points?
if (!scalar(@StopPoints)) {
    my $stoppoints_ref;
    ($CPUs,$stoppoints_ref) = GetStopPoints($ARGV[0],lc($ARGV[3]),$CPUs);
    @StopPoints = @{$ARGV[0]};
}


# Let the timing begin!
my @TimingData;
my ($timeA,$timeB);
if ($timing) {
    # The most-importantest stuff
    $TimingData[0]  = [Time::HiRes::gettimeofday()]; ###### START
    $TimingData[1]  = 0.0; # Sum of FastDiagonals-associated times
    $TimingData[2]  = 0.0; # Sum of FastDiagonals system (only) times
    $TimingData[3]  = 0;   # Number of FastDiagonals system calls
    $TimingData[4]  = 0.0; # Sum of SPALN-associated times
    $TimingData[5]  = 0.0; # Sum of SPALN system (only) times
    $TimingData[6]  = 0;   # Number of SPALN system calls
    $TimingData[7]  = 0.0; # BLAT system time
    $TimingData[8]  = 0.0; # BLAT+SPALN-associated time    
    $TimingData[9]  = 0.0; # sum of SPALN system (only) times from BLAT
    $TimingData[10] = 0;   # Number of SPALN system calls from BLAT
    $TimingData[11] = 0;   ################################## END
    # Less most-importantest stuff
    $TimingData[12] = 0.0; # Indexing time
}


# Check if 'SPALN' is installed (and, if so, check if BLAST is happy)
open(my $spalncheck,'which spaln |');
my $spaln = <$spalncheck>;
close($spalncheck);


# Well... is it? (Also, make a file to catch SPALN deaths)
my %ChrLengths;
if (!$spaln) {
    die "\n  ERROR:  Quilter can't run '--fast' option without SPALN installation\n\n" if ($noFD);
    $spaln = 0;
    print "\n  Warning:  Without installation of 'spaln' accuracy may suffer.\n\n";
} else {

    # Do we know how long each of these chromosomes are?
    my $chrlenindex = $ARGV[1].'.chr_length_index';
    my $chrlenfile;
    if (-e $chrlenindex) {
	
	# Rip chromosome lengths
	open($chrlenfile,'<',$chrlenindex);
	while (my $chrline = <$chrlenfile>) {
	    $chrline =~ s/\n|\r//g;
	    next if (!$chrline); # Catch garbage
	    $chrline =~ /(\S+) (\d+)/;
	    $ChrLengths{$1} = $2;
	}
	close($chrlenfile);

    } else {

	# Identify the lengths of the chromosomes in the genome file
	print "  Indexing chromosome lengths... ";
	my $chrlenref = GetChromosomeLengths($ARGV[1]);
	%ChrLengths   = %{$chrlenref};

	# Write these out to a file, so we don't have to do this again
	open($chrlenfile,'>',$chrlenindex);
	foreach $i (sort keys %ChrLengths) {
	    print $chrlenfile "$i $ChrLengths{$i}\n";
	}
	close($chrlenfile);

	print "done\n";
    }

    # Make sure that those damn environment variables get set!
    $spaln =~ s/\n|\r//g;
    $spaln =~ s/spaln$//;
    $ENV{'ALN_TAB'} = $spaln.'../table';
    $ENV{'ALN_DBS'} = $spaln.'../seqdb';

    # Now become a boolean!
    $spaln = 1;

}

# Next, check if we can blat things up in here
open(my $blatcheck,'which blat |');
my $blat = $blatcheck;
close($blatcheck);


# If we need an easel index on the genome, then we need an easel index on the genome
if (!(-e $ARGV[1].'.ssi')) {
    print "  Building easel index... ";
    if (system("esl\-sfetch \-\-index $ARGV[1]")) { die "\n\tFailed to build easel index on $ARGV[1]\n\n"; }
    print "done\n";
}


# Report that we're getting started
my $progmsg = '  Quilter.pl:  Preparing index for '.lc($ARGV[3]);
while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
print "$progmsg\r" if (!$verbose);


# Grab the time
$timeA = [Time::HiRes::gettimeofday()] if ($timing);


# Construct an index on the genes in the isoform file, extracted from the gtf file.
# Additionally, get a list of all genes that did not have any matches in the gtf file.
my $reqSpecies = uc($ARGV[3]);
my %GeneIndex;
my %FamChrs;
my %Blacklist;
my $AutoBLAT = 0;
if ($ARGV[2] ne '-') {
    my ($indexRef,$famchrsRef,$blacklistRef) 
	= BuildGeneIndex($ARGV[0],$ARGV[2],$reqSpecies,$debug);
    %GeneIndex = %{$indexRef};
    %FamChrs   = %{$famchrsRef};
    %Blacklist = %{$blacklistRef};
} elsif ($blat) {
    $AutoBLAT = 1;
} else {
    die "\n  ERROR: Without blat species '$reqSpecies' requires valid GTF index.\n\n";
}


# Record how long it took to build our index
$TimingData[12] = Time::HiRes::tv_interval($timeA) if ($timing);


# Now that we have our indices for each of the genes, we can actually get to work!
# This is done by reading in an isoform from the FASTA file, finding all of its
# indices (searching the gene name), and trying to find full extensions along some
# chromosome.


# We want a temporary folder to stuff everything in
my $foldername = 'tmp_Quilter/';
if (!$overw) {
    $i = 1;
    while (-e $foldername) {
	$foldername = 'tmp_Quilter_'.$i.'/';
	$i++;
    }
    if (system("mkdir $foldername")) { die "\n\tFailed to generate temporary folder '$foldername'\n\n"; }
} else {
    if (-e $foldername && system("rm $foldername/\*; rmdir $foldername")) { 
	die "\n\tFailed to clear existing temporary folder '$foldername'\n\n"; 
    }
    if (system("mkdir $foldername")) { die "\n\tFailed to generate temporary folder '$foldername'\n\n"; }
}


# To guarantee that our threads don't overlap, we try to find a
# good mid-point in the isoform file to split around.


# Non-0 threads report results to a file for thread 0
my $resultsfilehead = $foldername.'temp_results.'; 
my $resultsfile;


# Filename used to track progress of threads
my $progressbase = $foldername.'quilter.thread_progress.';


# Thread 0 spits out all sorts of happy friends, just like in the nature documentaries
my $processes = 1;
my $threadID  = 0; # In case CPUs = 1
my $pid;
while ($processes < $CPUs) {

    # Create a new thread 
    if ( $pid = fork ) { 
	# If things went wrong, get out while you can
	if (not defined $pid) { die "\n\tFork failed\n\n"; }
	# Giving the 'master' thread an ID of '0'
	$threadID = 0;
    } else { 
	# Child process picks up the last number the parent counted
	$threadID = $processes; 
	last;
    }
    $processes++;
}


# The name of various temporary files we'll be writing to
my $nuclfilename    = $foldername.'nucl_tempfile_'.$threadID.'.Quilter.fa';
my $protbasename    = $foldername.'prot_tempfile_'.$threadID.'-';
my $protfile_ext    = '.Quilter.fa';
my $hitfilename     = $foldername.'hit_tempfile_'.$threadID.'.Quilter.out';
my $missfilename    = $foldername.'miss_tempfile_'.$threadID.'.Quilter.out';
my $blatfilename    = $foldername.'blat_tempfile_'.$threadID.'.Quilter.fa'; 
my $resultsfilename = $resultsfilehead.$threadID.'.Quilter.out';


# Each thread writes its progress (number of completed families)
# to a file containing
my $progtime = time();
my $progfilename = $progressbase.$threadID;
my $num_complete = 0;


# We want to make sure we're working with a clean slate
system("rm $nuclfilename") if (-e $nuclfilename);
#system("rm $protfilename") if (-e $protfilename);
system("rm $hitfilename")  if (-e $hitfilename);
system("rm $missfilename") if (-e $missfilename);
system("rm $blatfilename") if (-e $blatfilename);


# A log file used to track how well spaln performs
open(my $SpalnLog,'>',$foldername.'spaln_logfile_'.$threadID.'.Quilter.out');


# A file to hold all our hits and a file to hold our unhit sequence info.
open(my $HitFile,'>',$hitfilename);
open(my $MissFile,'>',$missfilename);


# Used to track how things went down across the whole database
my @HitStats;
$HitStats[0] = 0; # Not indexed
$HitStats[1] = 0; # Wrong species
$HitStats[2] = 0; # Indexed, but no hits from FastDiagonals.c
$HitStats[3] = 0; # FastDiagonals.c gave us enough to score
$HitStats[4] = 0; # SPALN saves the day
$HitStats[5] = 0; # Failure -- not enough to score
$HitStats[6] = 0; # How many were just 1 or 2 end AAs off?
$HitStats[7] = 0; # How many needed blast assistance?

# Repetitiveness Confusion Matrix
my @RepConf;
$RepConf[0][0] = 0;
$RepConf[0][1] = 0;
$RepConf[1][0] = 0;
$RepConf[1][1] = 0;


# Where in the protein database do we want to start / stop?
my $startpoint = 1;
if ($threadID) { $startpoint = $StopPoints[$threadID-1]; }
my $stoppoint  = $StopPoints[$threadID];


# Open the file containing all of the isoforms
open (my $isoformfile,'<',$ARGV[0]) || die "\n\tCould not open isoform(s) file '$ARGV[0]'\n\n";
my $lineNum = 0;


# Walk your way up to the starting line (<= because we want to prime with new seqname)
my $line;
while ($lineNum < $startpoint) {
    $line = <$isoformfile>;
    $lineNum++;
}


# RAD! Time to find happy alignments!
my $max_gene_seqs = 0;
my $numHits = 0;
my $numIsos = 0;
while (!eof($isoformfile) && $lineNum < $stoppoint) {
    
    # Eat any comment lines
    while ($lineNum < $stoppoint && $line !~ m/^\>/) {
	$line = <$isoformfile>;
	$lineNum++;
    }
    #last if ($lineNum >= $stoppoint || eof($isoformfile)); # <-- Should be obsolete

    # Strip away line breaks
    $line =~ s/\r|\n//g;

     # Record sequence name information
    $line =~ /^\>[^\|]+\|[^\|]+\|([^\|]+)\|\S+\|([^\|]+)$/;
    if (!($1 && $2)) {
	die "$reqSpecies: $threadID:  $line ($lineNum: $startpoint,$stoppoint)\n";
    }
    my $gene      = uc($2);
    my $species   = uc($1);
    my $nameline  = $line;
    $nameline     =~ s/\n|\r//g;
    $nameline     =~ s/\>//;
    
    
    # If this is the wrong species, skip it << Should be obsolete-ish
    if (uc($species) ne $reqSpecies) {
	
	# Record that we had a wrong species
	$HitStats[1]++;
	
	# Next up!
	$line = <$isoformfile>;
	$lineNum++;
	next;
	
    }

    
    # We write each isoform to a file, which is then passed to FastDiagonals.
    # We'll also hold onto a copy of the protein sequence for later reference.
    #
    # Because gene families should be clustered within species,
    # we can press forward until we've either (i.) hit the end of
    # this family, or (ii.) hit the end of the file.
    #
    # Note that we ought to still be checking species. <-- eventually, maybe...
    #
    my $num_gene_seqs = 0;
    my @ProtSeqNames;
    my @ProtFileNames;
    my @ProtSequences;
    my @ProtSeqLengths;
    my $next_gene = $gene;

    while ($next_gene eq $gene && !eof($isoformfile)) {

	$num_gene_seqs++;
	$numIsos++;
	$max_gene_seqs = $num_gene_seqs if ($num_gene_seqs > $max_gene_seqs);

	my $protfilename = $protbasename.$num_gene_seqs.$protfile_ext;
	push(@ProtFileNames,$protfilename);

	open(my $protfile,'>',$protfilename);
	print $protfile "$line\n";

	$line =~ s/^\>//;
	push(@ProtSeqNames,$line);

	my $Protein = '';
	$line =  <$isoformfile>;
	$lineNum++;

	while (!eof($isoformfile) && $line !~ /^\>/) {
	    $line =~ s/\n|\r//g;
	    print $protfile "$line\n";
	    $Protein = $Protein.$line;
	    $line =  <$isoformfile>;
	    $lineNum++;
	}
	if (eof($isoformfile) && $line) { # Special catch for final line
	    $line =~ s/\n|\r//g;
	    print $protfile "$line\n";
	    $Protein = $Protein.$line;
	    $lineNum++;
	}
	my $ProtLength = length($Protein);
	close($protfile);

	push(@ProtSequences,$Protein);
	push(@ProtSeqLengths,$ProtLength);

	if (!eof($isoformfile)) {
	    $line =~ /\|([^\|]+)\s*$/;
	    $next_gene = uc($1);
	}

    }    

    # How we'll track each protein sequence's progress in finding a happy hit
    my @CurrentStats;
    my @HitStats;
    my @Repetitive;
    for ($i=0; $i<$num_gene_seqs; $i++) {
	$CurrentStats[$i] = 0;
	$Repetitive[$i]   = LooksRepetitive($ProtFileNames[$i]) unless ($spalner);
    }
    
    
    # If this is a blacklisted gene or we're "AutoBLAT"-ing, cut to the
    # chase with using BLAST
    if (($AutoBLAT || $Blacklist{$gene}) && $spaln && $blat) {
	for (my $i=0; $i<$num_gene_seqs; $i++) {
	    if (system("cat $ProtFileNames[$i] >> $blatfilename")) { die "\n  Concatenation of seq. failed\n\n"; }
	}
    } else {
	
	# In the event of multiple chromosomes having strong hits,
	# we'll grab the best `fast diagonals' hit (trusting manual
	# annotation over SPALN) or best `SPALN' hit if we can't
	# get a good fast diags hit.
	my @BestFastScores;
	my @BestFastStrings;
	my @BestSpalnScores;
	my @BestSpalnStrings;
	for (my $i=0; $i<$num_gene_seqs; $i++) {
	    $BestFastScores[$i]   = 0;
	    $BestFastStrings[$i]  = '';
	    $BestSpalnScores[$i]  = 0;
	    $BestSpalnStrings[$i] = '';
	}
	my $next_score;
	my $next_string;
	
	# Let's do some actual work!  Go through the gene index, considering each
	# chromosome that the gene has been indexed to, creating 'exon' objects for
	# each pair of start and end points, and then trying to find a full match to
	# the protein.
	my $chr_list_str = $FamChrs{$gene};
	$chr_list_str =~ s/^\#//;
	$chr_list_str =~ s/\#$//;

	foreach my $chromosomeName (sort(split(/\#/,$chr_list_str))) {


	    # We don't want to bother trying to align to a chromosome that
	    # doesn't correspond to any of the sequences in the genome file
	    # that we've been given, but we're going to miss quite a few if
	    # we don't correct for 'revcomp' indicator
	    my $raw_chr_name = $chromosomeName;
	    $raw_chr_name =~ s/\[revcomp\]//g;
	    next if (!$ChrLengths{$raw_chr_name});

	    
	    # Add this chromosome to the chromosome set (hash).
	    my @Chromosomes;
	    for ($i=0; $i<$num_gene_seqs; $i++) {
		my $Chromosome = New DiagonalSet($chromosomeName);
		push(@Chromosomes,$Chromosome);
	    }
	    

	    # Time how long we're playing in the land of Fast Diagonals
	    $timeA = [Time::HiRes::gettimeofday()] if ($timing && !$noFD);


	    # Iterate over all entries for this gene on this chromosome.
	    if (!$noFD) {

		# NOTE ON GENEINDEX ENTRIES (2018/02/26)
		#
		# These will be sorted according to direction, so a forward GeneIndex
		# entry will look like '3-4,6-7' whereas a revcomp entry will look
		# like '7-6,4-3'
		#

		# Find the start and end of the whole shebang -- note that these will
		# respect revcomp-iness
		my $gene_index_list = $GeneIndex{$gene.'|'.$chromosomeName};

		# Whether or not we're considering the reverse complement will be
		# identifiable by whether or not we appended '[revcomp]' to the
		# chromosome name.
		my $revcomp = 0;
		if ($chromosomeName =~ /\[revcomp\]/) {
		    $revcomp = 1;
		    $chromosomeName =~ s/\[revcomp\]//;
		}
		
		# Find the highest and lowest nucleotide indices in the coding 
		# region for this gene
		my $minNucl = 100000000000; # 100 GB genome? unlikely
		my $maxNucl = -1;
		my ($start,$end);
		foreach my $start_end (split(/\,/,$gene_index_list)) {
		    $start_end =~ /^(\d+)\-(\d+)$/;
		    $start = $1;
		    $end   = $2;
		    if ($revcomp) {
			$minNucl = MIN($end,$minNucl);
			$maxNucl = MAX($start,$maxNucl);
		    } else {
			$minNucl = MIN($start,$minNucl);
			$maxNucl = MAX($end,$maxNucl);
		    }
		}
		if ($revcomp) {
		    $start = $maxNucl;
		    $end = $minNucl;
		} else {
		    $start = $minNucl;
		    $end = $maxNucl;
		}

		# The system call to generate a FASTA file using the given index info.
		my $eslsfetchCmd;
		$eslsfetchCmd = 'esl-sfetch -c '.$start.'..'.$end;               # Range (revcomp-ready)
		$eslsfetchCmd = $eslsfetchCmd.' -o '.$nuclfilename;              # Output file
		$eslsfetchCmd = $eslsfetchCmd.' '.$ARGV[1].' '.$chromosomeName;  # Input file and sequence name
		$eslsfetchCmd = $eslsfetchCmd." > /dev/null 2>&1";               # '-o' still prints some info, but we don't care.

		# For now we bail altogether if this step goes wrong (for obv. reasons)
		if (system($eslsfetchCmd)) {
		    die "\n\tERROR: Command '$eslsfetchCmd' failed (AE) - $gene\n\n"; 
		}
		
		# Generate a new exon object
		for ($i=0; $i<$num_gene_seqs; $i++) {
		    RunFastDiagonals(\$Chromosomes[$i],$ProtFileNames[$i],$nuclfilename,
				     $ProtSeqLengths[$i],$start,$end,$revcomp,
				     $gene_index_list,$timing,\@TimingData,
				     $debug);
		}
		
	    } else { # We're just going to go ahead and skip FastDiagonals and locate the plausible coding range

		for ($i=0; $i<$num_gene_seqs; $i++) {
		    
		    my $Chromosome = $Chromosomes[$i];
			
		    my $minNucl = 100000000000; # 100 GB genome? unlikely
		    my $maxNucl = -1;
		    my $some_index = 0;
		    foreach my $start_end (split(/\,/,$GeneIndex{$gene.'|'.$chromosomeName})) {
			$start_end =~ /(\d+)\-(\d+)/;
			my $start = $1;
			my $end   = $2;
			$minNucl  = MIN($start,$minNucl);
			$maxNucl  = MAX($end,$maxNucl);
			$some_index = 1;
		    }
		    if ($some_index) {
			$Chromosome->{NumHits}          = 1;
			${$Chromosome->{NuclStarts}}[0] = $minNucl;
			${$Chromosome->{NuclEnds}}[0]   = $maxNucl;
		    }
		}

	    }
	    
	    
	    # If we have hits for this "Chromosome" (really, protein/chromosome-segment pair)
	    # then we try to stitch them together in a cool and fun way.
	    for ($i=0; $i<$num_gene_seqs; $i++) {
		
		my $Chromosome = $Chromosomes[$i];
		if ($Chromosome->{NumHits}) {
		
		    if (!$noFD) {
			
			# Sort the chromosome's exons by their start positions
			$Chromosome->SortHitsByField('ProtStarts',1); 
			
			# Search for any full diagonals (that cover entire protein), if it makes sense.
			if (@{$Chromosome->{ProtStarts}} && ${$Chromosome->{ProtStarts}}[0] == 0 && !$spalner) {
			    
			    ($next_score,$next_string) 
				= FindFullAlignments(\$Chromosome,$ProtSequences[$i],
						     $ProtSeqNames[$i],$HitFile,$debug);
			    
			    # GOT FISH!
			    if ($next_score > $BestFastScores[$i]) {
				$BestFastScores[$i]  = $next_score;
				$BestFastStrings[$i] = $next_string;
				$CurrentStats[$i]    = 3;
			    }
			    
			}
			
		    }

		    # This will be the end of this FD-associated stuff
		    $TimingData[1] += Time::HiRes::tv_interval($timeA) if ($timing && !$noFD);
		    
		    # If we don't have a better hit on this chromosome punt over to SPALN
		    if ($Chromosome->{NumHits} && $CurrentStats[$i] != 3 && $spaln) {
			
			$timeA = [Time::HiRes::gettimeofday()] if ($timing);
			
			($next_score,$next_string) = ExonAssistedSPALN(\$Chromosome,\%ChrLengths,$ProtSeqNames[$i],
								       $ARGV[1],$ProtFileNames[$i],$nuclfilename,
								       $HitFile,$ProtSeqLengths[$i],$SpalnLog,$spalner,
								       $timing,\@TimingData);
			
			# SPALN saves the day!!!
			if ($next_score > $BestSpalnScores[$i]) {
			    $BestSpalnScores[$i]  = $next_score;
			    $BestSpalnStrings[$i] = $next_string;
			    $CurrentStats[$i]     = 4;
			}
			
			$TimingData[4] += Time::HiRes::tv_interval($timeA) if ($timing);
		    }

		} elsif ($timing) {
		    $TimingData[1] += Time::HiRes::tv_interval($timeA);
		}
		
	    }
		
	    #### Time for another chromosome! ####
	    
	}


	# Try setting aside anything that still hasn't given us a hit for BLAT
	for ($i=0; $i<$num_gene_seqs; $i++) {

	    if ($CurrentStats[$i] == 0) {
	    
		# Do we have the technology?
		if ($spaln && $blat) {		
		    if (system("cat $ProtFileNames[$i] >> $blatfilename")) { die "\n  Concatenation of seq. failed\n\n"; }
		
		} else {
		
		    # What a terrible loss :'(
		    print $MissFile ">$ProtSeqNames[$i]\n"; # NOTE: "GN:" ARTIFACT EXCISION
		    #print $GeneMissFile "$gene\n";
		    $HitStats[5]++;
		    $CurrentStats[$i] = 5;

		}
		
	    } else {

		$numHits++;
		$HitStats[$CurrentStats[$i]]++;
		if ($CurrentStats[$i] == 3) { print $HitFile "$BestFastStrings[$i]";  }
		else                        { print $HitFile "$BestSpalnStrings[$i]"; }
		
	    }
	    
	}
	    
    }

    # (FOR CHECKING)
    # We report how this isoform went
    my ($repX,$repY);
    for ($i=0; $i<$num_gene_seqs; $i++) {
	if ($CurrentStats[$i]) { # Because we might be blat-ing some, 0 is 'skip'
	    if ($CurrentStats[$i] == 3) {
		$CurrentStats[$i] = "diagonals  " if ($verbose);
		$repX = 1;
	    } elsif ($CurrentStats[$i] == 6) {
		$CurrentStats[$i] = "close diag " if ($verbose);
		$repX = 1;
	    } elsif ($CurrentStats[$i] == 4) {
		$CurrentStats[$i] = "spaln      " if ($verbose);
		$repX = 1;
	    } elsif ($CurrentStats[$i] == 5) {
		$CurrentStats[$i] = "bummerzone " if ($verbose);
		$repX = 0;
	    }
	    if ($Repetitive[$i]) {
		$CurrentStats[$i] = $CurrentStats[$i].' [    repetitive]' if ($verbose);
		$repY = 0;
	    } else {
		$CurrentStats[$i] = $CurrentStats[$i].' [NOT repetitive]' if ($verbose);
		$repY = 1;
	    }
	    my $nameline = "  >$ProtSeqNames[$i]"; # NOTE: "GN:" ARTIFACT EXCISION
	    while (length($nameline) < 50) {
		$nameline = $nameline.' ';
	    }
	    print "$nameline $CurrentStats[$i]\n" if ($verbose);
	    $RepConf[$repX][$repY]++;
	}
    }

    
    # Record progress, check if it's time to post a status update
    if (!$verbose) {
	$num_complete++;
	if (time()-$progtime > 15) {
	    
	    if ($threadID) {
		open(my $progfile,'>',$progfilename);
		print $progfile "$num_complete\n";
		close($progfile);
		
	    } else {
		
		my $progressCmd = './ProgressTimer.pl '.$progressbase.' '.$CPUs;
		$progressCmd = $progressCmd.' '.$num_complete.' |';
		
		open(my $progfile,$progressCmd);
		my $overall_progress = readline($progfile);
		close($progfile);
		$overall_progress =~ s/\n|\r//g;
		
		$progmsg = '  Quilter.pl:  '.$overall_progress.' isoforms examined ('.lc($reqSpecies).')';
		while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
		print "$progmsg\r";
		
	    }
	    
	    $progtime = time();
	    
	}
    }
    
    #### Another isoform bites the dust!  BRING ME THE NEXT ONE! ####
    
}
close($isoformfile);


# Non-0 threads print their stats to a temporary file
if ($threadID) {
    open($resultsfile,'>',$resultsfilename) || die "\n  ERROR:  Failed to open temporary results file '$resultsfilename'\n\n";
    if ($numIsos) {	
	print $resultsfile "$RepConf[0][0],$RepConf[0][1],$RepConf[1][0],$RepConf[1][1]\n";	
	print $resultsfile "$numHits,$numIsos\n";
	print $resultsfile "$HitStats[0]\n";
	print $resultsfile "$HitStats[1]\n";
	print $resultsfile "$HitStats[2]\n";
	print $resultsfile "$HitStats[3]\n";
	print $resultsfile "$HitStats[4]\n";
	print $resultsfile "$HitStats[5]\n";
	print $resultsfile "$HitStats[6]\n";
    } else {
	print $resultsfile "0,0,0,0\n";
	print $resultsfile "0,0\n";
	print $resultsfile "0\n";
	print $resultsfile "0\n";
	print $resultsfile "0\n";
	print $resultsfile "0\n";
	print $resultsfile "0\n";
	print $resultsfile "0\n";
	print $resultsfile "0\n";
    }

    if ($timing) {
	for ($i=1; $i<7; $i++) { print $resultsfile "$TimingData[$i]\n"; }
    }

    close($resultsfile);

}


# No need for these files anymore.
system("rm $nuclfilename") if (-e $nuclfilename);
for ($i=1; $i<=$max_gene_seqs; $i++) {
    my $protfilename = $protbasename.$i.$protfile_ext;
    system("rm $protfilename") if (-e $protfilename);
}

close($MissFile);
close($HitFile);

# de-fork
if ($threadID) { exit(0); } 
while (wait() > -1) { }


# Remove progress files
for ($i = 1; $i < $CPUs; $i++) {
    $progfilename = $progressbase.$i;
    if (-e $progfilename && system("rm $progfilename")) { die "\n  Failed to remove file $progfilename\n\n"; }
}


# Thread 0 reads everyone's results in a brazen attempt to steal the
# credit!  But we know better!  Good thing everyone loves reading I/O
# comments, so the truth can get out there!
foreach my $j (1..$CPUs-1) {

    $resultsfilename = $resultsfilehead.$j.'.Quilter.out';
    open($resultsfile,'<',$resultsfilename) || die "\n  ERROR:  Failed to open results file '$resultsfilename'\n\n";
    my @linesplitter;
    $line = <$resultsfile>;
    $line =~ s/\n|\r//g;
    @linesplitter = split(',',$line);
    $RepConf[0][0] += $linesplitter[0];
    $RepConf[0][1] += $linesplitter[1];
    $RepConf[1][0] += $linesplitter[2];
    $RepConf[1][1] += $linesplitter[3];
    $line = <$resultsfile>;
    $line =~ s/\n|\r//g;
    @linesplitter = split(',',$line);
    $numHits += $linesplitter[0];
    $numIsos += $linesplitter[1];
    foreach $i (0..6) {
	$line = <$resultsfile>;
	$line =~ s/\n|\r//g;
	$HitStats[$i] += int($line);
    }
    if ($timing) {
	for ($i=1; $i<7; $i++) {
	    $line = <$resultsfile>;
	    $line =~ s/\n|\r//g;
	    $TimingData[$i] += (0.0 + $line);
	}
    }
    close($resultsfile);
    system("rm $resultsfilename");
}


# Try to find a suitable name for the final output
# (we don't want to overwrite past outputs).
my $finalHits   = 'Hits.Quilter.out';
my $finalMisses = 'Misses.Quilter.out';
if ($resdir) {
    if ($ARGV[$resdir] =~ /\/$/) {
	$finalHits   = $ARGV[$resdir].$finalHits;
	$finalMisses = $ARGV[$resdir].$finalMisses;
    } else {
	$finalHits   = $ARGV[$resdir].'/'.$finalHits;
	$finalMisses = $ARGV[$resdir].'/'.$finalMisses;
    }
}
if (!$overw) {
    $i = 2;
    while (-e $finalHits || -e $finalMisses) {
	$finalHits   = 'Hits.'.$i.'.Quilter.out';
	$finalMisses = 'Misses.'.$i.'.Quilter.out';
	$i++;
    }
} else {
    if (-e $finalHits)   { system("rm $finalHits");   }
    if (-e $finalMisses) { system("rm $finalMisses"); }
}


# Try running BLAT on anything that we missed
if ($spaln && $blat) {

    # Start the big BLAT+SPALN timer
    $timeB = [Time::HiRes::gettimeofday()] if ($timing);

    # A mapping from numbers to sequence names
    my %BlatNameIndex;
    my %SeqLengths;
    my $BlatID = 1;

    # Grab each thread's misses (if there were any)
    my $bigBlat  = 'Quilter.BLATseqs.fa';
    my $tempProt = 'Quilter.prot.tmp.fa';
    system("rm $bigBlat") if (-e $bigBlat);

    # For each of the processes, gather the list of sequences
    # that we failed to get hits for
    my $j = 0;
    while ($j < $CPUs) {

	# Construct the name of the file where we've stashed
	# all our sequences
	$blatfilename = $foldername.'blat_tempfile_'.$j.'.Quilter.fa';

	# If there is something in this file? (if not, skip)
	if (-s $blatfilename) {

	    # Open up the file
	    open(my $blatIn,'<',$blatfilename);

	    # Prime the system by reading in the first line
	    my $line = <$blatIn>;
	    while (!eof($blatIn)) {
		
		# Strip out any whitespace and keep running
		# until we've found a line that starts a new
		# sequence (this is mainly an empty-line catch)
		$line =~ s/\n|\r//g;
		while (!eof($blatIn) && $line !~ /^\>(\S+)$/) {
		    $line = <$blatIn>;
		    $line =~ s/\n|\r|\s//g;
		}

		# If we hit the file, jump ship
		last if (eof($blatIn));

		# Open up a file to store this sequence individually.
		# We do this so that we can run esl-seqstat to (a)
		# independently verify the length, and (b) make sure
		# the result makes easel happy
		open(my $ptfile,'>',$tempProt);
		$line =~ /^\>(\S+)$/;
		my $seqname = $1;
		$BlatNameIndex{$BlatID} = $seqname;
		print $ptfile ">$BlatID\n";
		die "  $line\n" if ($BlatID eq 'GN:');
		$BlatID++;

		# Read in the actual sequence content
		$line = <$blatIn>;
		$line =~ s/\n|\r//g;
		while ($line && $line !~ /^\>/) {
		    print $ptfile "$line\n";
		    $line = <$blatIn>;
		    $line =~ s/\n|\r//g if ($line);
		}
		close($ptfile);

		# Run esl-seqstat and stash the sequence length data
		open(my $eslseqstat,"esl-seqstat --amino $tempProt |") || die "\n  ERROR: esl-seqstat failed on sequence '$seqname'\n\n";

		# For debugging: If we didn't get anything out of our seqstat, print
		# sequence name.
		while (my $statline = <$eslseqstat>) {
		    if ($statline =~ /Total \# residues\:\s+(\d+)/) {
			$SeqLengths{$seqname} = $1;
			last;
		    }
		}
		close($eslseqstat);

		# Concatenate the sequence to our big BLAT-worthy database
		system("cat $tempProt >> $bigBlat")
		
	    }

	    close($blatIn);

	}

	system("rm $blatfilename") if (-e $blatfilename);

	$j++;
	
    }
    system("rm $tempProt") if (-e $tempProt);

    # Run blat on the whole shebang (if there's a shebang to
    # speak of)
    if (-s $bigBlat) {

	# Inform user that we're trying BLAT
	$progmsg = '  Quilter.pl:  Running BLAT on missed '.lc($reqSpecies).' isoforms';
	while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
	print "$progmsg\r";
	    
	# Construct and run a happy BLAT command
	my $BlatStdOut  = 'Quilter.BLAT.std.out';
	my $BlatStdErr  = 'Quilter.BLAT.std.err';
	my $BlatResults = lc($ARGV[3]).'.Quilter.BLAT.out';
	my $blatCmd = 'blat -tileSize=7 -minIdentity=90 -maxIntron=1';
	$blatCmd    = $blatCmd.' -t=dnax -q=prot -out=blast8 -minScore=80';
	$blatCmd    = $blatCmd.' 1>'.$BlatStdOut.' 2>'.$BlatStdErr;
	$blatCmd    = $blatCmd.' '.$ARGV[1].' '.$bigBlat.' '.$BlatResults;
	$timeA = [Time::HiRes::gettimeofday()] if ($timing);
	if (system($blatCmd)) { die "\n  ERROR: BLAT command '$blatCmd' failed\n\n"; }
	$TimingData[7] = Time::HiRes::tv_interval($timeA) if ($timing);

	# As long as BLAT ran successfully we aren't actually
	# interested in what it has to say.
	if (-e $BlatStdOut) { system("rm $BlatStdOut"); }
	if (-e $BlatStdErr) { system("rm $BlatStdErr"); }

	# Inform user that we're aligning with BLAT results
	$progmsg = '  Quilter.pl:  Using SPALN on BLAT output ('.lc($reqSpecies).')';
	while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
	print "$progmsg\r";
	    
	# Run SPALN on the results of the BLAT search
	open(my $Results,'>',$finalHits);
	open(my $Misses,'>',$finalMisses);
	BLATAssistedSPALN(\%ChrLengths,\%SeqLengths,\%BlatNameIndex,$ARGV[1],$ARGV[0],
			  $BlatResults,$Results,$Misses,$SpalnLog,$ARGV[3],$CPUs,
			  $spalner,$timing,\@TimingData);
	close($Results);
	close($Misses);

	# Clean up
	system("rm $BlatResults") if (-e $BlatResults);
	system("rm $bigBlat");

    }

    $TimingData[8] = Time::HiRes::tv_interval($timeB) if ($timing);

}


# Print out that we're all done
$progmsg = "  Quilter.pl:  ".lc($reqSpecies)." complete";
while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
print "$progmsg\r" if (!$verbose);


# Close the "SpalnLog" file
close($SpalnLog);


# Now thread 0 can print out its results, if that's something you're into
if ($verbose) {
    print "\n\n";
    print "+-------------------------------------------------------+\n";
    print "                   QUILTER STATISTICS                    \n";
    print "+-------------------------------------------------------+\n\n";
    print "                        repetitive    non-repetitive\n";
    print "  Confusion: non-hitting     $RepConf[0][0]          $RepConf[0][1]\n";
    print "               hitting       $RepConf[1][0]          $RepConf[1][1]\n\n";
    print "   OVERALL : $numHits / $numIsos isoforms hit\n";
    print "           : Not Indexed   : $HitStats[0]\n";
    print "           : Wrong Species : $HitStats[1]\n";
    print "           : No Diagonals  : $HitStats[2]\n";
    print "           : Fast Score    : $HitStats[3]\n";
    print "           : (close shave) : $HitStats[6]\n";
    print "           : SPALN Score   : $HitStats[4]\n";
    print "           : BLAT Helped   : $HitStats[7]\n";
    print "           : No Score      : $HitStats[5]\n\n";
}


# Report the total number of hits we had during our search (where a
# hit is a full extension through an isoform).  Note that if we're
# doing BLAT/SPALN search, we're going to have our final list of misses
# from that
foreach $i (0..$CPUs-1) {
    $hitfilename  = $foldername.'hit_tempfile_'.$i.'.Quilter.out';
    $missfilename = $foldername.'miss_tempfile_'.$i.'.Quilter.out';
    system("cat $hitfilename >> $finalHits");
    system("cat $missfilename >> $finalMisses") unless ($spaln && $blat);
    system("rm $hitfilename");
    system("rm $missfilename");
}


# Clean up and close out temporary directory
system("rm -rf $foldername");


# If we're timing, this is where we'll make a ruckus about how long everything took
if ($timing) {
    $TimingData[11] = Time::HiRes::tv_interval($TimingData[0]);
    print "\n\n";
    print "  $ARGV[3] Quilter runtime data ('*' indicates sum over $CPUs processes):\n";
    print "\n";
    print "      TOTAL RUNTIME (wallclock) : $TimingData[11] seconds\n";
    print "\n";
    print "      Time to build index...... : $TimingData[12] seconds\n";
    print "\n";
    print "    * Time associated with FastDiagonals..... : $TimingData[1] seconds\n";
    print "    * System time spent running FastDiagonals : $TimingData[2] seconds\n";
    print "    * Number of FastDiagonals system calls... : $TimingData[3]\n";
    print "\n";
    print "    * Time associated with SPALN (non-BLAT)..... : $TimingData[4] seconds\n";
    print "    * System time spent running SPALN (non-BLAT) : $TimingData[5] seconds\n";
    print "    * Number of SPALN (non-BLAT) system calls... : $TimingData[6]\n";
    print "\n";
    print "      Time associated with BLAT+SPALN........ : $TimingData[8] seconds\n";
    print "      System time spent running BLAT......... : $TimingData[7] seconds\n";
    print "      System time spent running SPALN (+BLAT) : $TimingData[9] seconds\n";
    print "      Number of SPALN (+BLAT) system calls... : $TimingData[10]\n";
    print "\n\n";
}


1; # Terrific work, team!  You guys rock!






#########################
####                 ####
####   SUBROUTINES   ####
####                 ####
#########################





###############################################
#
# Return the maximum of two values
#
sub MAX 
{
    my $a = shift;
    my $b = shift;
    return $a if ($a > $b);
    return $b;
}



###############################################
#
# Return the minimum of two values
#
sub MIN
{
    my $a = shift;
    my $b = shift;
    return $a if ($a < $b);
    return $b;
}






#########################################################################
#
#  Function Name: TranslateCodon
#
#  About: Convert a DNA triple to an amino acid. ('x' for 'stop')
#
sub TranslateCodon
{
    
    my $codonref   = shift;
    my @codonArray = @{$codonref};
    
    # In case we didn't get a hold of a full codon (fair game) or
    # something weird slipped in (less fair, but maybe fair),
    # spit an 'X'
    return 'X' if (!(@codonArray && $codonArray[0] && $codonArray[1] && $codonArray[2]));
    my $codon = $codonArray[0].$codonArray[1].$codonArray[2];
    $codon    =  uc($codon);
    $codon    =~ s/U/T/g;
    return 'X' if ($codon =~ /[^ACGT]/);

    
    # Guarantees uppercase
    @codonArray = split('',$codon);

    
    if ($codonArray[0] eq 'A') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "AAA") { return 'K'; }
	    if ($codon eq "AAC") { return 'N'; }
	    if ($codon eq "AAG") { return 'K'; }
	    if ($codon eq "AAT") { return 'N'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "ACA") { return 'T'; }
	    if ($codon eq "ACC") { return 'T'; }
	    if ($codon eq "ACG") { return 'T'; }
	    if ($codon eq "ACT") { return 'T'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "AGA") { return 'R'; }
	    if ($codon eq "AGC") { return 'S'; }
	    if ($codon eq "AGG") { return 'R'; }
	    if ($codon eq "AGT") { return 'S'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "ATA") { return 'I'; }
	    if ($codon eq "ATC") { return 'I'; }
	    if ($codon eq "ATG") { return 'M'; }
	    if ($codon eq "ATT") { return 'I'; }
	}
    }

    if ($codonArray[0] eq 'C') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "CAA") { return 'Q'; }
	    if ($codon eq "CAC") { return 'H'; }
	    if ($codon eq "CAG") { return 'Q'; }
	    if ($codon eq "CAT") { return 'H'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "CCA") { return 'P'; }
	    if ($codon eq "CCC") { return 'P'; }
	    if ($codon eq "CCG") { return 'P'; }
	    if ($codon eq "CCT") { return 'P'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "CGA") { return 'R'; }
	    if ($codon eq "CGC") { return 'R'; }
	    if ($codon eq "CGG") { return 'R'; }
	    if ($codon eq "CGT") { return 'R'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "CTA") { return 'L'; }
	    if ($codon eq "CTC") { return 'L'; }
	    if ($codon eq "CTG") { return 'L'; }
	    if ($codon eq "CTT") { return 'L'; }
	}
    }
    
    if ($codonArray[0] eq 'G') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "GAA") { return 'E'; }
	    if ($codon eq "GAC") { return 'D'; }
	    if ($codon eq "GAG") { return 'E'; }
	    if ($codon eq "GAT") { return 'D'; }
	}
	if ($codonArray[1] eq 'C') {	    
	    if ($codon eq "GCA") { return 'A'; }
	    if ($codon eq "GCC") { return 'A'; }
	    if ($codon eq "GCG") { return 'A'; }
	    if ($codon eq "GCT") { return 'A'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "GGA") { return 'G'; }
	    if ($codon eq "GGC") { return 'G'; }
	    if ($codon eq "GGG") { return 'G'; }
	    if ($codon eq "GGT") { return 'G'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "GTA") { return 'V'; }
	    if ($codon eq "GTC") { return 'V'; }
	    if ($codon eq "GTG") { return 'V'; }
	    if ($codon eq "GTT") { return 'V'; }
	}
    }
    if ($codonArray[0] eq 'T') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "TAA") { return 'x'; }
	    if ($codon eq "TAC") { return 'Y'; }
	    if ($codon eq "TAG") { return 'x'; }
	    if ($codon eq "TAT") { return 'Y'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "TCA") { return 'S'; }
	    if ($codon eq "TCC") { return 'S'; }
	    if ($codon eq "TCG") { return 'S'; }
	    if ($codon eq "TCT") { return 'S'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "TGA") { return 'x'; }
	    if ($codon eq "TGC") { return 'C'; }
	    if ($codon eq "TGG") { return 'W'; }
	    if ($codon eq "TGT") { return 'C'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "TTA") { return 'L'; }
	    if ($codon eq "TTC") { return 'F'; }
	    if ($codon eq "TTG") { return 'L'; }
	    if ($codon eq "TTT") { return 'F'; }
	}
    }

    # Weird codon is weird. TO THE BIN WITH YOU!
    return 'X';
    
}





#########################################################################
#
#  Function Name:  GetStopPoints
#
#  About: Yo, you wanna know where your threads are supposed to focus
#         while they're mapping proteins to the genome?  Yeah you do,
#         so check this function out for reals.  We scan the protein
#         database and try to break things down in a way that's divides
#         work evenly across threads.
#
sub GetStopPoints
{
    my $db_name  = shift;
    my $req_spec = shift; # required species
    my $num_CPUs = shift;

    # We actually do a preliminary scan to make sure that there are a
    # sufficient number of sequences to warrant the requested number
    # of processes.
    #
    # Note that we assume that the database is organized by family,
    # although use of '-n' (rather than '-nplus') indicates that this
    # isn't a Mirage usage...
    #
    open(my $DB,'<',$db_name) || die "\n  ERROR:  Failed to open protein database '$db_name' (GSP1)\n\n";
    my $num_fams = 0;
    my $last_fam = 0;
    while (my $line = <$DB>) {
	if ($line =~ /^\>(\S+)/) {
	    my $seqname = $1;
	    $seqname =~ /^[^\|]+\|[^\|]+\|([^\|]+)\|/;
	    my $species = lc($1);
	    if ($species eq $req_spec) {
		$seqname =~ /\|([^\|]+)$/;
		my $next_fam = lc($1);
		if ($last_fam && $last_fam ne $next_fam) {
		    $num_fams++;
		}
		$last_fam = $next_fam;
	    }
	}
    }
    close($DB);

    my @StopPoints;
    if ($num_fams == 0) { ################################# No entries -- I-like-a-do-WHAAAAT?!

	die "\n  ERROR:  Species '$req_spec' doesn't have any entries in '$db_name'\n\n";

    } elsif ($num_fams <= $num_CPUs) { #################### Need to just assign one fam per CPU

	$num_CPUs = $num_fams;
	open(my $DB,'<',$db_name) || die "\n  ERROR:  Failed to open protein database '$db_name' (GSP2)\n\n";
	$last_fam = 0;
	my $line_num = 0;
	while (my $line = <$DB>) {

	    $line_num++;

	    if ($line =~ /^\>(\S+)/) {
		
		my $seqname = $1;
		$seqname =~ /^[^\|]+\|[^\|]+\|([^\|]+)\|/;
		my $species = lc($1);

		if ($species eq $req_spec) {
		    $seqname =~ /\|([^\|]+)$/;
		    my $next_fam = lc($1);
		    if ($last_fam && $last_fam ne $next_fam) { push(@StopPoints,$line_num); }
		    $last_fam = $next_fam;
		}

	    }
	}
	close($DB);
	
    } else { ############################################## Ripped from Mirage.pl

	my $wc_cmd = "wc -l \"$db_name\" \|";
	open(my $WC,$wc_cmd) || die "\n  ERROR:  Failed to run linecount command '$wc_cmd'\n\n";
	my $wc_line = <$WC>;
	close($WC);

	my $tot_line_count;
	if ($wc_line =~ /^\s*(\d+)/) {
            $tot_line_count = $1;
        } else {
            die "\n  ERROR:  Failed to get a line count for file '$db_name'\n\n";
        }

	open(my $DB,'<',$db_name) || die "\n  ERROR:  Failed to open protein database '$db_name' (GSP3)\n\n";

	my $avg_frac = $tot_line_count / $num_CPUs;
	
	my $fams_in_block = 0;

	my $lastfam;
	my $linenum = 0;
	my $j=0;
	while ($j<$num_CPUs) {

	    my $line;
	    
	    while ($linenum < $tot_line_count && $linenum < ($j+1)*$avg_frac) {
		$line = <$DB>;
		$linenum++;
		if ($line =~ /^\>\S+\|([^\|]+)$/) {
		    $lastfam = lc($1);
		    $fams_in_block = 1;
		}
	    }
	    
	    # This is a catch for the special case where we might have
	    # (say) numProcesses+1 sequences where one of those is gigantic
	    # and goofs up our best laid plans.
	    #
	    if ($fams_in_block) {
		
		while ($linenum < $tot_line_count && $line !~ /^\>/) {
		    $line = <$DB>;
		    $linenum++;
		    if ($line =~ /^\>\S+\|([^\|]+)$/) {
			my $nextfam = lc($1);
			if ($nextfam eq $lastfam) {
			    $line = <$DB>;
			    $linenum++;
			}
		    }
		}
		
		push(@StopPoints,$linenum);
		$j++;
		
	    } else {

		$num_CPUs--;

	    }

	}
	close($DB);
	
    }

    return($num_CPUs,\@StopPoints);

}





#########################################################################
#
#  Function Name: BuildGeneIndex
#
#  About: BuildGeneIndex generates a hash the known locations of genes
#         in a species' genome, as provided in the .gtf index file.
#         This hash (rather, hash of hashes) has the gene name as its
#         primary key and the chromosome as its secondary key.
#         The values stored in the hash are arrays corresponding to pairs
#         of start and end positions of exons known to code for the
#         gene (these positions being nucleotide positions in the chromosome).
#         The function also tracks which genes have no entries in the
#         index, storing these in a simple 'Blacklist' hash.
#         The gene index hash and blacklist hash are both returned at the
#         end of the function.
#
sub BuildGeneIndex
{
    my $dbname  = shift;
    my $gtfname = shift;

    my $species = shift;

    my $debug = shift;

    my $verbose = shift;

    print "\n" if ($debug);
    
    my %FamNames;  # Final Entry
    my %GeneNames; # First Entry (minus "GN:") -- a list of synonyms

    open(my $db,'<',$dbname) || die "\n\tCouldn't open '$dbname'\n\n";
    while (my $line = <$db>) {

	$line =~ s/\n|\r|\s//g;
	next if (!$line || $line !~ /\>/);

	$line =~ /^\>([^\|]+)\|[^\|]+\|([^\|]+)\|\S+\|([^\|]+)$/;
	my $fam  = uc($3);
	my $spec = uc($2);
	my $gene = uc($1);

	if ($spec eq $species) {

	    $FamNames{$fam} = $fam;

	    # NOTE: We can't totally get away from the 'GN:' thing in our UniProt
	    #       database, so we need to preserve a teensy little fragment
	    #       of our prior excessive 'GN:'-consciousness
	    if ($gene ne '-') {
		$gene =~ s/^GN\://;
		$FamNames{$gene} = $fam;
	    }

	}
	
    }
    close($db);

    my %IndexedFams;
    
    my %FamChrs;
    my %GeneIndex;
    my %RepeatTracker;
    open(my $gtf,'<',$gtfname) || die "\n\tCouldn't open '$gtfname'\n\n";
    while (my $line = <$gtf>) {

	next if ($line =~ /^\#/ || lc($line) !~ /\S+\s+\S+\s+exon/); # This is dumb, Alex.  You know why (wild boy).
	$line =~ /^(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+([\+\-])/;
	
	# Note that startPos is always less than endPos, even when we're
	# indexing into the reverse complement.
	my $chrName  = $1;
	my $startPos = $2;
	my $endPos   = $3;
	my $revcomp  = 0;
	if ($4 eq '-') {
	    $revcomp = 1;
	    $chrName = $chrName.'[revcomp]';
	    my $tempPos = $startPos;
	    $startPos = $endPos;
	    $endPos = $tempPos;
	}

	$line =~ /gene\_name \"([^\"]+)\"\;/;
	my $geneName = uc($1);
	my $famName  = $FamNames{$geneName};
	
	next if (!($chrName && $startPos && $endPos && $famName) || abs($endPos-$startPos)<4);

	my $repeat_id = $chrName.'|'.$startPos.'|'.$endPos.'|'.$famName;
	if (!$RepeatTracker{$repeat_id} && $FamNames{$famName}) {

	    if ($FamChrs{$famName}) {
		my $already_there = 0;
		for my $chrs_for_fam (split(/\#/,$FamChrs{$famName})) {
		    if ($chrs_for_fam eq $chrName) {
			$already_there = 1;
		    }
		}
		if ($already_there == 0) {
		    $FamChrs{$famName} = $FamChrs{$famName}.'#'.$chrName;
		}
	    } else {
		$FamChrs{$famName} = $chrName;
	    }

	    my $key = $famName.'|'.$chrName;
	    if ($GeneIndex{$key}) { $GeneIndex{$key} = $GeneIndex{$key}.','.$startPos.'-'.$endPos; }
	    else                  { $GeneIndex{$key} = $startPos.'-'.$endPos;                      }
	    $RepeatTracker{$repeat_id} = 1;
	    $IndexedFams{$famName}     = 1;
	    
	}
	
    }
    close($gtf);

    # We'll want to have the exons sorted according to nucleotide positions (ascending if
    # forward strand, descending if reverse complement).
    foreach my $fam_chr_pair (keys %GeneIndex) {

	$fam_chr_pair =~ /\S+\|(\S+)/;
	my $chr = $1;
	my $revcomp = 0;
	$revcomp = 1 if ($chr =~ /\[revcomp\]/);

	# This is a fun algorithm called "lazy sort"
	my %StartPosHash;
	foreach my $start_end_pair (split(/\,/,$GeneIndex{$fam_chr_pair})) {
	    $start_end_pair =~ /(\d+)\-(\d+)/;
	    my $start = $1;
	    my $end   = $2;
	    if ($StartPosHash{$start}) { $StartPosHash{$start} = $StartPosHash{$start}.','.$end; }
	    else                       { $StartPosHash{$start} = $end;                           }	    
	}

	# Keep that lazy train a-rolling
	my $sorted_entry = '';
	if ($revcomp) {

	    foreach my $start (sort {$b <=> $a} keys %StartPosHash) {
	        foreach my $end (sort {$b <=> $a} split(/\,/,$StartPosHash{$start})) {
		    $sorted_entry = $sorted_entry.','.$start.'-'.$end;
		}
	    }
	    $sorted_entry =~ s/^\,//;
	    
	} else {
	    
	    foreach my $start (sort {$a <=> $b} keys %StartPosHash) {
	        foreach my $end (sort {$a <=> $b} split(/\,/,$StartPosHash{$start})) {
		    $sorted_entry = $sorted_entry.','.$start.'-'.$end;
		}
	    }
	    $sorted_entry =~ s/^\,//;
	    
	}

	# Wooot! What a nice, lazy sort!
	$GeneIndex{$fam_chr_pair} = $sorted_entry;

    }

    # Sanity check
    my %MissingFams;
    print "Names of missing gene families:\n\n" if ($debug);
    foreach my $fam (sort keys %FamNames) {
	if (!$IndexedFams{$fam}) {
	    print "  $fam\n" if ($debug);
	    $MissingFams{$fam} = 1;
	}
    }
    print "\n" if ($debug);

    return (\%GeneIndex,\%FamChrs,\%MissingFams);
    
}





#########################################################################
#
#  Function Name: GetChromsomeLengths
#
#  About: Because we sometimes search outside of the strict range
#         provided by tools like blat and spaln we need to be
#         safe about where our end points are.
#
sub GetChromosomeLengths
{
    my $nuclfilename = shift;
    my $seqstatCmd = 'esl-seqstat --dna -a '.$nuclfilename.' |';

    my %ChrLengths;    
    open(my $seqstats,$seqstatCmd);
    while (my $stat = <$seqstats>) {
	if ($stat =~ /\= (\S+)\s+(\d+)/) {     
	    $ChrLengths{$1} = $2;
	}
    }
    close($seqstats);

    return \%ChrLengths;
}







#########################################################################
#
#  Function Name: RunFastDiagonals
#
#  About: RunFastDiagonals performs the system call to search the protein
#         against some section of the genome (the particulars of this 'section'
#         corresponds to some entry in the .gtf index).  It then uses the
#         output of the program to add to our current archive of hits for
#         this particular gene and chromosome.         
#
sub RunFastDiagonals
{

    # Get reference to the 'DiagonalSet' object and lock it down
    my $chrRef     = shift;
    my $Chromosome = $$chrRef;

    # Get the names of the protein and nucleotide files
    my $protfile = shift;
    my $nuclfile = shift;

    # How long is our protein?
    my $proteinLength = shift;

    # Get the start and end points of the chromosome
    my $RegionStart = shift;
    my $RegionEnd   = shift;

    # Is this search on the chromosome's reverse complement?
    my $revcomp = shift;

    # The full list of gene index entries for this gene family
    my $gene_index_list = shift;

    # Are we timing?
    my $timing     = shift;
    my $timingdata = shift;
    my $timeA;

    # Are we debugging?
    my $debug = shift;

    # Record how many diagonals each protein sequence has
    my $diagNum = 0;

    # PRINT THAT SUPER-FLY FastDiagonals.c DEBUGGING OUTPUT
    # If we want to peek at how FastDiagonals is doing, we can run it in debug mode,
    # SEPARATELY from running it in parsable-format
    #if ($debug && system("\.\/FastDiagonals $protfile $nuclfile \-debug")) { die "\n\tFastDiagonals failed\n\n"; }

    # Break our list of gene index entries into a proper list
    my @TrueIndexList  = split(/\,/,$gene_index_list);
    my $gene_index_len = scalar(@TrueIndexList);
    my @RelativeList;
    foreach my $start_end_pair (@TrueIndexList) {

	$start_end_pair =~ /(\d+)\-(\d+)/;
	my $start = $1;
	my $end   = $2;

	if ($revcomp) {
	    $start = $RegionStart - $start;
	    $end   = $RegionStart - $end;
	} else {
	    $start = $start - $RegionStart;
	    $end   = $end   - $RegionStart;
	}

	push(@RelativeList,$start);
	push(@RelativeList,$end);

    }
    
    my $nucl_seq_len = abs($RegionStart-$RegionEnd)+1;

    # Construct our FastDiagonals command
    my $FDcmd = './FastDiagonals "'.$protfile.'" "'.$nuclfile.'" ';
    $FDcmd = $FDcmd.' '.$nucl_seq_len.' '.$gene_index_len;
    foreach my $relative_pos (@RelativeList) { $FDcmd = $FDcmd.' '.$relative_pos; }
    $FDcmd = $FDcmd.' |';

    # Run FastDiagonals
    $timeA = [Time::HiRes::gettimeofday()] if ($timing);
    open(my $stdoutput,$FDcmd) || die "\n\tFastDiagonals failed (command: FDcmd)\n\n";
    if ($timing) {
	@{$timingdata}[2] += Time::HiRes::tv_interval($timeA);
	@{$timingdata}[3]++;
    }

    # We should be able to just iterate over our IndexList entries,
    # thanks to the way that FastDiagonals provides output...
    foreach my $index_entry (@TrueIndexList) {

	$index_entry  =~ /(\d+)\-(\d+)/;
	my $NuclStart = $1;
	my $NuclEnd   = $2;

	my $startoffsetChars;
	my $numDiagonals;
	my $fullDiagLength;
	my @DiagProtStarts;
	my @DiagProtEnds;
	my @DiagScores;
	my $numTerms;
	my @TermProtStarts;
	my @TermNuclEnds;
	my @TermScores;
	my $stopCodonFound;
	my $endoffsetLength;
	my $endoffsetChars;
	my $line;
	my @genTuple; # Generic tuple
	
	# Considering offsets of 0, 1, and 2
	my $startoffsetLength = 0;
	while ($startoffsetLength < 2) { # NOTE: Incrementing happens at the end, so 2 DOES get counted
	    
	    # Resetting our arrays
	    @DiagProtStarts   = ();
	    @DiagProtEnds     = ();
	    @DiagScores       = ();
	    @TermProtStarts   = ();
	    @TermNuclEnds     = ();
	    @TermScores       = ();
	    
	    # 1.
	    # Read in the current starting offset and record the
	    # characters this corresponds to (if relevant)
	    # Try to get ahead of failures:
	    $startoffsetLength = <$stdoutput>;
	    while ($startoffsetLength !~ /\S/) { $startoffsetLength = <$stdoutput>; }
	    if ($startoffsetLength =~ /ERROR/) { die "\n  $startoffsetLength (FD on $protfile and $nuclfile)\n\n"; }
	    last if ($startoffsetLength !~ /[012]/);
	    
	    $startoffsetLength = int($startoffsetLength);
	    
	    if ($startoffsetLength) {
		$line = <$stdoutput>;
		$line =~ s/\n|\r//g;
		$startoffsetChars = $line;
	    } else {
		$startoffsetChars = 0;
	    }
	    
	    # 2.
	    # How long is the end offset? This is also where
	    # we check whether a stop codon was hit, which
	    # is recorded.  Eventually we might want to just not
	    # record non-terminal hits that end in a stop codon...
	    $endoffsetLength = <$stdoutput>;
	    $endoffsetLength =~ s/\r|\n//g;
	    if ($endoffsetLength eq 'X') { # Stop codon!
		
		$stopCodonFound = 1;
		
	    } else {
		
		$stopCodonFound = 0;
		
		$endoffsetLength = int($endoffsetLength);
		if ($endoffsetLength) {
		    $line = <$stdoutput>;
		    $line =~ s/\n|\r//g;
		    $endoffsetChars = $line;
		} else {
		    $endoffsetChars = 0;
		}
	    }
	    
	    
	    # 3.
	    # Read in how many diagonals there are, and how many
	    # amino acids a full diagonal contains.  Then, read
	    # in the diagonal starting positions and scores.
	    $line = <$stdoutput>;
	    $line =~ s/\r|\n//g;
	    @genTuple = split(' ',$line);
	    
	    $numDiagonals   = int($genTuple[0]);
	    $fullDiagLength = int($genTuple[1]);
	    
	    my $codonStart = $NuclStart;
	    if ($revcomp) { $codonStart -= $startoffsetLength; }      
	    else          { $codonStart += $startoffsetLength; }
	    
	    my $codonEnd = $NuclEnd; # In the event of a stop codon, go long
	    if (!$stopCodonFound) {
		if ($revcomp) { $codonEnd = ($codonStart-((3*$fullDiagLength)-1))-$endoffsetLength; }
		else          { $codonEnd = ($codonStart+((3*$fullDiagLength)-1))+$endoffsetLength; }
	    }
	    
	    foreach my $i (0..$numDiagonals-1) {
		
		$line = <$stdoutput>;
		$line =~ s/\r|\n//g;
		
		# We don't need to record non-terminal hits that ended in
		# a stop codon.
		if (!$stopCodonFound) {
		    
		    @genTuple = split(' ',$line);
		    
		    ${$Chromosome->{ProtStarts}}[$diagNum]  = $genTuple[0];
		    ${$Chromosome->{ProtEnds}}[$diagNum]    = $genTuple[0]+$fullDiagLength-1;
		    ${$Chromosome->{ProtStrings}}[$diagNum] = ${$Chromosome->{ProtStarts}}[$diagNum].'-'.${$Chromosome->{ProtEnds}}[$diagNum];
		    ${$Chromosome->{NuclStarts}}[$diagNum]  = $NuclStart;
		    ${$Chromosome->{NuclEnds}}[$diagNum]    = $NuclEnd;
		    
		    ${$Chromosome->{StartOffsetLengths}}[$diagNum] = $startoffsetLength;
		    ${$Chromosome->{StartOffsets}}[$diagNum]       = $startoffsetChars;
		    ${$Chromosome->{EndOffsetLengths}}[$diagNum]   = $endoffsetLength;
		    ${$Chromosome->{EndOffsets}}[$diagNum]         = $endoffsetChars;
		    
		    ${$Chromosome->{NuclStrings}}[$diagNum] = $codonStart.'-'.$codonEnd;
		    
		    ${$Chromosome->{Scores}}[$diagNum]     = $genTuple[1];
		    ${$Chromosome->{StopCodon}}[$diagNum]  = 0;
		    ${$Chromosome->{IsTerminal}}[$diagNum] = 0;
		    
		    $diagNum++;
		}
		
	    }
	    
	    # 4.
	    # Read in how many diagonals we started and then
	    # terminated because they reached the end of the 
	    # protein sequence.
	    $numTerms = int(<$stdoutput>);
	    foreach my $i (0..$numTerms-1) {
		
		$line = <$stdoutput>;
		$line =~ s/\r|\n//g;
		@genTuple = split(' ',$line);
		
		# Add all our information into the Chromosome (DiagonalSet object)
		${$Chromosome->{ProtStarts}}[$diagNum]  = $genTuple[0];
		${$Chromosome->{ProtEnds}}[$diagNum]    = $proteinLength-1;
		${$Chromosome->{ProtStrings}}[$diagNum] = ${$Chromosome->{ProtStarts}}[$diagNum].'-'.${$Chromosome->{ProtEnds}}[$diagNum];
		
		# There's no actual need to know the end nucleotide on terminals -- just walk it in!
		${$Chromosome->{NuclStarts}}[$diagNum] = $NuclStart;
		if ($revcomp) {
		    ${$Chromosome->{NuclEnds}}[$diagNum] = $NuclStart-(3*(1+$proteinLength-$genTuple[0])-1);
		} else {
		    ${$Chromosome->{NuclEnds}}[$diagNum] = $NuclStart+(3*(1+$proteinLength-$genTuple[0])-1);
		}
		
		${$Chromosome->{StartOffsetLengths}}[$diagNum] = $startoffsetLength;
		${$Chromosome->{StartOffsets}}[$diagNum]       = $startoffsetChars;
		${$Chromosome->{EndOffsetLengths}}[$diagNum]   = 0;
		${$Chromosome->{EndOffsets}}[$diagNum]         = 0;
		
		${$Chromosome->{NuclStrings}}[$diagNum] = $codonStart.'-X';
		
		${$Chromosome->{Scores}}[$diagNum]     = $genTuple[1];
		${$Chromosome->{StopCodon}}[$diagNum]  = $stopCodonFound;
		${$Chromosome->{IsTerminal}}[$diagNum] = 1;
		
		$Chromosome->{NumTerminals}++;
		$diagNum++;
	    }
	    
	    $Chromosome->{NumHits} = $diagNum;
	}

    }

    # Exon down!
    close($stdoutput);
    
}








#########################################################################
#
#  Function Name: FindFullAlignments
#
#  About: FindFullAlignments is the top-level subroutine for stitching
#         together diagonals (runs of near-identity between protein and
#         genome segments), trying to find some set of diagonals that
#         covers the entire protein.  WARNING:  This subroutine might
#         too much for immature developers to handle (it has a lot of
#         variables starting with '$hit').
#
sub FindFullAlignments
{
    my ($i,$j,$k);
    
    my $chromosomeRef = shift;
    my $Chromosome = $$chromosomeRef;
    
    my $Protein    = shift;
    my $protlength = length($Protein);

    my $seqname = shift;

    my $HitFile = shift;

    my $debug = shift;
    my $happy = 0;
    
    # We're only concerned, at this point, about hits that start at position
    # zero of the protein.
    my $Hit = 0;

    # There are rare cases where HUGE numbers of hits with identical start and
    # end points (in the protein sequence) cause us to perform absurb amounts
    # of tree traversal, so we can use this to track red herrings.
    my %RedHerrings;
    my %AlreadyRecd;
    
    # Are we working on the reverse complement?
    my $revcomp;
    if ($Chromosome->{ChrName} =~ /\[revcomp\]/) { $revcomp = 1; }
    else                                         { $revcomp = 0; }

    # Let's nab us some alignments!
    my @FinalNuclStrs;
    my @FinalProtStrs;
    my @FinalScores;
    while ($Hit < $Chromosome->{NumHits} && ${$Chromosome->{ProtStarts}}[$Hit] == 0) {
	
	# We have (at least the start of) a hit! Check if this is a full match.
	# Otherwise, walk down the exons and try to get to a terminal.
	my $hitLocated = ${$Chromosome->{IsTerminal}}[$Hit];

	# Looks like we got our work cut out for us...
	if (!$hitLocated) {

	    # Lock and load!
	    push(my @hits,$Hit);
	    push(my @nuclStrings,${$Chromosome->{NuclStrings}}[$Hit]);
	    push(my @protStrings,${$Chromosome->{ProtStrings}}[$Hit]);
	    push(my @cumulativeScores,${$Chromosome->{Scores}}[$Hit]);

	    # This looks a little ugly, but gets us out of deep recursion.
	    while (@hits) {
		
		# Pop the last hit position (and its full string) off the stack
		my $prevHit         = pop(@hits);
		my $prevNuclString  = pop(@nuclStrings);
		my $prevProtString  = pop(@protStrings);
		my $prevScore       = pop(@cumulativeScores);
		my $prevProtStart   = ${$Chromosome->{ProtStarts}}[$prevHit];
		my $prevProtEnd     = ${$Chromosome->{ProtEnds}}[$prevHit];
		my $prevStartOffset = ${$Chromosome->{StartOffsetLengths}}[$prevHit];
		my $prevEndOffset   = ${$Chromosome->{EndOffsetLengths}}[$prevHit];

		# Next hit position must be:
		my $nextHit       = $prevHit + 1;
		my $nextProtStart = $prevProtEnd + 1;
		my $nextOffset    = (3 - $prevEndOffset) % 3;
		$nextProtStart++ if ($nextOffset);
		
		# We assume that we're going to be disappointed, but might get
		# proved wrong later down the line.
		if (!$RedHerrings{$prevStartOffset.','.$prevProtStart.','.$prevProtEnd.','.$prevEndOffset}){
		    if ($prevProtStart) {
			$RedHerrings{$prevStartOffset.','.$prevProtStart.','.$prevProtEnd.','.$prevEndOffset} = 1
		    }
		} else {
		    next;
		}
		
		if ($nextProtStart >= $protlength && !$AlreadyRecd{$prevNuclString.$prevProtString}) {

		    # As a special case, if we only have one more amino acid to go
		    # and we're part-way in (i.e., nextOffset > 0) we goof a tiny
		    # bit.
		    $prevProtString =~ /\-(\d+)$/;
		    my $prot_final_end = int($1);
		    if ($prot_final_end != $protlength-1) {

			$AlreadyRecd{$prevNuclString.$prevProtString} = 1;

			$prot_final_end =  $protlength-1;
			$prevProtString =~ s/\-\d+$/\-$prot_final_end/;

			$prevNuclString =~ /\-(\d+)/;
			my $nucl_final_end = int($1);
			if ($revcomp) { $nucl_final_end -= 3; }
			else          { $nucl_final_end += 3; }
			$prevNuclString =~ s/\-\d+$/\-$nucl_final_end/;

		    }
		    
		    push(@FinalNuclStrs,$prevNuclString);
		    push(@FinalProtStrs,$prevProtString);
		    push(@FinalScores,$prevScore);
		    $AlreadyRecd{$prevNuclString.$prevProtString} = 1;
		    next;
		}
		
		# Walk along until we find a suitable starting point -- eventually convert to binary search
		while ($nextHit < $Chromosome->{NumHits} && ${$Chromosome->{ProtStarts}}[$nextHit] < $nextProtStart) {
		    $nextHit++;
		}
		
		# Check if we can get a slick extension going!
		# If we can't, blacklist this dang position.
		while ($nextHit < $Chromosome->{NumHits} && ${$Chromosome->{ProtStarts}}[$nextHit] == $nextProtStart) {
		    
		    # Making sure this candidate makes sense with what we've seen
		    if (($revcomp && ${$Chromosome->{NuclStarts}}[$nextHit] >= ${$Chromosome->{NuclEnds}}[$prevHit])
			|| (!$revcomp && ${$Chromosome->{NuclStarts}}[$nextHit] <= ${$Chromosome->{NuclEnds}}[$prevHit])
			|| ${$Chromosome->{StartOffsetLengths}}[$nextHit] != $nextOffset) 
		    {
			$nextHit++;
			next;
		    }
		    
		    # Note that we don't put the end position on, yet, since we
		    # don't know whether we have an offset to worry about or not
		    my $nextNuclString = $prevNuclString.','.${$Chromosome->{NuclStrings}}[$nextHit];
		    my $nextProtString = $prevProtString.','.${$Chromosome->{ProtStrings}}[$nextHit];			
		    
		    # Try to achieve tranquility:
		    if (${$Chromosome->{IsTerminal}}[$nextHit]) {
			
			# Tranquiltiy achieved!
			$hitLocated = 1;
			my $finalScore = $prevScore + ${$Chromosome->{Scores}}[$nextHit]; # Picking up 'cut-off points'...
			if (!$AlreadyRecd{$nextNuclString.$nextProtString}) {
			    push(@FinalNuclStrs,$nextNuclString);
			    push(@FinalProtStrs,$nextProtString);
			    push(@FinalScores,$finalScore);
			    $AlreadyRecd{$nextNuclString.$nextProtString} = 1;
			}
			
		    } else {
			    
			# Not, yet, but we'll keep working at it
			push(@hits,$nextHit);
			push(@nuclStrings,$nextNuclString);
			push(@protStrings,$nextProtString);
			push(@cumulativeScores,($prevScore+${$Chromosome->{Scores}}[$nextHit]));

		    }

		    $nextHit++;

		}
	    }
	    
	} else {
	    
	    # 2 EZ
	    if (!$AlreadyRecd{${$Chromosome->{NuclStrings}}[$Hit].${$Chromosome->{ProtStrings}}[$Hit]}) {
		push(@FinalNuclStrs,${$Chromosome->{NuclStrings}}[$Hit]);
		push(@FinalProtStrs,${$Chromosome->{ProtStrings}}[$Hit]);
		push(@FinalScores,${$Chromosome->{Scores}}[$Hit]);
		$AlreadyRecd{${$Chromosome->{NuclStrings}}[$Hit].${$Chromosome->{ProtStrings}}[$Hit]} = 1;
	    }
	    
	}
	
	# Did we ever get around to finding a hit?
	$happy = 1 if ($hitLocated);
	
	# Next up!
	$Hit++;
    }

    # If we didn't hit, bail here
    return (0,0) if (!$happy);
    
    # Now that we've got our great big list of hits, print them off!
    # Note that we could sort the scores and select only the (best? best n?)
    my $max_score = $FinalScores[0];
    my $max_pos   = 0;
    foreach $Hit (1..@FinalNuclStrs-1) {
	if ($FinalScores[$Hit] > $max_score) {
	    $max_pos   = $Hit;
	    $max_score = $FinalScores[$Hit];
	}
    }

    # Convert to a position list, rather than an array of ranges.
    my $FinalList  = ListCodonCenters($FinalNuclStrs[$max_pos],$FinalProtStrs[$max_pos],$revcomp);

    # Put everything together in a nice string and pass on up the chain of command!
    my $hitstring = "Isoform ID : $seqname\n";
    $hitstring    = $hitstring."Method     : FastDiagonals\n";
    $hitstring    = $hitstring."Chromosome : $Chromosome->{ChrName}\n";
    $hitstring    = $hitstring."Match Pos.s: $FinalList\n\n";
    
    # Yay!
    return ($max_score,$hitstring);
    
}





#########################################################################
#
#  Function Name: ListCodonCenters
#
#  About: This function takes a string representing a range of genome
#         indices and returns a list corresponding to the positions in
#         the given ranges that would be the centers of codons.
#         It is important that this is only correct when each triple is
#         a codon
#
sub ListCodonCenters
{
    my $rangestr = shift;
    my $protlocs = shift;
    my $revcomp  = shift;

    #print "\n$protlocs\n$rangestr\n\n";
    
    $protlocs    =~ /\-(\d+)$/;
    my $numprots = $1 + 1;
    my $obsprots = 0;
    my @ProtRanges = split(',',$protlocs);
    my @Ranges   = split(',',$rangestr);
    my $numexons = @Ranges;

    
    my $CodonList = '';

    my $pickup = 0;
    
    foreach my $i (0..$numexons-1) {
	
	my @nuclpair = split('-',$Ranges[$i]);
	my @protpair = split('-',$ProtRanges[$i]);

	# Be sure to start in the center of an exon
	my $nuclpos = $nuclpair[0];
	if ($revcomp) { $nuclpos--; }
	else          { $nuclpos++; }

	# Do we need to pick up a "filler" codon?
	if ($pickup) {
	    if ($revcomp) { $nuclpos += 3; }
	    else          { $nuclpos -= 3; }
	}
	
	# If we're on the final exon we can just "walk it in"
	if ($i == $numexons-1) {

	    while ($obsprots < $numprots) {

		$CodonList = $CodonList.','.$nuclpos;
		$obsprots++;
		
		if ($revcomp) { $nuclpos -= 3; }
		else          { $nuclpos += 3; }

	    } 
	    
	} else {

	    # Stroll through the exon
	    while (($revcomp && $nuclpos >= $nuclpair[1]) 
		   || (!$revcomp && $nuclpos <= $nuclpair[1])) {
		
		$CodonList = $CodonList.','.$nuclpos;
		$obsprots++;
		
		if ($revcomp) { $nuclpos -= 3; }
		else          { $nuclpos += 3; }
		
	    }

	    # Mark the splice site
	    $CodonList = $CodonList.',*';

	    # Calculate the next offset
	    if (($revcomp && $nuclpair[1]-$nuclpos == 1) 
		|| (!$revcomp && $nuclpos-$nuclpair[1] == 1)) {
		$pickup = 1;
	    } else {
		$pickup = 0;
	    } 

	}
	
    }

    # There will be extra commas at the start and end to shear off
    $CodonList =~ s/^\,//;
    $CodonList =~ s/\,$//;

    return $CodonList;

}





#########################################################################
#
#  Function Name: ExonAssistedSPALN
#
#  About: This function takes a region of sequence (determined through
#         considerations of the gtf index) and searches for a hit in
#         that region using spaln.
#
sub ExonAssistedSPALN
{
    my $ChromosomeRef = shift;
    my $Chromosome    = $$ChromosomeRef;

    # In case we're working with the reverse complement, we want to know
    my $revcomp = 0;
    my $chrname = $Chromosome->{ChrName};
    if ($chrname =~ /\[revcomp\]/) {
	$chrname =~ s/\[revcomp\]//;
	$revcomp = 1;
    }

    my $ChrLengthsRef = shift;
    my %ChrLengths    = %{$ChrLengthsRef};

    my $seqname = shift;

    my $genomefilename = shift;

    my $protfilename = shift;
    my $nuclfilename = shift;
    
    my $HitFile = shift;

    my $ProtLength = shift;
    
    my $SpalnMisses = shift;

    # Are we specifically looking for SPALN output?
    my $spalner = shift;

    # Are we timing?
    my $timing     = shift;
    my $timingdata = shift;
    my $timeA;
    $timing = 5 if ($timing);

    # Find min. and max. nucleotide positions based on exons
    my $minNucl = MIN(${$Chromosome->{NuclStarts}}[0],${$Chromosome->{NuclEnds}}[0]);
    my $maxNucl = MAX(${$Chromosome->{NuclStarts}}[0],${$Chromosome->{NuclEnds}}[0]);
    foreach my $exon (1..$Chromosome->{NumHits}-1) {
	$minNucl = MIN($minNucl,MIN(${$Chromosome->{NuclStarts}}[$exon],${$Chromosome->{NuclEnds}}[$exon]));
	$maxNucl = MAX($maxNucl,MAX(${$Chromosome->{NuclStarts}}[$exon],${$Chromosome->{NuclEnds}}[$exon]));
    }

    # Pull in a little extra sequence, just to be safe
    $minNucl = MAX($minNucl-100000,1);
    $maxNucl = MIN($maxNucl+100000,$ChrLengths{$chrname});

    # Toss the selected sequence into our file
    my $eslsfetchCmd;
    $eslsfetchCmd = 'esl-sfetch -c '.$minNucl.'..'.$maxNucl;
    $eslsfetchCmd = $eslsfetchCmd.' -o '.$nuclfilename;
    $eslsfetchCmd = $eslsfetchCmd.' '.$genomefilename.' '.$chrname;
    $eslsfetchCmd = $eslsfetchCmd." > /dev/null 2>&1";
    if (system($eslsfetchCmd)) { die "\n\tERROR: Command '$eslsfetchCmd' failed (EAS)\n\n"; }    
    
    # Assemble the sytem call to run SPALN
    my $spalnCmd = 'spaln -Q3 -O1 -S3 -ya3 '.$nuclfilename.' '.$protfilename;
    $spalnCmd = $spalnCmd.' 2>/dev/null';# if (!$verbose);
    $spalnCmd = $spalnCmd.' |';

    # Run and parse SPALN's output
    my ($hitscore,$hitline) = ParseSPALNOutput($spalnCmd,$revcomp,$minNucl-1,0,$Chromosome->{ChrName},
					       $seqname,$ProtLength,$protfilename,$SpalnMisses,
					       $spalner,$timing,$timingdata);
    @{$timingdata}[6]++ if ($timing);


    # Return your bliss if we got a clean hit, bail if things look bad
    if    ($hitscore)               { return ($hitscore,$hitline); }
    elsif (!$hitscore && !$hitline) { return (0,0);                }

    
    #
    # Trying a bigger chunk of sequence
    #

    
    # If ParseSPALNOutput returns (0,1), it means that we hit at >90%, so
    # how about we expand our borders and see if we can't get up to 97%?

    # Pull in a LOT of extra sequence
    $minNucl = MAX($minNucl-1000000,1);
    $maxNucl = MIN($maxNucl+1000000,$ChrLengths{$chrname});

    # Toss the selected sequence into our file
    $eslsfetchCmd = 'esl-sfetch -c '.$minNucl.'..'.$maxNucl;
    $eslsfetchCmd = $eslsfetchCmd.' -o '.$nuclfilename;
    $eslsfetchCmd = $eslsfetchCmd.' '.$genomefilename.' '.$chrname;
    $eslsfetchCmd = $eslsfetchCmd." > /dev/null 2>&1";
    if (system($eslsfetchCmd)) { die "\n\tERROR: Command '$eslsfetchCmd' failed (EAS)\n\n"; }    
    
    # Assemble the sytem call to run SPALN
    $spalnCmd = 'spaln -Q3 -O1 -S3 -ya3 '.$nuclfilename.' '.$protfilename;
    $spalnCmd = $spalnCmd.' 2>/dev/null';# if (!$verbose);
    $spalnCmd = $spalnCmd.' |';

    # Run and parse SPALN's output
    ($hitscore,$hitline) = ParseSPALNOutput($spalnCmd,$revcomp,$minNucl-1,0,$Chromosome->{ChrName},
					    $seqname,$ProtLength,$protfilename,$SpalnMisses,
					    $spalner,$timing,$timingdata);
    @{$timingdata}[6]++ if ($timing);


    # How's that bliss looking now?
    if ($hitscore) { return ($hitscore,$hitline); }
    else           { return (0,0);                }

}







#########################################################################
#
#  Function Name: BLATAssistedSPALN
#
#  About: Instead of using the gtf index to come up with our search
#         region (like we did in ExonAssistedSPALN) here we use blat
#         to find a likely region for our hit before applying SPALN.
#
sub BLATAssistedSPALN
{
    my $ChrLengthsRef = shift;
    my %ChrLengths    = %{$ChrLengthsRef};

    my $SeqLengthsRef = shift;
    my %SeqLengths    = %{$SeqLengthsRef};

    my $NameIndexRef  = shift;
    my %BlatNameIndex = %{$NameIndexRef};

    my $genomefilename  = shift;
    my $proteinfilename = shift;

    my $BlatFile = shift;

    my $HitFile  = shift;
    my $MissFile = shift;

    my $CmdLog = shift;

    my $InputSpecies = shift;

    my $CPUs = shift;

    # Are we just looking at SPALN output?
    my $spalner = shift;

    # Are we timing?
    my $timing     = shift;
    my $timingdata = shift;
    $timing = 9 if ($timing);

    # Record that we haven't confirmed any hits for any seqs
    my %ConfirmedHits;
    foreach my $seqID (keys %BlatNameIndex) {
	$ConfirmedHits{$BlatNameIndex{$seqID}} = 0;
    }

    # Now we can run blast and parse its output
    my %Chromosomes;
    my %NuclStarts;
    my %NuclEnds;
    my %ProtStarts;
    my %ProtEnds;
    my %Revcomp;
    my %Scores;
    my %BestScore;
    my %BestChr;
    my %TargetChr;
    my %TargetRev;

    # Tracks how many hits scored better than 1e-20 (many indicates possible TE)
    open(my $TEfile,'>',$InputSpecies.'.TEs.Quilter.out');
    my %TEtracker;
    
    # Parse the blast output
    open(my $blatout,'<',$BlatFile);
    while ($line = <$blatout>) {

	# Strip out linebreaks and try to catch junk lines
	$line =~ s/\n|\r//g;
	next if (!$line);
	
	# Get the sequence, chromosome, and the start/end positions in prot. and nucl. seqs
	my $chrname;
	if ($line =~ /^(\S+)\s+(\S+)\s+(\d+)\S*\s+\S+\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	    
	    my $seqID      = $BlatNameIndex{$1};
	    die "  Bad BLAT line: $line ($1 -> '$seqID')\n" if (!$seqID);
	    my $chrname    = $2;
	    my $hitpercent = $3;
	    my $protstart  = $4;
	    my $protend    = $5;
	    my $nuclstart  = $6;
	    my $nuclend    = $7;

	    my $rev = 0;
	    $rev = 1 if ($nuclstart > $nuclend);
	    
	    # Nice hit!
	    if ($hitpercent >= 90) {
		
		$ConfirmedHits{$seqID} = 1;

		if ($Chromosomes{$seqID}) { push(@{$Chromosomes{$seqID}},$chrname); }
		else                      { @{$Chromosomes{$seqID}} = ($chrname);   }
	    
		if ($rev == 0) { #################################### Forward Hit
		    if ($NuclStarts{$seqID}) {
			push(@{$Revcomp{$seqID}},0);
			push(@{$NuclStarts{$seqID}},$nuclstart);
			push(@{$NuclEnds{$seqID}},$nuclend);
		    } else {
			@{$Revcomp{$seqID}}    = (0);
			@{$NuclStarts{$seqID}} = ($nuclstart);
			@{$NuclEnds{$seqID}}   = ($nuclend);
		    }
		    
		} else { ############################################ Reverse Hit
		    if ($NuclStarts{$seqID}) { 
			push(@{$Revcomp{$seqID}},1);
			push(@{$NuclStarts{$seqID}},$nuclend);
			push(@{$NuclEnds{$seqID}},$nuclstart);	
		    } else {
			@{$Revcomp{$seqID}}    = (1);
			@{$NuclStarts{$seqID}} = ($nuclend);
			@{$NuclEnds{$seqID}}   = ($nuclstart);
		    }
		    
		}
		
		# Record the start and end positions in the protein
		if ($ProtStarts{$seqID}) {
		    push(@{$ProtStarts{$seqID}},$protstart);
		    push(@{$ProtEnds{$seqID}},$protend);
		} else {
		    @{$ProtStarts{$seqID}} = ($protstart);
		    @{$ProtEnds{$seqID}}   = ($protend);
		}
		

		# Get the score
		if ($line =~ /(\S+)\s+(\S+)\s*$/) {
		    
		    my $e_value = $1;
		    my $score   = $2;

		    if ($Scores{$seqID}) { 
			push(@{$Scores{$seqID}},$score);
		    } else { 
			@{$Scores{$seqID}} = ($score);
			$BestScore{$seqID} = $score;
			$TargetChr{$seqID} = $chrname;
			$TargetRev{$seqID} = $rev;
		    }

		    if ($score > $BestScore{$seqID}) {
			$BestScore{$seqID} = $score;
			$TargetChr{$seqID} = $chrname;
			$TargetRev{$seqID} = $rev;
		    }

		    if ($e_value =~ /e\-(\d+)/) {
			if ($1 > 20) {
			    if ($TEtracker{$seqID}) {
				$TEtracker{$seqID}++;
			    } else {
				$TEtracker{$seqID} = 1;
			    }
			}
		    } elsif ($e_value =~ /^0\.0/) {
			if ($TEtracker{$seqID}) {
			    $TEtracker{$seqID}++;
			} else {
			    $TEtracker{$seqID} = 1;
			}			
		    }
		}
	    }
	}
    }
    close($blatout);

    # We'll sort these just to make sure things stay constant
    # across our threads.
    my @SeqNames;
    foreach my $seq_num (sort keys %BlatNameIndex) { 
	my $seqname = $BlatNameIndex{$seq_num};
	push(@SeqNames,$seqname);
    }
    my $num_seq_ids = scalar(@SeqNames);

    # If we have fewer SeqIDs than CPUs make that not the case
    if ($num_seq_ids < $CPUs) { $CPUs = $num_seq_ids; }

    # Once again, Thread 0 calls on its friends to save the day
    my $processes = 1;
    my $threadID  = 0;
    my $pid;
    while ($processes < $CPUs) {
	if ($pid=fork) {
	    if (not defined $pid) { die "\n\tFork failed\n\n"; }
	    $threadID=0;
	} else {
	    $threadID=$processes;
	    last;
	}
	$processes++;
    }

    $timing = 0 unless ($threadID == 0);

    # Name all the files each thread will have access to
    my $protfilename  = 'Quilter.BAS.'.$threadID.'.prot.fa';
    my $nuclfilename  = 'Quilter.BAS.'.$threadID.'.nucl.fa';
    my $thread_hits   = 'Quilter.BAS.hits.'.$threadID.'.out';
    my $thread_misses = 'Quilter.BAS.misses.'.$threadID.'.out';
    my $thread_tes    = 'Quilter.BAS.TEs.'.$threadID.'.out';
    open(my $ThreadHitFile,'>',$thread_hits)    || die "\n  ERROR:  Failed to open BLAT hit file '$thread_hits'\n\n";
    open(my $ThreadMissFile,'>',$thread_misses) || die "\n  ERROR:  Failed to open BLAT miss file '$thread_misses'\n\n";
    open(my $ThreadTEFile,'>',$thread_tes)      || die "\n  ERROR:  Failed to open BLAT TE file '$thread_misses'\n\n";

    # Where do you want to start?  End?
    my $division_size = int($num_seq_ids / $CPUs);
    my $first_seq_id  = $threadID * $division_size;
    my $last_seq_id   = (($threadID+1) * $division_size) - 1;
    if ($threadID == $CPUs-1) { $last_seq_id = $num_seq_ids - 1; }

    # Iterate over sequences, SPALNing 'em up and down
    #
    # UPDATE -- Increasingly wide ranges as hits fail to connect
    #
    for (my $seqID_num = $first_seq_id; $seqID_num <= $last_seq_id; $seqID_num++) {

	my $seqname = $SeqNames[$seqID_num];

	# Did we find any suitable hits?
	if ($ConfirmedHits{$seqname}) {

	    # We're going to need that sequence
	    my $eslsfetchCmd = 'esl-sfetch  '.$proteinfilename." '".$seqname."' > ".$protfilename;
	    if (system($eslsfetchCmd)) { die "\n  ERROR:  esl-sfetch command '$eslsfetchCmd' failed\n\n"; }
	    
	    $seqname =~ /([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|(\S*\|)?([^\|\s]+)\s*$/;
	    my $orig_gene = $1;
	    my $isoform   = $2;
	    my $species   = $3;
	    my $iso_id    = $4;
	    my $group_id  = $6;

	    my $seqlength = $SeqLengths{$seqname};

	    # Record this sequence (and number of significant e-values) if
	    # it meets our threshold (100)
	    if ($TEtracker{$seqname} && $TEtracker{$seqname} >= 100) {
		print $ThreadTEFile "$seqname\n";
	    }
	    
	    # Radical dude! Now we run through each chromosome and
	    # generate a list of places we might be interested in
	    # looking.
	    my $BestScore = $BestScore{$seqname};
	    my $chrname   = $TargetChr{$seqname};
	    my $revcomp   = $TargetRev{$seqname};	    

	    # Generate an index that sorts the start positions
	    my @NStarts;
	    my @NEnds;
	    my @PStarts;
	    my @PEnds;
	    my @SeqScores;
	    my $numhits = 0;
	    foreach my $blathit (0..@{$Chromosomes{$seqname}}-1) {

		if (${$Chromosomes{$seqname}}[$blathit] eq $chrname 
		    && ${$Revcomp{$seqname}}[$blathit] == $revcomp) {
		    
		    $NStarts[$numhits]   = ${$NuclStarts{$seqname}}[$blathit];
		    $NEnds[$numhits]     = ${$NuclEnds{$seqname}}[$blathit];
		    $PStarts[$numhits]   = ${$ProtStarts{$seqname}}[$blathit];
		    $PEnds[$numhits]     = ${$ProtEnds{$seqname}}[$blathit];
		    $SeqScores[$numhits] = ${$Scores{$seqname}}[$blathit];

		    $numhits++;

		}
		
	    }
	    
	    
	    my $sortIndexRef = GenSortIndex(\@NStarts);
	    my @sortIndex    = @{$sortIndexRef};


	    # Sort those dang arrays.  Note that start is always less than
	    # end, even if we're in revcomp-land.
	    my @NuclStart;
	    my @NuclEnd;
	    my @ProtStart;
	    my @ProtEnd;
	    my @Score;
	    my $i = 0;
	    while ($i < @sortIndex) {
		$NuclStart[$i] = $NStarts[$sortIndex[$i]];
		$NuclEnd[$i]   = $NEnds[$sortIndex[$i]];
		$ProtStart[$i] = $PStarts[$sortIndex[$i]];
		$ProtEnd[$i]   = $PEnds[$sortIndex[$i]];
		$Score[$i]     = $SeqScores[$sortIndex[$i]];
		$i++;
	    }
	    

	    # Now, what we'll do is scan along the hits looking for breaks of
	    # more than 200,000 characters between the end of one and the start
	    # of the next and search each of these clusters using SPALN, looking
	    # for the best hit we can get in this chromosome (which will be 
	    # equal to BestScore)
	    my @SpalnQueue;
	    my $original_chr =  $chrname;
	    $original_chr    =~ s/\[revcomp\]//;
	    $i = 0;
	    while ($i < @sortIndex) {
		
		# Prime our cluster
		my $minNucl = $NuclStart[$i];
		my $maxNucl = $NuclEnd[$i];
		my $minProt = $ProtStart[$i];
		my $maxProt = $ProtEnd[$i];
		my $AddHit  = 0;
		$AddHit = 1 if ($Score[$i] == $BestScore);

		$i++;
		while ($i < @sortIndex && $NuclStart[$i] < $NuclEnd[$i-1]+200000) {
		    $maxNucl = MAX($maxNucl,$NuclEnd[$i]);
		    $minProt = MIN($minProt,$ProtStart[$i]);
		    $maxProt = MAX($maxProt,$ProtEnd[$i]);
		    $AddHit  = 1 if ($Score[$i] == $BestScore);
		    $i++;
		}

		# Only add this cluster for consideration if its max. score is
		# at least 80% of the max hit to the chromosome
		if ($AddHit) {
		    
		    # Figure out how much of our protein is covered
		    my $coverage = $maxProt-$minProt;
		    
		    # Put it all in a string and add the string to our list
		    push(@SpalnQueue,"$minNucl#$maxNucl#$original_chr#$chrname#$coverage");
		    
		}
	    }


	    # How much extra sequence we'll pull in (successively larger and
	    # slower, but if we're getting any sort of hit at all this should
	    # find a full translated alignment).
	    #
	    # 3:  100 KB
	    # 2:  1   MB
	    # 1:  10  MB
	    # 0:  Bail
	    #
	    my $highscore = 0;
	    my $hitline;
	    my $remaining_attempts = 3;
	    while ($remaining_attempts) {

		$remaining_attempts--;
		my $extra_seq = 100000 * (10 ** (2-$remaining_attempts));
	
		# Run through the list of items
		my $numQueries = @SpalnQueue;
		foreach my $query (@SpalnQueue) {
		    
		    my @QueryArray   = split('#',$query);
		    my $minNucl      = $QueryArray[0];
		    my $maxNucl      = $QueryArray[1];
		    my $original_chr = $QueryArray[2];
		    my $chrname      = $QueryArray[3];
		    my $coverage     = $QueryArray[4];
		    
		    # Pull in a little extra sequence, just to be safe
		    $minNucl = MAX($minNucl-$extra_seq,1);
		    $maxNucl = MIN($maxNucl+$extra_seq,$ChrLengths{$original_chr}-1);
		    
		    # Toss the selected sequence into our file (need to be a
		    # little careful about reverse complement)
		    $eslsfetchCmd = 'esl-sfetch -c '.$minNucl.'..'.$maxNucl;
		    $eslsfetchCmd = $eslsfetchCmd.' -o '.$nuclfilename;
		    $eslsfetchCmd = $eslsfetchCmd.' '.$genomefilename.' '.$original_chr;
		    $eslsfetchCmd = $eslsfetchCmd." > /dev/null 2>&1";
		    if (system($eslsfetchCmd)) { die "\n\tERROR: Command '$eslsfetchCmd' failed (BAS)\n\n"; } 

		    # Assemble the sytem call to run SPALN
		    my $spalnCmd = 'spaln -Q3 -O1 -S3 -ya3 '.$nuclfilename.' '.$protfilename;
		    $spalnCmd = $spalnCmd.' 2>/dev/null';# if (!$verbose);
		    $spalnCmd = $spalnCmd.' |';	
		    
		    print $CmdLog "> (".localtime()." $spalnCmd\n\n";

		    if ($revcomp) { $chrname = $chrname.'[revcomp]'; }
		    
		    # Parse SPALN's output
		    my ($nextScore,$nextLine) = ParseSPALNOutput($spalnCmd,$revcomp,$minNucl-1,$highscore,
								 $chrname,$seqname,$seqlength,$protfilename,
								 $CmdLog,$spalner,$timing,$timingdata);
		    @{$timingdata}[10]++ if ($timing);

		    
		    # Is this better than any hit we've seen?
		    if ($nextScore > $highscore) {
			$highscore = $nextScore;
			$hitline = $nextLine;
		    }

		}

		last if ($highscore);
		#die if ($highscore);

	    }
	    
	    
	    # Print out the hit line and head on home
	    if ($highscore) {
		$hitline =~ s/Method     \: SPALN/Method     \: SPALN\+BLAT/;
		print $ThreadHitFile  "$hitline";
	    } else { 
		print $ThreadMissFile "$seqname\n"; 
	    }
	    
	} else {

	    # No hits from blat for this sequence
	    print $ThreadMissFile "$seqname\n";

	}
    }

    # Clean up after yourself, man!
    #
    if (-e $nuclfilename) { system("rm $nuclfilename"); }
    if (-e $protfilename) { system("rm $protfilename"); }

    # Close up the thread-specific files
    #
    close($ThreadHitFile);
    close($ThreadMissFile);
    close($ThreadTEFile);

    # de-fork
    if ($threadID) { exit(0); }
    while (wait() > -1) { }

    # Now it's up to thread 0 to tie up the loose ends (by iterating
    # through the thread-specific hit and miss files and transferring
    # the data over to the big-deal HitFile and MissFile -- don't forget
    # old TEfile, either!)
    #
    for ($threadID=0; $threadID<$CPUs; $threadID++) {

	$thread_hits = 'Quilter.BAS.hits.'.$threadID.'.out';
	if (-e $thread_hits) {
	    open($ThreadHitFile,'<',$thread_hits) || die "\n  ERROR:  Failed to open BLAT hit file '$thread_hits' (input)\n\n";
	    while (my $line = <$ThreadHitFile>) {
		print $HitFile "$line";
	    }
	    close($ThreadHitFile);
	    system("rm \"$thread_hits\"");
	}

	$thread_misses = 'Quilter.BAS.misses.'.$threadID.'.out';
	if (-e $thread_misses) {
	    open($ThreadMissFile,'<',$thread_misses) || die "\n  ERROR:  Failed to open BLAT miss file '$thread_misses' (input)\n\n";
	    while (my $line = <$ThreadMissFile>) {
		print $MissFile "$line";
	    }
	    close($ThreadMissFile);
	    system("rm \"$thread_misses\"");
	}
	
	$thread_tes = 'Quilter.BAS.TEs.'.$threadID.'.out';
	if (-e $thread_tes) {
	    open($ThreadTEFile,'<',$thread_tes) || die "\n  ERROR:  Failed to open BLAT TE file '$thread_misses' (input)\n\n";
	    while (my $line = <$ThreadTEFile>) {
		print $TEfile "$line";
	    }
	    close($ThreadTEFile);
	    system("rm \"$thread_tes\"");
	}

    }

}






#########################################################################
#
#  Function Name: ParseSPALNOutput
#
#  About: The name says it all.  Run a SPALN command and make sense of
#         its wisdom.
#
#    1.  Don't trust percent identities.
#    2.  J = S.
#
sub ParseSPALNOutput
{
    my ($i,$j,$k);

    my $spalnCmd  = shift;
    my $revcomp   = shift;
    my $offset    = shift;
    my $highscore = shift;

    my $ChrName = shift;

    my $seqname = shift;

    my $prot_len     = shift;
    my $protfilename = shift;
    
    my $SpalnMisses = shift;

    # NOTE: When we're just looking for SPALN output we can't just
    #       run the command and jump ship, due to how the various
    #       functions depend on score comparison (well, BLAT).
    #       For this reason, we'll set the 'hitstring' to be the full
    #       output AFTER we get the score.
    #
    my $spalner = shift;

    # Are we timing?
    my $timing     = shift;
    my $timingdata = shift;
    my $timeA;

    my ($line,$hitstring);

    # DEBUGGING
    #my $printspalnCmd = $spalnCmd;
    #$printspalnCmd =~ s/\|$//;
    #system($printspalnCmd);
    # DEBUGGING

    
    $timeA = [Time::HiRes::gettimeofday()] if ($timing);
    open(my $stdout,$spalnCmd) || die "\n  ERROR: Spaln command '$spalnCmd' failed\n\n";
    @{$timingdata}[$timing] += Time::HiRes::tv_interval($timeA) if ($timing);


    # Look for where SPALN has called the start and end of the region.
    my ($range_low,$range_high);
    $line = readline($stdout);
    while (!eof($stdout) && $line !~ /\S+\/(\d+)\-(\d+)/) {
	$line = readline($stdout);
    }

    # If we've hit the end of the file, no point continuing
    if (eof($stdout)) {
	close($stdout);
	return(0,0); 
    }


    # Grab the range of positions on the chromosome that we hit in.
    # Note that we want these to be formatted low=start, high=end
    $line =~ /\S+\/(\d+)\-(\d+)/;
    my $start_pos = $1;
    my $end_pos   = $2;


    # Go to the section that lays out the score and identity
    #
    # NOTE: Even though we generally know what's going on with the
    #       complementarity going into our adventures with SPALN,
    #       it at least one case there's been a flip.
    #
    my $saw_complement = 0;
    my ($hitscore,$percent,$true_num_chars);
    while ($line = <$stdout>) {

	# Complementarity marker!
	if ($line =~ /\;C complement\(/) {
	    $saw_complement = 1;
	}

	if ($line =~ /Score \= (\d+)/) {
	    
	    # SPALN's score
	    #
	    $hitscore = $1;

	    # It turns out we need to manually confirm that every
	    # character is included, because SPALN loves pranking
	    # me around.
	    #
	    if ($line =~ /Score \= \S+ \S+\, (\d+)\.\d+ \S+\, (\d+)\.\d+ \S+\, (\d+)\.\d+ \S+\, (\d+)\.\d+/) {
		$true_num_chars = $1 + $2 + $3 + $4;
	    } else {

		# It seems that in cases where SPALN gives BAD hits
		# it breaks away from the expected format, so we'll
		# just assume that we want to try again (hence the
		# (0,1) return).
		#
		close($stdout);
		return(0,1);
		
	    }

	    last;
	}
    }


    # Check for any unexpected complementarity stuff
    if ($revcomp) {
	if (!$saw_complement) {
	    $revcomp = 0;
	    $ChrName =~ s/\[revcomp\]//;
	}
    } else {
	if ($saw_complement) {
	    $revcomp = 1;
	    $ChrName = $ChrName.'[revcomp]';	    
	}
    }

    
    # If we've hit the end of the file or have worse than 97% identity
    # jump ship.
    if (eof($stdout)) {
	close($stdout); 
	return(0,0); 
    } elsif ($true_num_chars < $prot_len) {
	close($stdout);
	return(0,1); # <- This '1' might give us a second chance (pull in more sequence)
    }
    
    
    # NOW we can go to the real business (the alignment lines)
    $line = readline($stdout); # 'ALIGNMENT'
    while (!eof($stdout) && $line !~ /ALIGNMENT/) {
	$line = readline($stdout);
    }


    # Still not totally stoked on seeing an eof
    if (eof($stdout)) {
	close($stdout);
	return(0,0);
    }


    #
    # Here's where we'll gather the SPALN output when we're looking
    # to see the actual SPALN output, for debugging.
    #
    if ($spalner) {

	close($stdout);

	open($stdout,$spalnCmd) || die "\n  ERROR:  Spaln command '$spalnCmd' failed\n\n";
	
	$hitstring = 'SEQ: '.$seqname."\n";
	$hitstring = $hitstring.'CMD: '.$spalnCmd."\n";
	while ($line = <$stdout>) { $hitstring = $hitstring.$line; }
	$hitstring = $hitstring."\nEND\n\n";

	close($stdout);

	return($hitscore,$hitstring);
	
    }


    #
    # ALRIGHT, GENTS AND LADIES, it looks like we're actually
    # doing this thing.
    #


    # Scan through the SPALN output constructing 4 arrays:
    #
    #     - An array with the DNA sequence characters
    #     - An array with the protein sequence characters
    #     - An array mapping protein sequence indices to
    #       DNA sequence indices
    #     - An array mapping DNA sequence characters to
    #       positions in the genome
    #
    # The combination of these arrays will allow us to do
    # all of the work we need in terms of generating our
    # final mapping of the protein to the genome, correcting
    # for micro-exons.
    #
    my $mismatches = 0;
    my @FullNuclSeq;
    my @FullProtSeq;
    my @AAPositions;
    my @NuclPositions;    
    my $trans_line;
    my $full_length = 0;
    my $num_aas     = 0; 
    my $current_pos = $offset;
    my $first_jump  = 0;
    while (!eof($stdout)) {

	# Grab the next line, clean it up, make sure it's
	# meaningful.
	#
	$line = readline($stdout);
	$line =~ s/\n|\r//g;
	next if (!$line);


	# If we're skipping around in the genome, adjust
	# the position variable accordingly.
	#
	if ($line =~ /skip\s+(\d+)\s+nt/) {
	    my $skiplen = int($1);
	    if ($revcomp) { $current_pos -= $skiplen; }
	    else          { $current_pos += $skiplen; }
	    next;
	}


	# If we've made it this far but don't match this
	# terminal format then we're looking at SPALN's
	# translation of the nucleotide sequence
	#
	if ($line !~ /\| \S+\/\d+\-\d+/) {
	    $trans_line = $line;
	    next;
	}

	my $nn_line = $line;

	# EXCELLENT!  This is the next line of DNA characters
	#
	my @NextNucls = split(//,$line);

	# We also split up what would be the translation, so that
	# we can evaluate percent ID for ourselves.
	#
	my @TransLine = split(//,$trans_line);

	# The final check we need to do is to add/subtract
	# the initial 'jump' into the hit (but only if this
	# is the very beginning)
	#
	if ($first_jump == 0) {

	    if ($line =~ /\s+(\d+)/) {

		$first_jump = $1;

		# Special catch for BLAT parsing
		if ($revcomp==3) { $current_pos  = ($start_pos+1)-$first_jump; }
		else             { $current_pos += $first_jump;                }

		$first_jump = 1;

	    } else {
		return(0,0);
	    }

	}

	
	# ... which means that the following line is the next
	# line of protein characters.
	#
	$line = readline($stdout);
	$line =~ s/\n|\r//g;
	my @NextAAs = split('',$line);


	# Advance to the actual DNA sequence.
	# NOTE:  it's possible for SPALN to have a line of
	#        '-' characters, which we'd want to throw
	#        away.  BUT we still need to count up the
	#        position indices when we encounter these. << Do we?
	#
	$i = 0;
	#while ($i < @NextNucls && $NextNucls[$i] !~ /[a-zA-Z]|\|/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] =~ /\s/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] !~ /\s/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] =~ /\s/) { $i++; }
	


	# Some sort of SPALN weirdness...
	#
	next if ($i == @NextNucls || $NextNucls[$i] eq '|');


	# Run along the two lines storing all of the relevant
	# information that we can get out of them.
	#
	while ($NextNucls[$i] !~ /\|/) {
	    
	    push(@FullNuclSeq,$NextNucls[$i]);
	    push(@NuclPositions,$current_pos);

	    # If we're centered on a codon, take note.
	    #
	    if ($NextAAs[$i] =~ /\w/) {		
		push(@FullProtSeq,$NextAAs[$i]);
		push(@AAPositions,$full_length);
		$num_aas++;
	    }
	    
	    # Also, note a *SUBSTITUTION* mismatch (sorry for yelling)
	    $mismatches++ 
		if ($NextAAs[$i] =~ /\S/ 
		    && $NextAAs[$i] ne $TransLine[$i] 
		    && !(uc($TransLine[$i]) eq 'J' && uc($NextAAs[$i]) eq 'S') # Sometimes SPALN calls an 'S' a 'J'
		);
	    # I'm using that format because we might see more delightful tricks


	    # Increment the position in the genome
	    #
	    if ($NextNucls[$i] =~ /\w/) {
		if ($revcomp) { $current_pos--; }
		else          { $current_pos++; }
	    }

	    $full_length++;
	    $i++;
	    
	}
	
    }
    

    # Done reading input!  Now we just need to verify that
    # any and all micro-exons get cleaned up!
    #
    close($stdout);


    # If we aren't satisfied with the overall percent identity of the
    # hit go back (possibly to pull in more sequence and try again).
    #
    return (0,1) if ($mismatches/$prot_len > 0.05);
    

    # Selenocysteine does some crazy stuff, man.  I translated some once
    # and I'm still coming down.
    #
    if ($num_aas != $prot_len) {
	
	my ($fps_ref,$aap_ref) = SelenocysteineCheck(\@FullNuclSeq,\@FullProtSeq,\@AAPositions,$protfilename);
	@FullProtSeq = @{$fps_ref};
	@AAPositions = @{$aap_ref};
	$num_aas = scalar(@FullProtSeq);

	return (0,0) if ($num_aas != $prot_len);

    }


    # ********************************************************************
    # *** 'IF ZERO' WILL TOGGLE MICRO EXON SEARCHING OFF *****************
    # ********************************************************************
    #
    #
    #if (0) {

    # Run along SPALN's output looking for 'exons' with
    # improbably low amino acid counts.  Currently, we
    # say that an exon should have at least 4 amino acids.
    #
    my @MicroExonRuns;
    my $micro_start = -1;
    my $micro_end   = 0;
    
    $i = 0;
    while ($i < $num_aas) {


	# Figure out the genome position of the current amino
	# acid, so that we can infer exon boundaries (large
	# changes in nucleotide position).
	#
	my $aa_loc   = $NuclPositions[$AAPositions[$i]];
	my $aa_start = $i;
	my $exon_len = 1;

	
	# Run through this exon, counting how many amino acids
	# it has.
	#
	$i++;
	while ($i < $num_aas && abs($NuclPositions[$AAPositions[$i]] - $NuclPositions[$AAPositions[$i-1]]) < 5) {
	    $exon_len++;
	    $i++;
	}


	# If we have too few amino acids, then we're working
	# with a micro exon, baby!  Observe that we chain together
	# runs of contiguous micro-exons.
	#
	if ($exon_len < 4) {
	    
	    if ($micro_start == -1) { $micro_start = $aa_start; }
	    $micro_end = $i;
	    
	} elsif ($micro_start != -1) {
	    
	    push(@MicroExonRuns,$micro_start.'-'.$micro_end);
	    $micro_start = -1;
	    
	}
	
    }
    
    
    # Next we'll run through each of the indicated micro-exon
    # regions ('regions' because they can be chained together)
    # and try to re-align the amino acids to the ends of the
    # 'major' exons that flank them.
    #
    # NOTE: Be careful that we don't step over a breakpoint by
    #       checking whether NuclPositions[i] and NP[i+1] are
    #       1 off from one another.
    #
    foreach my $micro_run_str (@MicroExonRuns) {

	
	# Split the indices out of the string
	#
	my @MicroRun = split(/\-/,$micro_run_str);
	$micro_start = $MicroRun[0];
	$micro_end   = $MicroRun[1]; # One past actual end, so we can do strict '<'


	# Identify the amino acid sequence that constitutes
	# this micro exon (or chain of micro exons)
	#
	my @MicroAAs;
	for ($j = $micro_start; $j < $micro_end; $j++) {
	    push(@MicroAAs,$FullProtSeq[$j]);
	}

	# Figure out how long the micro exon is (keeping in mind that '_end' is
	# one position past the actual end index)
	#
	my $micro_len = $micro_end-$micro_start;

	
	# Look at the preceding exon
	#
	my @RearExt;
	my $rear_ext_len = 0;
	if ($micro_start > 0) {

	    
	    # Get situated so that we can start trying to
	    # transfer amino acids from our micro exons to
	    # the back end of the preceding exon.
	    #
	    # NOTE (for making sense of this):
	    # --------------------------------
	    # NuclPositions are the chromosome indices, AAPositions are the
	    # array indices that SPALN has identified as codon centers -- thus
	    # looking up the NuclPosition of an AAPosition gets you the genome
	    # index for the AAPosition's codon center.
	    #
	    my $rear_codon_center = $NuclPositions[$AAPositions[$micro_start-1]];
	    my $nuclseq_pos = $AAPositions[$micro_start-1]+3; # Which array index are we centered around?
	    while ($rear_ext_len < $micro_len) {
		
		# Make sure we don't overstep (by checking if we've passed over
		# an stretch of genome indices indicating an intron)
		#
		my $step_check = $NuclPositions[$nuclseq_pos];
		if ($revcomp) { $step_check += 3 * ($rear_ext_len+1); }
		else          { $step_check -= 3 * ($rear_ext_len+1); }
		
		last if ($step_check != $rear_codon_center);


		# Figure out what the next codon is
		#
		my @NextCodon = ();
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos-1]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos+1]);

		
		# If we have a match we are very happy.  If
		# we do not have a match we are very sad.
		#
		# Why would we do 'micro_start+(rear_ext_len+1)'?
		#
		if (uc($FullProtSeq[$micro_start+$rear_ext_len]) eq uc(TranslateCodon(\@NextCodon))) {
		    push(@RearExt,$nuclseq_pos);
		    $nuclseq_pos += 3;
		    $rear_ext_len++;
		} else {
		    last;
		}
		
	    }
	    
	}
	
	# If we found a complete run, we're sooooo stoked!
	#
	if ($rear_ext_len == $micro_len) {

	    # Record our triumph
	    #
	    for ($i=0; $i<$micro_len; $i++) { $AAPositions[$micro_start+$i] = $RearExt[$i]; }
	    next;
	    
	}

	
	# Look at the exon ahead
	#
	my @FwdExt;
	my $fwd_ext_len = 0;
	if ($micro_end < $num_aas) {
	    
	    # Try to find a good extension to the front of
	    # the upcoming exon.
	    #
	    my $fwd_codon_center = $NuclPositions[$AAPositions[$micro_end]];
	    my $nuclseq_pos = $AAPositions[$micro_end]-3;
	    while ($fwd_ext_len < $micro_len) {

		
		# Make sure we don't overstep
		#
		my $step_check = $NuclPositions[$nuclseq_pos];
		if ($revcomp) { $step_check -= 3 * ($fwd_ext_len+1); }
		else          { $step_check += 3 * ($fwd_ext_len+1); }
		
		last if ($step_check != $fwd_codon_center);


		# What does the next codon encode?
		#
		my @NextCodon = ();
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos-1]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos+1]);
		

		# If we have a match, in terms of what the codon
		# encodes, then we can continue our considerations.
		# Otherwise we cannot.
		#
		# Note that the first nucleotide position in the
		# 'FwdExt' array corresponds to the last amino acid
		# in the micro-exon region.
		#
		if (uc($FullProtSeq[$micro_end-($fwd_ext_len+1)]) eq uc(TranslateCodon(\@NextCodon))) {
		    push(@FwdExt,$nuclseq_pos);
		    $nuclseq_pos -= 3;
		    $fwd_ext_len++;
		} else {
		    last;
		}
		
	    }
	    
	}

	
	# Did we find a complete run?
	#
	if ($fwd_ext_len == $micro_len) {
	    
	    for ($i=0; $i<$micro_len; $i++) { $AAPositions[$micro_end-($i+1)] = $FwdExt[$i]; }
	    next;
	    
	}
	

	# Append as much as you can to the preceding
	# major exon first.  We don't put any work into
	# evaluating whether amino acids that could go
	# either way should be fitted to one or the other
	# side.
	#
	$i = 0;
	while ($i < $rear_ext_len) {
	    $AAPositions[$micro_start+$i] = $RearExt[$i];
	    $i++;
	}

	
	# Skip ahead to the stuff we're appending to the
	# upcoming major exon, if you need to.
	#
	if ($i < $micro_len - $fwd_ext_len) { $i = $micro_len - $fwd_ext_len; }


	# Append this stuff onto the forward friend.
	#
	$j = $fwd_ext_len - 1;
	while ($i < $micro_len) {
	    $AAPositions[$micro_end-($micro_len-$i)] = $FwdExt[$j];
	    $i++;
	    $j--;
	}
	
    }
    

    #}
    #
    #
    # ********************************************************************
    # *** 'IF ZERO' TOGGLE END *******************************************
    # ********************************************************************

        
    # Radical, dude!  Now all we need to do is just convert
    # our list of amino acid positions into a list of
    # their corresponding middle amino acids (with '*'s to
    # denote splice site boundaries).
    #
    # We'll try to be careful about how we place codons that are
    # called as insertions into the genome.
    #
    my $prev_codon_c = $NuclPositions[$AAPositions[0]];
    my $CodonCenters = $prev_codon_c;
    for ($i = 1; $i < $num_aas; $i++) {

	# A gap of more than 4 is taken to deserve a splice site marker
	#
	if (abs($NuclPositions[$AAPositions[$i]] - $NuclPositions[$AAPositions[$i-1]]) > 4) {
	    $CodonCenters = $CodonCenters.',*';
	}

	# Prepping for a fresh character!
	#
	$CodonCenters = $CodonCenters.',';
	
	# Are we in an insertion (wrt genome)?
	#
	if ($FullNuclSeq[$AAPositions[$i]] !~ /\w/) {
	    $CodonCenters = $CodonCenters.$prev_codon_c;
	} else {
	    $prev_codon_c = $NuclPositions[$AAPositions[$i]];
	    $CodonCenters = $CodonCenters.$prev_codon_c;
	}
	
    }
    

    # Now that we've got our purty little hit, we write it on out and tell
    # 'em how happy that made us feel.
    my $hitline = "Isoform ID : $seqname\n";
    $hitline    = $hitline."Method     : SPALN\n";
    $hitline    = $hitline."Chromosome : $ChrName\n";
    $hitline    = $hitline."Match Pos.s: $CodonCenters\n\n";


    return ($hitscore,$hitline);

}





#########################################################################
#
#  Function Name: SelenocysteineCheck
#
#  About: This function checks if a selonocyteine (U, encoded by TGA)
#         explains an observed absence of aminos in SPALN output.
#
sub SelenocysteineCheck
{
    my $nucl_seq_ref = shift;
    my $prot_seq_ref = shift;
    my $aa_pos_ref   = shift;
    my $protfilename = shift;

    my @NuclSeq = @{$nucl_seq_ref};
    my @ProtSeq = @{$prot_seq_ref};
    my @AAPos   = @{$aa_pos_ref};

    my @TrueProtSeq;
    open(my $protfile,'<',$protfilename) || die "\n  ERROR:  Failed to open file '$protfilename' during selenocysteine check\n\n";
    while (my $line = <$protfile>) {
	$line =~ s/\n|\r//g;
	next if ($line =~ /\>/);
	foreach my $amino (split(//,$line)) { push(@TrueProtSeq,uc($amino)); }
    }
    close($protfile);


    my $aa_pos = 0;
    while ($aa_pos < scalar(@TrueProtSeq)) {

	# Looks like we've discovered a selenocysteine in our midst!
	if ($TrueProtSeq[$aa_pos] eq 'U' && $aa_pos < scalar(@ProtSeq) && $ProtSeq[$aa_pos] ne 'U') {

	    # Is this consistent with stepping forward 3 from the previous
	    # amino's position?  3 back from forthcoming amino?
	    my $nailed_it = 0;
	    if ($aa_pos) {
		
		# Should be to the right of the last position...
		my $nucl_pos = $AAPos[$aa_pos-1] + 2;
		my $codon = uc($NuclSeq[$nucl_pos].$NuclSeq[$nucl_pos+1].$NuclSeq[$nucl_pos+2]);
		if ($codon eq 'TGA') {

		    # Note:  Splice moves everything at the indicated position to
		    #        the right and sticks in what you want there.
		    splice(@ProtSeq,$aa_pos,0,'U');
		    splice(@AAPos,$aa_pos,0,$nucl_pos+1);
		    $nailed_it = 1;
		} 

	    } 

	    if (!$nailed_it && $aa_pos < scalar(@ProtSeq)-1) {
		
		# ... Or left of the next position
		my $nucl_pos = $AAPos[$aa_pos] - 4;
		my $codon = uc($NuclSeq[$nucl_pos].$NuclSeq[$nucl_pos+1].$NuclSeq[$nucl_pos+2]);
		if ($codon eq 'TGA') {
		    splice(@ProtSeq,$aa_pos-1,0,'U');
		    splice(@AAPos,$aa_pos-1,0,$nucl_pos+1);
		} 

	    }

	}

	$aa_pos++;

    }

    return (\@ProtSeq,\@AAPos);

}





#########################################################################
#
#  Function Name: LooksRepetitive
#
#  About: This function is used to check proteins that failed to align
#         well to the genome for signs of high repetition.  This might
#         need fine-tuning in the future (-r value / threshold).
#
sub LooksRepetitive
{
    my $threshold = 30;
    my $protfilename = shift;
    
    open(my $stdoutput,"./TransSW $protfilename \- \-self \|");
    
    my $longest_run = <$stdoutput>;
    $longest_run =~ s/\r|\n//g;
    if ($longest_run >= $threshold) { return 1; }
    else                            { return 0; }

}






#########################################################################
#
#  Function Name:  GenSortIndex
#
#  About:  Figure out how an array should be sorted (note that this does
#          not sort the array, but instead provides the order to access
#          elements in that array).
#
sub GenSortIndex
{

    my $targetRef = shift;
    my @Target    = @{$targetRef};

    my @Index;
    my $length = scalar(@Target);
    for (my $i=0; $i<$length; $i++) { $Index[$i] = $i; }

    my $group_size = 1;
    while ($group_size < $length) {

	my $num_groups = POSIX::ceil($length / $group_size);
	for (my $group_id = 0; $group_id < $num_groups; $group_id += 2) {

	    my $start1 = $group_size * $group_id;
	    my $start2 = $group_size * ($group_id+1);

	    next if ($start2 > $length);

	    my $end1 = $start2;
	    my $end2 = ($group_size * ($group_id+2));

	    if ($end2 > $length) { $end2 = $length; }

	    my $run1 = $start1;
	    my $run2 = $start2;

	    my @Temp;
	    while ($run1 < $end1 && $run2 < $end2) {
		if ($Target[$Index[$run1]] < $Target[$Index[$run2]]) {
		    push(@Temp,$Index[$run1]);
		    $run1++;
		} else {
		    push(@Temp,$Index[$run2]);
		    $run2++;
		}
	    }

	    while ($run1<$end1) {
		push(@Temp,$Index[$run1]);
		$run1++;
	    }

	    while ($run2<$end2) {
		push(@Temp,$Index[$run2]);
		$run2++;
	    }

	    for (my $i=0; $i<scalar(@Temp); $i++) {
		$Index[$start1+$i] = $Temp[$i];
	    }

	}

	# Increase the group size and repeat
	$group_size *= 2;

    }
    
    # Great! Now pass that bad-boy back!
    return \@Index;

}




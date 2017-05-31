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
use File::Basename;
use lib dirname (__FILE__);
use Cwd;
use DiagonalSets;



sub MAX;
sub MIN;
sub TranslateCodon;
sub BuildGeneIndex;
sub AddExon;
sub GetChromosomeLengths;
sub RunFindDiagonals;
sub FindFullAlignments;
sub ListCodonCenters;
sub GetBestExons;
sub AddToHitString;
sub CopyHit;
sub ExonAssisteSPALN;
sub BLATAssistedSPALN;
sub ParseSPALNOutput;
sub AdjustSPALNOutput;
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
    print "  USAGE:  \$ perl  Quilter.pl  <isoforms>  <genome>  <gtf index>  <species>  [options]\n\n"; 
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
my $timer   = 0; # Set a timer.
my $overw   = 0; # Don't overwrite (default)
my $resdir  = 0; # Use a user provided results directory
my $CPUs    = 2; # Default: use 2 cores
my $verbose = 0; # Use verbose printing
my $spalner = 0; # Are we just running this as a means for getting SPALN output?


# Check for additional options
my $i = $expectedArgs;
while ($i < @ARGV) {
    if ($ARGV[$i] eq '-debug') {
	$debug = 1;
    } elsif ($ARGV[$i] eq '-overwrite') {
	$overw = 1;
    } elsif ($ARGV[$i] eq '-o') {
	$i++;
	$resdir = $i;
	$overw  = 1;
    } elsif ($ARGV[$i] eq '-setcwd') {
	$i++;
	chdir($ARGV[$i]);
    } elsif ($ARGV[$i] eq '-v') {
	$verbose = 1;
    } elsif ($ARGV[$i] eq '-n') {
	$i++;
	$CPUs = int($ARGV[$i]);
	if ($CPUs <= 0) {
	    print "\n\tUnsupported number of CPUs requested ($CPUs)\n";
	    print "\tReverting to default (2)\n\n";
	    $CPUs = 2;
	}
    } elsif ($ARGV[$i] eq '-SPALN') {
	$spalner = 1;
    } else {
	print "\n\tUnrecognized option $ARGV[$i] ignored.\n\n";
    }
    $i++;
}


# Check if 'SPALN' is installed (and, if so, check if BLAST is happy)
open(my $spalncheck,'which spaln |');
my $spaln = <$spalncheck>;
close($spalncheck);


# Well... is it? (Also, make a file to catch SPALN deaths)
my %ChrLengths;
if (!$spaln) {
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


# If we're running a timer, start it now.  If we aren't, start it anyways.
my $startTime = time;
my $lastTime;
my $timeDiffMins;
my $timeDiffSecs;


# Report that we're getting started
my $progmsg = '  Quilter.pl:  Preparing index for '.lc($ARGV[3]);
while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
print "$progmsg\r" if (!$verbose);


# Construct an index on the genes in the isoform file, extracted from the gtf file.
# Additionally, get a list of all genes that did not have any matches in the gtf file.
my $reqSpecies = uc($ARGV[3]);
my %GeneIndex;
my %Blacklist;
my $AutoBLAT = 0;
if ($ARGV[2] ne '-') {
    my ($indexRef,$blacklistRef) = BuildGeneIndex($ARGV[0],$ARGV[2],$reqSpecies,$debug);
    %GeneIndex = %{$indexRef};
    %Blacklist = %{$blacklistRef};
} elsif ($blat) {
    $AutoBLAT = 1;
} else {
    die "\n  ERROR: Without blat species '$reqSpecies' requires valid GTF index.\n\n";
}

# TIMER SET TO 'ON'
# It might be worth knowing how long it took to get our index squared away.
$lastTime = time;
$timeDiffMins = int(($lastTime-$startTime) / 60);
$timeDiffSecs = ($lastTime-$startTime) % 60;
print "\n  Time to build gene index: $timeDiffMins minutes $timeDiffSecs seconds\n\n" if ($verbose);


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
    if (!(-e $foldername) && system("mkdir $foldername")) { die "\n\tFailed to generate temporary folder '$foldername'\n\n"; }
}


# To guarantee that our threads don't overlap, we try to find a
# good mid-point in the isoform file to split around.
open(my $lineCounter,"wc -l $ARGV[0] \|") || die "\n\tLine counting failed.  Is '$ARGV[0]' the correct isoform file location?\n\n";
my $lineCountLine = <$lineCounter>;
if (! ($lineCountLine =~ /^\s*(\d+).*/)) {
    die "\n\tUnexpected output format from line counter: $lineCountLine\n\n";
}
my $LineCount = $1;
if ($LineCount < $CPUs) { $CPUs = $LineCount; }
my $portion   = $LineCount/$CPUs;


# Non-0 threads report results to a file for thread 0
my $resultsfilehead = $foldername.'temp_results.'; 
my $resultsfile;


# Filename used to track progress of threads
my $progressbase = 'quilter.thread_progress.';


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
my $protfilename    = $foldername.'prot_tempfile_'.$threadID.'.Quilter.fa';
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
system("rm $protfilename") if (-e $protfilename);
system("rm $hitfilename")  if (-e $hitfilename);
system("rm $missfilename") if (-e $missfilename);
system("rm $blatfilename") if (-e $blatfilename);


# A log file used to track how well spaln performs
open(my $SpalnLog,'>',$foldername.'spaln_logfile_'.$threadID.'.Quilter.out');


# A file to hold all our hits and a file to hold our unhit sequence info.
open(my $HitFile,'>',$hitfilename);
open(my $MissFile,'>',$missfilename);


#my $genehits   = 'GENE_HITS_'.$threadID;
#my $genemisses = 'GENE_MISSES_'.$threadID;
#open(my $GeneHitFile,'>',$genehits);
#open(my $GeneMissFile,'>',$genemisses);


# Open the file containing all of the isoforms
open (my $isoformfile,'<',$ARGV[0]) || die "\n\tCould not open isoform(s) file '$ARGV[0]'\n\n";
my $line = <$isoformfile>;
my $lineNum = 1;


# If we aren't thread 0, we get at least 'midpoint' lines into the
# file before we start reading anything.
while ($lineNum < ($threadID*$portion)) {
    $line = <$isoformfile>;
    $lineNum++;
}


# Used to track how things went down
my $CurrentStat;
my @HitStats;
$HitStats[0] = 0; # Not indexed
$HitStats[1] = 0; # Wrong species
$HitStats[2] = 0; # Indexed, but no hits from FindDiagonals.c
$HitStats[3] = 0; # FindDiagonals.c gave us enough to score
$HitStats[4] = 0; # SPALN saves the day
$HitStats[5] = 0; # Failure -- not enough to score
$HitStats[6] = 0; # How many were just 1 or 2 end AAs off?
$HitStats[7] = 0; # How many needed blast assistance?
my $repetitive;
my @RepConf; # Repetitiveness Confusion Matrix
$RepConf[0][0] = 0;
$RepConf[0][1] = 0;
$RepConf[1][0] = 0;
$RepConf[1][1] = 0;


# RAD! Time to find happy alignments!
my $numHits = 0;
my $numIsos = 0;
my $stoppoint = ($threadID+1) * $portion;
if ($threadID == $CPUs-1) { $stoppoint = $LineCount; }
while (!eof($isoformfile) && $lineNum < $stoppoint) {
    
    # Reset stat tracker
    $CurrentStat = 0;

    # Eat any comment lines
    while (!eof($isoformfile) && !($line =~ m/^\>/)) {
	$line = <$isoformfile>;
	$lineNum++;
    }
    last if ($lineNum >= $stoppoint || eof($isoformfile));
    

    # Strip away line breaks
    $line =~ s/\r|\n//g;

    # Make grabbing info. easy
    $line =~ m/^\>GN\:([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|(\S*\|)?([^\|\s]+)\s*$/;
    
    # If we couldn't match one or more fields, let me know.
    if (!$1 || !$2 || !$3 || !$4 || !$6) {
	print "  Strange line: $line\n";
	$line = <$isoformfile>;
	$lineNum++;
	next;
    }

    # Record sequence name information
    my $orig_gene = $1;
    my $gene      = uc($1);
    my $species   = uc($3);
    my $seqline   = $line;
    $seqline      =~ s/\n|\r//g;
    $seqline      =~ s/\>GN\://;
    $numIsos++;
    
    
    # Change any whitespaces to '_'s
    $gene    =~ s/\s/\_/g;
    $species =~ s/\s/\_/g;
    $seqline =~ s/\s/\_/g;
    
    
    # If this is the wrong species, skip it
    if (uc($species) ne $reqSpecies) {
	
	# Record that we had a wrong species
	$CurrentStat = 1;
	$HitStats[1]++;
	
	# Next up!
	$line = <$isoformfile>;
	$lineNum++;
	next;
	
    }
    
    
    # Individual timer
    my $isotime = time;
    
    
    # We write each isoform to a file, which is then passed to FindDiagonals.
    # We'll also hold onto a copy of the protein sequence for later reference.
    open(my $protfile,'>',$protfilename);
    print $protfile "$line\n";
    my $Protein = '';
    $line =  <$isoformfile>;
    $lineNum++;
    while ($line && $line !~ /^\>/) {
	$line =~ s/\n|\r//g;
	print $protfile "$line\n";
	$Protein = $Protein.$line;
	$line =  <$isoformfile>;
	$lineNum++;
    }
    my $ProtLength = length($Protein);
    close($protfile);    
    
    
    # If this is a blacklisted gene or we're "AutoBLAT"-ing, cut to the
    # chase with using BLAST
    if (($AutoBLAT || $Blacklist{$gene}) && $spaln && $blat) {
	if (system("cat $protfilename >> $blatfilename")) { die "\n  Concatenation of seq. failed\n\n"; }
	
    } else {
	
	# Does this seq. look repetitive?
	$repetitive = LooksRepetitive($protfilename) unless ($spalner);

	# In the event of multiple chromosomes having strong hits,
	# we'll grab the best `fast diagonals' hit (trusting manual
	# annotation over SPALN) or best `SPALN' hit if we can't
	# get a good fast diags hit.
	my $best_fast_score   = 0;
	my $best_fast_string  = '';
	my $best_SPALN_score  = 0;
	my $best_SPALN_string = '';
	my $next_score;
	my $next_string;
	
	# Let's do some actual work!  Go through the gene index, considering each
	# chromosome that the gene has been indexed to, creating 'exon' objects for
	# each pair of start and end points, and then trying to find a full match to
	# the protein.
	foreach my $chromosomeName (sort keys %{$GeneIndex{$gene}}) {

	    # We don't want to bother trying to align to a chromosome that
	    # doesn't correspond to any of the sequences in the genome file
	    # that we've been given, but we're going to miss quite a few if
	    # we don't correct for 'revcomp' indicator
	    my $raw_chr_name = $chromosomeName;
	    $raw_chr_name =~ s/\[revcomp\]//g;
	    next if (!$ChrLengths{$raw_chr_name});

	    #print "CHROMOSOME: $chromosomeName ($gene)\n";
	    
	    # Add this chromosome to the chromosome set (hash).
	    my $Chromosome = New DiagonalSet($chromosomeName);
	    
	    
	    # Iterate over all entries for this gene on this chromosome.
	    my $i = 0;
	    while ($i < @{$GeneIndex{$gene}{$chromosomeName}}) {
		
		# Grab the next start and end pair
		my $start = $GeneIndex{$gene}{$chromosomeName}[$i];
		my $end   = $GeneIndex{$gene}{$chromosomeName}[$i+1];		
		
		# Generate a new exon object
		AddExon(\$Chromosome,$protfilename,$ARGV[1],$nuclfilename,
			$gene,$chromosomeName,$ProtLength,$start,$end,$debug);
		
		# Next pair!
		$i += 2;
		
	    }
	    
	    
	    # If none of the indices yielded hits, we do something EXTRA special (blat)
	    if ($Chromosome->{NumHits}) {
		
		# Sort the chromosome's exons by their start positions
		$Chromosome->SortHitsByField('ProtStarts',1); 
		
		# Search for any full diagonals (that cover entire protein), if it makes sense.
		if (@{$Chromosome->{ProtStarts}} && ${$Chromosome->{ProtStarts}}[0] == 0 && !$spalner) {
		    ($next_score,$next_string) = FindFullAlignments(\$Chromosome,$Protein,$seqline,$HitFile,$debug);
		    
		    #print $GeneHitFile "$gene\n";
		    
		    # GOT FISH!
		    if ($next_score > $best_fast_score) {
			$best_fast_score  = $next_score;
			$best_fast_string = $next_string;
			$CurrentStat = 3;
		    }
		    
		}
		
		# If we don't have a better hit on this chromosome punt over to SPALN
		if ($Chromosome->{NumHits} && $CurrentStat != 3 && $spaln) {
		    
		    ($next_score,$next_string) = ExonAssistedSPALN(\$Chromosome,\%ChrLengths,$seqline,
								   $ARGV[1],$protfilename,$nuclfilename,
								   $HitFile,$ProtLength,$SpalnLog,$spalner);
			
		    # SPALN saves the day!!!
		    if ($next_score > $best_SPALN_score) {
			$best_SPALN_score  = $next_score;
			$best_SPALN_string = $next_string;
			$CurrentStat = 4 if ($CurrentStat == 0);
		    }
		    
		}
		
	    } 
	    
	    #### Time for another chromosome! ####
	    
	}


	# Try setting aside anything that still hasn't given us a hit for BLAT
	if ($CurrentStat == 0) {
	    
	    # Do we have the technology?
	    if ($spaln && $blat) {		
		if (system("cat $protfilename >> $blatfilename")) { die "\n  Concatenation of seq. failed\n\n"; }
		
	    } else {
		
		# What a terrible loss :'(
		print $MissFile ">GN:$seqline\n";
		#print $GeneMissFile "$gene\n";		    
		$HitStats[5]++;
		$CurrentStat = 5;
		
	    }
	    
	} else {

	    $numHits++;
	    $HitStats[$CurrentStat]++;
	    if ($CurrentStat == 3) { print $HitFile "$best_fast_string";  }
	    else                   { print $HitFile "$best_SPALN_string"; }
	    
	}
	    
    }

    # TIME!
    $isotime = time - $isotime;

    # (FOR CHECKING)
    # We report how this isoform went
    my ($repX,$repY);
    if ($CurrentStat) { # Because we might be blat-ing some, 0 is 'skip'
	if ($CurrentStat == 3) {
	    $CurrentStat = "diagonals  " if ($verbose);
	    $repX = 1;
	} elsif ($CurrentStat == 6) {
	    $CurrentStat = "close diag " if ($verbose);
	    $repX = 1;
	} elsif ($CurrentStat == 4) {
	    $CurrentStat = "spaln      " if ($verbose);
	    $repX = 1;
	} elsif ($CurrentStat == 5) {
	    $CurrentStat = "bummerzone " if ($verbose);
	    $repX = 0;
	}
	if ($repetitive) {
	    $CurrentStat = $CurrentStat.' [    repetitive]' if ($verbose);
	    $repY = 0;
	} else {
	    $CurrentStat = $CurrentStat.' [NOT repetitive]' if ($verbose);
	    $repY = 1;
	}
	my $nameline = "  >GN:$seqline";
	while (length($nameline) < 50) {
	    $nameline = $nameline.' ';
	}
	print "$nameline $CurrentStat -- $isotime seconds\n" if ($verbose);
	$RepConf[$repX][$repY]++;
    }

    
    # Record progress, check if it's time to post a status update
    if (!$verbose) {
	$num_complete++;
	if (time()-$progtime > 2) {
	    
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

#close($GeneHitFile);
#close($GeneMissFile);

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
    close($resultsfile);
}


# No need for these files anymore.
system("rm $nuclfilename") if (-e $nuclfilename);
system("rm $protfilename") if (-e $protfilename);

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


# How long did this take?
$lastTime = time;
$timeDiffMins = int(($lastTime-$startTime) / 60);
$timeDiffSecs = ($lastTime-$startTime) % 60;
print "\n  1st Phase Quilter runtime: $timeDiffMins minutes $timeDiffSecs seconds\n\n" if ($verbose);


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
	$line = int($line);
	$HitStats[$i] += $line;
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
	if (system($blatCmd)) { die "\n  ERROR: BLAT command '$blatCmd' failed\n\n"; }

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
	BLATAssistedSPALN(\%ChrLengths,\%SeqLengths,\%BlatNameIndex,$ARGV[1],$ARGV[0],$BlatResults,$Results,$Misses,$SpalnLog,$ARGV[3]);
	close($Results);
	close($Misses);

	# Clean up
	system("rm $BlatResults") if (-e $BlatResults);
	system("rm $bigBlat");

    }
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


# TIMER SET TO 'ON'
# Report total runtime
$lastTime = time;
$timeDiffMins = int(($lastTime-$startTime) / 60);
$timeDiffSecs = ($lastTime-$startTime) % 60;
print "\n  Overall Quilter runtime: $timeDiffMins minutes $timeDiffSecs seconds\n\n" if ($verbose);


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
#
sub BuildGeneIndex
{
    # Get the names of the files we'll be working with
    my $isoformfilename = shift;
    my $indexfilename   = shift;

    # What species do we require?
    my $species = shift;

    # Are we debugging?
    my $debug = shift;

    # Are we using verbose output?
    my $verbose = shift;
    
    # Open the file containing the isoforms and note all gene names in
    # a hash.
    print "\n" if ($debug);
    my %ProtNames;
    my %GeneNames;
    #my %ChromosomeNames; # CHROMOSOME UNIQUENESS CHECK

    open(my $isoformfile,'<',$isoformfilename) || die "\n\tCouldn't open '$isoformfilename'\n\n";
    while (my $line = <$isoformfile>) {

	next unless ($line =~ /^\>GN:([^\|]+)\|([^\|]+)\|([^\|]+)\|/);

	# Grab the pertinent details
	my $seqfam   = $1;
	my $seqiso   = $2;
	my $seqspec  = $3;
	$line =~ s/\n|\r//g;
       
	# Does this entry match the species that we're looking for?
	if (uc($seqspec) eq $species) {

	    # NOTE: switching to all upper-case entries and changing whitespaces to '_'
	    my $pName =  uc($seqfam);
	    my $gName =  uc($seqiso);
	    $pName    =~ s/\s/\_/g;
	    $gName    =~ s/\s/\_/g;

	    # Add to the hash of protein names
	    if (!$ProtNames{$pName}) { $ProtNames{$pName} = [1,$gName];    }
	    else                     { push(@{$ProtNames{$pName}},$gName); }
	    
	    # Add to the hash of gene names
	    if (!$GeneNames{$gName}) { $GeneNames{$gName} = [$pName];      }
	    else                     { push(@{$GeneNames{$gName}},$pName); }
	    
	}
	
    }
    close($isoformfile);
    print "\n" if ($debug);

    
    # Open the gtf index and add all entries that match gene names from
    # our isoform file.  If we're debugging, record all isoforms that we've
    # found index entries for.
    my %GeneIndex;
    my %RepeatTracker;
    open(my $indexfile,'<',$indexfilename);
    while (my $line = <$indexfile>) {
	
	# We don't care about non-exon indices.
	next if $line =~ /^#/;	
	next unless $line =~ /\S+\s+\S+\s+exon/;
	$line =~ /(\S+).+?(\d+)\s+(\d+).+?([\+\-]).+?gene_name "(\S+)"/;

	# I think all of the special-character cases have been covered, but it's
	# good to be sure.
	if (!$5 || !$4 || !$3 || !$2 || !$1) {
	    #print "\n\tWARNING: Line failed to match expected pattern ";
	    #print "(usually caused by special characters in 'gene_name'):\n$line\n\n";
	    next;
	}
	
	my $chrName   = $1;
	my $endPos    = $3;
	my $startPos  = $2;
	my $revcomp   = 0;
	if ($4 eq '-') {
	    $revcomp = 1;
	    $chrName = $chrName.'[revcomp]';
	}
	my $indexName = uc($5);


	# ACTUALLY BUILDING THE HASH:
	#
	# As you can see, the GeneIndex is a hash of hashes, with the purported protein
	# name (Phosphosite) as a primary key and chromosome name as a secondary key.
	#
	# Entries should be read in pairs.	
	#
	#
	if (!$RepeatTracker{$chrName}{$startPos}{$endPos}{$indexName} && $endPos-$startPos > 3) {
	    if ($ProtNames{$indexName}) {
		
		push(@{$GeneIndex{$indexName}{$chrName}},($startPos,$endPos));
		$RepeatTracker{$chrName}{$startPos}{$endPos}{$indexName} = 1;
		$ProtNames{$indexName}[0] = 2;
		
	    } elsif ($GeneNames{$indexName}) { 

		# NOTE: Because we end up searching on what Phosphosite calls the
		#       protein name (even though these tend to align with the gene_name
		#       entries in .gtf) we don't worry about indexing on what
		#       Phosphosite calls the gene name. (Hence the following commenting-out.)
		##push(@{$GeneIndex{$indexName}{$chrName}},($startPos,$endPos));
		$RepeatTracker{$chrName}{$startPos}{$endPos}{$indexName} = 1;
		
		foreach my $protname (@{$GeneNames{$indexName}}) {
		    $ProtNames{$protname}[0] = 2;
		    push(@{$GeneIndex{$protname}{$chrName}},($startPos,$endPos)) if ($protname ne $indexName);
		}
		
		$ProtNames{$indexName}[0] = 2 if ($ProtNames{$indexName});
		
	    }
	}
    }
    close($indexfile);


    # As a sanity check, look up all entries in the GeneNames hash and
    # make some noise about any that didn't appear in the gtf index.
    my %MissedProts;
    my $numMissing = 0;
    my $protCount = keys %ProtNames;    
    print "Names of missing proteins:\n\n" if ($debug);
    foreach my $prot (sort keys %ProtNames) {
	if ($ProtNames{$prot}[0] != 2) {
	    print "  $prot\n" if ($debug);
	    $MissedProts{$prot} = 1;
	}
	$protCount++;
    }
    print "\n" if ($debug);
    
    
    # We'll alert the user in any case if we have proteins that aren't hitting.
    $numMissing = keys %MissedProts;
    if ($numMissing) {
	if ($verbose) { 
	    print "\n  WARNING: $numMissing (of $protCount) proteins from file '$isoformfilename' not found in index.\n";	
	}
    } elsif ($debug) {
	print "Corresponding genome index entries found for all sequences.";
    }

    
    # Return the indexed genes (and the list of blacklisted ones!)!
    return (\%GeneIndex,\%MissedProts);

}






#########################################################################
#
#  Function Name: AddExon
#
#  About: Once we've built an index on the genes we're interested in, we
#         need to rip each recorded exon off (as a start-index/end-index
#         pair) and actually compare that exon against the provided protein
#         sequence, using the translated-similarity-detection program
#         FindDiagonals.  AddExon funnels all the relevant information
#         towards the actual run of FindDiagonals (in RunFindDiagonals).
#         This primarily involves running esl-sfetch to get the relevant
#         part of the genome on-hand.
#
#
sub AddExon
{
    # The chromosome we're adding to
    my $chromosomeRef = shift;
    my $Chromosome    = $$chromosomeRef;

    # File names
    my $protfilename   = shift;
    my $genomefilename = shift;
    my $nuclfilename   = shift;

    # Information about the exon we intend to add
    my $gene           = shift;
    my $chromosomeName = shift;
    my $proteinLength  = shift;
    my $exonStart      = shift;
    my $exonEnd        = shift;

    # Are we debugging?
    my $debug = shift;
    
    # Whether or not we're considering the reverse complement will be
    # identifiable by whether or not we appended '[revcomp]' to the
    # chromosome name.
    my $revcomp = 0;
    if ($chromosomeName =~ /\[revcomp\]/) {
	$revcomp = 1;
	$chromosomeName =~ s/\[revcomp\]//;
	my $temp = $exonEnd;
	$exonEnd = $exonStart;
	$exonStart = $temp;
    }

    # The system call to generate a FASTA file using the given index info.
    my $eslsfetchCmd;
    $eslsfetchCmd = 'esl-sfetch -c '.$exonStart.'..'.$exonEnd;             # Range
    $eslsfetchCmd = $eslsfetchCmd.' -o '.$nuclfilename;                    # Output file
    $eslsfetchCmd = $eslsfetchCmd.' '.$genomefilename.' '.$chromosomeName; # Input file and sequence name

    # I'M MOVE ERROR OUTPUT TO GARBAGEVILLE
    #
    $eslsfetchCmd = $eslsfetchCmd." > /dev/null 2>&1";   # '-o' still prints some info, but we don't care.
    #my $garbage_file = $nuclfilename.'.GARBAGE_TIME.err.out';
    #$eslsfetchCmd = $eslsfetchCmd." > /dev/null 2>$garbage_file";   # '-o' still prints some info, but we don't care.
    
    # For now we bail altogether if this step goes wrong (for obv. reasons)
    if (system($eslsfetchCmd)) {
	# DEBUGGING
	#$eslsfetchCmd = 'esl-sfetch -c '.$exonStart.'..'.$exonEnd;             # Range
	#$eslsfetchCmd = $eslsfetchCmd.' -o '.$nuclfilename;                    # Output file
	#$eslsfetchCmd = $eslsfetchCmd.' '.$genomefilename.' '.$chromosomeName; # Input file and sequence name
	#system($eslsfetchCmd);
	# DEBUGGING
	die "\n\tERROR: Command '$eslsfetchCmd' failed (AE) - $gene\n\n"; }
    
    # Find all high-scoring diagonals for this exon (aligned with the given protein)
    RunFindDiagonals(\$Chromosome,$protfilename,$nuclfilename,$proteinLength,
		     $exonStart,$exonEnd,$revcomp,$debug);

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
    my $seqstatCmd = 'esl-seqstat -a '.$nuclfilename.' |';

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
#  Function Name: RunFindDiagonals
#
#  About: RunFindDiagonals performs the system call to search the protein
#         against some section of the genome (the particulars of this 'section'
#         corresponds to some entry in the .gtf index).  It then uses the
#         output of the program to add to our current archive of hits for
#         this particular gene and chromosome.         
#
sub RunFindDiagonals
{

    # Get reference to the 'DiagonalSet' object and lock it down
    my $chrRef     = shift;
    my $Chromosome = $$chrRef;

    # How many entries does this chromosome already have?
    my $diagNum = $Chromosome->{NumHits};

    # Get the names of the protein and nucleotide files
    my $protfile = shift;
    my $nuclfile = shift;

    # How long is our protein?
    my $proteinLength = shift;

    # Get the start and end points of the chromosome
    my $NuclStart = shift;
    my $NuclEnd   = shift;

    # Is this search on the chromosome's reverse complement?
    my $revcomp = shift;

    # Are we debugging?
    my $debug = shift;

    # PRINT THAT SUPER-FLY FindDiagonals.c DEBUGGING OUTPUT
    # If we want to peek at how FindDiagonals is doing, we can run it in debug mode,
    # SEPARATELY from running it in parsable-format
    #if ($debug && system("\.\/FindDiagonals $protfile $nuclfile \-debug")) { die "\n\tFindDiagonals failed\n\n"; }

    # Run FindDiagonals
    open(my $stdoutput,"\.\/FindDiagonals $protfile $nuclfile \|") || die "\n\tFindDiagonals failed\n\n";
    
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
	last if ($startoffsetLength !~ /[012]/);

	$startoffsetLength = int($startoffsetLength);

	if ($startoffsetLength) {
	    $line = <$stdoutput>;
	    $line =~ s/\n|\r//g;
	    $startoffsetChars = $line;
	} else {
	    $startoffsetChars = 0;
	}
	
	# 4.
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

	
	# 2.
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
	
	# 3.
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
			
			# Remove the chain of positions and offsets that got us here
			# from the blacklist
			#my @safePositions = split(/,/,$nextProtStart);
			#my @safeOffsets   = split(/,/,$nextOffset);
			#foreach $i (0..@safePositions-1) {
			#$RedHerrings{$safePositions[$i].'/'.$safeOffsets[$i]} = 0;
			#}
			
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
#  Function Name: FindPartialAlignments
#
#  About: FindPartialAlignments piggybacks on the diagonals that we've
#         already found by trimming off any offset characters from our
#         best extensions and then trying to stitch these together.
#         It's the last thing we do before turning to spaln and blat.
#
# >>> NORD: This function is currently deprecated.  We are no longer
#           interested in expending this effort, since our pipeline
#           of passing "Full" alignment misses to SPALN and BLAT is
#           currently fast enough.
#
#           I'm going to hold onto just in case it ever becomes useful
#           as a reference...
#
#sub FindPartialAlignments
#{
#    my ($i,$j,$k);
#
#    my $chromosomeRef = shift;
#    my $Chromosome = $$chromosomeRef;
#    
#    my $Protein    = shift;
#    my $protlength = length($Protein);
#
#    my $seqname = shift;
#
#    my $debug = shift;
#    
#    # The first contrast with FindFullAlignments is that we're
#    # interested in stitching together anything that's stitchable
#    # (which is to say, we're effectively doing a full traversal)
#    my $Hit = 0;
#
#    # Are we working on the reverse complement?
#    my $revcomp;
#    if ($Chromosome->{ChrName} =~ /\[revcomp\]/) { $revcomp = 1; }
#    else                                         { $revcomp = 0; }
#
#    # Because we don't want to end up with an explosion of redundant
#    # information (e.g., A->B->C  and  A->B  and  B->C), we don't
#    # allow any individual hit to be pulled from the Chromosome more
#    # than once.
#    my @SpentHits;
#    my %Pushed; # Recording full runs we've already encountered
#    foreach $i (0..$Chromosome->{NumHits}-1) { $SpentHits[$i] = 0; }
#
#    # We'll store the results of our search in a 'HitString' object
#    # (defined in the DiagonalSets module).  This can be thought of
#    # as a reduction of the Chromomsome to its longest implied runs.
#    my $ChrRedux = New DiagonalSet($Chromosome->{ChrName},1);
#    
#    # For this subroutine, we're doing the same sort of stack-
#    # based DFS as we used in FindFullAlignments, but instead of
#    # A. caring about whether the starting position is '0' in the
#    # Protein, or B. stopping if we can't get all the way to the
#    # end of the sequence, we're C. ID-ing ALL 'tendrils'
#    while ($Hit < $Chromosome->{NumHits}) {
#	
#	# Scan past anything we've already used
#	while ($Hit < $Chromosome->{NumHits} && $SpentHits[$Hit]) {
#	    $Hit++;	    
#	}	
#	last if ($Hit == $Chromosome->{NumHits});
#
#	# How is this looking?
#	my $hitLocated = ${$Chromosome->{IsTerminal}}[$Hit];
#	my $nuclstring = ${$Chromosome->{NuclStarts}}[$Hit].'-'.${$Chromosome->{NuclEnds}}[$Hit];
#	my $protstring = ${$Chromosome->{ProtStarts}}[$Hit].'-'.${$Chromosome->{ProtEnds}}[$Hit]; 
#	my $nuclstart  = ${$Chromosome->{StartOffsetLengths}}[$Hit];
#
#	# In the event that we have a hit starting at zero, we record it 'as is'
#	# Note that because we've failed to find a full extension, this will never
#	# also be a terminal hit.
#	if (${$Chromosome->{ProtStarts}}[$Hit] == 0 && !$Pushed{$nuclstring.$protstring}) {
#	    AddToHitString(\$ChrRedux,\$Chromosome,$nuclstring,$protstring,$Hit,$Hit);
#	}
#	
#	# Did we get lucky? Maybe...
#	if (!$hitLocated) { # MAYBE NOT!
#	    
#	    # Let's get stacking!
#	    push(my @hits,$Hit);
#	    push(my @nuclstrings,$nuclstring);
#	    push(my @protstrings,$protstring);
#	    push(my @nuclstarts,$nuclstart);
#	    push(my @starthits,$Hit);
#	    
#	    while (@hits) {
#		
#		# Pop the last hit position (and its string) off the stack
#		my $prevHit        = pop(@hits);
#		my $prevNuclString = pop(@nuclstrings);
#		my $prevProtString = pop(@protstrings);
#		my $originaloffset = pop(@nuclstarts);
#		my $originalhit    = pop(@starthits);
#
#		# Record that we've tried extending this hit
#		$SpentHits[$prevHit] = 1;
#		
#		# Setting requirements for the next hit position:
#		my $nextHit       = $prevHit+1;
#		my $nextProtStart = ${$Chromosome->{ProtEnds}}[$prevHit] + 1;
#		my $nextOffset    = (3 - ${$Chromosome->{EndOffsetLengths}}[$prevHit]) % 3;
#		$nextProtStart++ if ($nextOffset);
#		
#		# Special case consideration
#		if ($nextProtStart == $protlength && !$Pushed{$prevNuclString.$prevProtString}) {
#		    AddToHitString(\$ChrRedux,\$Chromosome,$prevNuclString,$prevProtString,$originalhit,$prevHit);
#		    $Pushed{$prevNuclString.$prevProtString} = 1;
#		    next;
#		}
#		
#		# Walk along looking for suitable extension points
#		my $hitlocated = 0;
#		while ($nextHit < $Chromosome->{NumHits}) {
#		    
#		    # Is this position viable?
#		    if (${$Chromosome->{ProtStarts}}[$nextHit] < $nextProtStart || $SpentHits[$nextHit]) {
#
#			$nextHit++;
#
#		    } elsif (${$Chromosome->{ProtStarts}}[$nextHit] == $nextProtStart) {
#			
#			# The final (BIG!) check for this possible next hit:
#			# Does it (1.) have the correct offset (2.) have a
#			# sensible DNA positioning (consistent with our direction)
#			if (${$Chromosome->{StartOffsetLengths}}[$nextHit] == $nextOffset
#			    && !(($revcomp && ${$Chromosome->{NuclStarts}}[$nextHit] >= ${$Chromosome->{NuclEnds}}[$prevHit]) 
#				 || (!$revcomp && ${$Chromosome->{NuclStarts}}[$nextHit] <= ${$Chromosome->{NuclEnds}}[$prevHit]))) {
#			    
#			    # WHAT A CINDERELLA STORY THIS IS!
#			    my $nextNuclString = $prevNuclString.','.${$Chromosome->{NuclStarts}}[$nextHit];
#			    my $nextProtString = $prevProtString.','.${$Chromosome->{ProtStarts}}[$nextHit];
#			    $nextProtString = $nextProtString.'-'.${$Chromosome->{ProtEnds}}[$nextHit];
#			    my $endoffset;
#			    
#			    # You'd better believe we located a hit!
#			    $hitlocated = 1;
#
#			    # Want to wrap this up to go?
#			    if (${$Chromosome->{IsTerminal}}[$nextHit]) {
#				
#				$nextNuclString = $nextNuclString.'-'.${$Chromosome->{NuclEnds}}[$nextHit];
#				
#				if (!$Pushed{$nextNuclString.$nextProtString}) {
#				    
#				    # Adding to our HitString object
#				    AddToHitString(\$ChrRedux,\$Chromosome,$nextNuclString,$nextProtString,$originalhit,$nextHit);
#				    $Pushed{$nextNuclString.$nextProtString} = 1;
#				    
#				}
#				
#				$SpentHits[$nextHit] = 1;
#				
#			    } else {
#
#				push(@hits,$nextHit);
#				$nextNuclString = $nextNuclString.'-'.${$Chromosome->{NuclEnds}}[$nextHit];
#				push(@nuclstrings,$nextNuclString);
#				push(@protstrings,$nextProtString);
#				push(@nuclstarts,$originaloffset);
#				push(@starthits,$originalhit);
#				
#			    }
#
#			    $nextHit++;
#			    
#			} else {
#			    
#			    $nextHit++;
#
#			}
#			
#		    } else {
#			
#			# The end of this run -- if we didn't push this string
#			# any farther, we record it as a standalone tendril.
#			if (!$hitlocated && !$Pushed{$prevNuclString.$prevProtString}) {
#			    AddToHitString(\$ChrRedux,\$Chromosome,$prevNuclString,$prevProtString,$originalhit,$prevHit);
#			    $Pushed{$prevNuclString.$prevProtString} = 1;
#			}
#			last;
#
#		    }
#
#		}
#		
#	    }
#		
#	} else {
#
#	    # This corresponds to finding a starter hit that happens to
#	    # be terminal (to be PC about what bigots would call "terminal hits")
#	    if (!$Pushed{$nuclstring.$protstring}) { # We might have seen this already if it starts at 0
#		AddToHitString(\$ChrRedux,\$Chromosome,$nuclstring,$protstring,$Hit,$Hit);
#		$Pushed{$nuclstring.$protstring} = 1;
#		$SpentHits[$Hit] = 1; # Not really necessary, but helps me sleep at night
#	    }
#	}
#	
#	$Hit++;
#	
#    }
#
#    # For debugging we might barf out what we're seeing.
#    #if ($debug) {
#    #print "\n\n";
#    #foreach $i (0..$ChrRedux->{NumHits}-1) {
#    #print "Start Offset: ${$ChrRedux->{StartOffsetLengths}}[$i]\n";
#    #print "Nucl Tendril: ${$ChrRedux->{NuclStrings}}[$i]\n";
#    #print "Prot Tendril: ${$ChrRedux->{ProtStrings}}[$i]\n";
#    #print "  End Offset: ${$ChrRedux->{EndOffsetLengths}}[$i]\n";
#    #print "------------+\n";
#    #}
#    #print "\n\n";
#    #}
#
#    # Finally, because we'll be running translated Smith-Waterman looking
#    # for hits to particular protein positions, we trim off all the
#    # offset info. for the reduced hits.
#    my ($offset,$old_pos,$new_pos);
#    foreach $i (0..$ChrRedux->{NumHits}-1) {
#
#	# Start offset trimming
#	$offset = ${$ChrRedux->{StartOffsetLengths}}[$i];
#	if ($offset) {
#
#	    $old_pos = ${$ChrRedux->{NuclStarts}}[$i];
#
#	    if ($revcomp) { $new_pos = $old_pos-$offset; }
#	    else          { $new_pos = $old_pos+$offset; }
#
#	    ${$ChrRedux->{NuclStrings}}[$i] =~ s/^\d+\-/\-/;
#	    ${$ChrRedux->{NuclStrings}}[$i] =  $new_pos.${$ChrRedux->{NuclStrings}}[$i];
#	    ${$ChrRedux->{NuclStarts}}[$i]  =  $new_pos; 
#
#	    ${$ChrRedux->{StartOffsetLengths}}[$i] = 0;
#	    ${$ChrRedux->{StartOffsets}}[$i]       = 0;
#
#	}
#
#
#	# End offset trimming
#	$offset = ${$ChrRedux->{EndOffsetLengths}}[$i];
#	if ($offset) {
#
#	    $old_pos = ${$ChrRedux->{NuclEnds}}[$i];
#
#	    if ($revcomp) { $new_pos = $old_pos+$offset; }
#	    else          { $new_pos = $old_pos-$offset; }
#
#	    ${$ChrRedux->{NuclStrings}}[$i] =~ s/\-\d+$/\-/;
#	    ${$ChrRedux->{NuclStrings}}[$i] =  ${$ChrRedux->{NuclStrings}}[$i].$new_pos;
#	    ${$ChrRedux->{NuclEnds}}[$i]    =  $new_pos; 
#
#	    ${$ChrRedux->{EndOffsetLengths}}[$i] = 0;
#	    ${$ChrRedux->{EndOffsets}}[$i]       = 0;
#
#	}
#
#    }
#
#    return (\$ChrRedux);
#    
#}





#########################################################################
#
#  Function Name: GetBestExons
#
#  About: GetBestExons pulls out the best hits from FindDiagonals
#         that are (1) internally consistent as possibly being part of
#         the final hit and (2) at least length 4 (unlikely to occur
#         purely by chance).  These exons are pulled in order of length.
#
sub GetBestExons
{

    my ($i,$j,$k);
    my $ChromosomeRef = shift;
    my $Chromosome    = $$ChromosomeRef;

    # Sort the hits by their length. We're primarily interested in the
    # longest hits, since long hits are unlikely to be chance alignments
    $Chromosome->SortHitsByField('Length',0);

    # Setting a lower bound on the length limit for a hit that we're
    # willing to marry ourselves to
    my $lowerlimit = 4;

    #foreach $i (0..$Chromosome->{NumHits}-1) {
    #print "${$Chromosome->{ProtStrings}}[$i]: ${$Chromosome->{NuclStarts}}[$i],${$Chromosome->{NuclEnds}}[$i]\n";
    #}

    my $ProtLength = shift;
    
    my $ChrRedux = New DiagonalSet($Chromosome->{ChrName});

    # Because the chromosome is currently sorted by descending order of
    # hit length, we can do a linear walk through our hits pulling out
    # all hits that fill in some part of the remaining gap
    my @AcceptableEnds;
    foreach $i (0..$ProtLength-1) {
	$AcceptableEnds[$i]   = $ProtLength-1;
    }
    
    my ($nextStart,$nextEnd);
    foreach $i (0..$Chromosome->{NumHits}-1) {

	# A really short hit (length 3 or less) could throw off 
	# the acceptable nucleotide range during our hits.
	last if (${$Chromosome->{ProtEnds}}[$i]-${$Chromosome->{ProtStarts}}[$i] < $lowerlimit);
	
	# Is this position still available?  Can we go as far as we'd need?
	if ($AcceptableEnds[${$Chromosome->{ProtStarts}}[$i]] >= ${$Chromosome->{ProtEnds}}[$i]) {
		    
	    CopyHit(\$ChrRedux,\$Chromosome,$i);

	    # Adjusting our list of available positions.
	    $j = ${$Chromosome->{ProtStarts}}[$i];
	    $j = 0 if ($j < 0);
	    while ($j < $ProtLength && $j <= ${$Chromosome->{ProtEnds}}[$i]) { 
		$AcceptableEnds[$j]   = 0;
		$j++;
	    }

	}

    }

    #foreach $i (0..$ChrRedux->{NumHits}-1) {
    #print "-------------.\n";
    #print " PROT STRING : ${$ChrRedux->{ProtStrings}}[$i]\n";
    #print " NUCL STRING : ${$ChrRedux->{NuclStrings}}[$i]\n";
    #print "             ;\n";
    #}
    #print "\n\n";
    
    # Typically, having our chromosome organized by order of amino
    # acid position is the most helpful way to structure it, so we
    # give it a nice quick sort before chucking it back up the line.
    $ChrRedux->SortHitsByField('ProtStarts',1);
    
    return \$ChrRedux;
}





#########################################################################
#
#  Function Name: AddToHitString
#
#  About: AddToHitString is effectively comparable to the 'AddExon'
#         subroutine, but is used to copy data from a non-HitString
#         DiagonalSet to a HitString DiagonalSet.  It copies start
#         points from one entry and terminal points from another to
#         give a picture of where the 'runs' are.
#
sub AddToHitString
{
    my $chrReduxRef = shift;
    my $ChrRedux    = $$chrReduxRef;

    my $chrRef     = shift;
    my $Chromosome = $$chrRef;

    my $NuclString = shift;

    my $ProtString = shift;

    my $startHit = shift;
    my $finalHit = shift;

    my $hitNum = $ChrRedux->{NumHits};
    $ChrRedux->{NumHits}++;

    my $isTerm = ${$Chromosome->{IsTerminal}}[$finalHit];
    ${$ChrRedux->{IsTerminal}}[$hitNum] = $isTerm;
    $ChrRedux->{NumTerminals}          += $isTerm;

    ${$ChrRedux->{NuclStrings}}[$hitNum] = $NuclString;
    ${$ChrRedux->{NuclStarts}}[$hitNum]  = ${$Chromosome->{NuclStarts}}[$startHit];
    ${$ChrRedux->{NuclEnds}}[$hitNum]    = ${$Chromosome->{NuclEnds}}[$finalHit];

    ${$ChrRedux->{ProtStrings}}[$hitNum] = $ProtString;
    ${$ChrRedux->{ProtStarts}}[$hitNum]  = ${$Chromosome->{ProtStarts}}[$startHit];
    ${$ChrRedux->{ProtEnds}}[$hitNum]    = ${$Chromosome->{ProtEnds}}[$finalHit];

    ${$ChrRedux->{StartOffsetLengths}}[$hitNum] = ${$Chromosome->{StartOffsetLengths}}[$startHit];
    ${$ChrRedux->{StartOffsets}}[$hitNum]       = ${$Chromosome->{StartOffsets}}[$startHit];
    ${$ChrRedux->{EndOffsetLengths}}[$hitNum]   = ${$Chromosome->{EndOffsetLengths}}[$finalHit];
    ${$ChrRedux->{EndOffsets}}[$hitNum]         = ${$Chromosome->{EndOffsets}}[$finalHit];

    # WRONG, BUT FIELD NEEDED VALUE
    ${$ChrRedux->{Scores}}[$hitNum] = ${$Chromosome->{Scores}}[$finalHit];

}






#########################################################################
#
#  Function Name: CopyHit
#
#  About: Unlike AddToHitString, which combines two distinct hits,
#         CopyHit is a straightforward copy of one DiagonalSet.
#
sub CopyHit
{
    my $chrReduxRef = shift;
    my $ChrRedux    = $$chrReduxRef;

    my $chrRef     = shift;
    my $Chromosome = $$chrRef;

    my $grabFrom = shift;
    my $placeAt  = $ChrRedux->{NumHits};
    $ChrRedux->{NumHits}++;

    ${$ChrRedux->{ProtStrings}}[$placeAt] = ${$Chromosome->{ProtStrings}}[$grabFrom];
    ${$ChrRedux->{ProtStarts}}[$placeAt]  = ${$Chromosome->{ProtStarts}}[$grabFrom];
    ${$ChrRedux->{ProtEnds}}[$placeAt]    = ${$Chromosome->{ProtEnds}}[$grabFrom];

    ${$ChrRedux->{StartOffsetLengths}}[$placeAt] = ${$Chromosome->{StartOffsetLengths}}[$grabFrom];
    ${$ChrRedux->{StartOffsets}}[$placeAt]       = ${$Chromosome->{StartOffsets}}[$grabFrom];
    ${$ChrRedux->{EndOffsetLengths}}[$placeAt]   = ${$Chromosome->{EndOffsetLengths}}[$grabFrom];
    ${$ChrRedux->{EndOffsets}}[$placeAt]         = ${$Chromosome->{EndOffsets}}[$grabFrom];

    ${$ChrRedux->{NuclStrings}}[$placeAt] = ${$Chromosome->{NuclStrings}}[$grabFrom];
    ${$ChrRedux->{NuclStarts}}[$placeAt]  = ${$Chromosome->{NuclStarts}}[$grabFrom];    
    ${$ChrRedux->{NuclEnds}}[$placeAt]    = ${$Chromosome->{NuclEnds}}[$grabFrom];

    ${$ChrRedux->{Scores}}[$placeAt]     = ${$Chromosome->{Scores}}[$grabFrom];
    ${$ChrRedux->{StopCodon}}[$placeAt]  = ${$Chromosome->{StopCodon}}[$grabFrom];
    ${$ChrRedux->{IsTerminal}}[$placeAt] = ${$Chromosome->{IsTerminal}}[$grabFrom];

    $ChrRedux->{NumTerminals} += ${$ChrRedux->{IsTerminal}}[$placeAt];
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

    # Find min. and max. nucleotide positions based on exons
    my $minNucl = MIN(${$Chromosome->{NuclStarts}}[0],${$Chromosome->{NuclEnds}}[0]);
    my $maxNucl = MAX(${$Chromosome->{NuclStarts}}[0],${$Chromosome->{NuclEnds}}[0]);
    foreach my $exon (1..$Chromosome->{NumHits}-1) {
	$minNucl = MIN($minNucl,MIN(${$Chromosome->{NuclStarts}}[$exon],${$Chromosome->{NuclEnds}}[$exon]));
	$maxNucl = MAX($maxNucl,MAX(${$Chromosome->{NuclStarts}}[$exon],${$Chromosome->{NuclEnds}}[$exon]));
    }

    # Pull in a little extra sequence, just to be safe
    $minNucl = MAX($minNucl-1000,1);
    $maxNucl = MIN($maxNucl+1000,$ChrLengths{$chrname});

    # In case we've been given an inordinate amount of sequence to investigate,
    # we'll jump ship
    return (0,0) if ($maxNucl - $minNucl > 1000000);

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
					       $seqname,$ProtLength,$SpalnMisses,$spalner);


    # Return your bliss if we got a clean hit, bail if things look bad
    if    ($hitscore)               { return ($hitscore,$hitline); }
    elsif (!$hitscore && !$hitline) { return (0,0);                }

    
    #
    # Trying a bigger chunk of sequence
    #

    
    # If ParseSPALNOutput returns (0,1), it means that we hit at >90%, so
    # how about we expand our borders and see if we can't get up to 97%?

    # Pull in a LOT of extra sequence
    $minNucl = MAX($minNucl-400000,1);
    $maxNucl = MIN($maxNucl+400000,$ChrLengths{$chrname});

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
					    $seqname,$ProtLength,$SpalnMisses,$spalner);


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

    # Are we just looking at SPALN output?
    my $spalner = shift;

    my $protfilename = 'Quilter.BAS.prot.fa';
    my $nuclfilename = 'Quilter.BAS.nucl.fa';
    
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
		else                      { @{$Chromosomes{$seqID}} = [$chrname];   }
	    
		if ($rev == 0) { #################################### Forward Hit
		    if ($NuclStarts{$seqID}) {
			push(@{$Revcomp{$seqID}},0);
			push(@{$NuclStarts{$seqID}},$nuclstart);
			push(@{$NuclEnds{$seqID}},$nuclend);
		    } else {
			@{$Revcomp{$seqID}}    = [0];
			@{$NuclStarts{$seqID}} = [$nuclstart];
			@{$NuclEnds{$seqID}}   = [$nuclend];
		    }
		    
		} else { ############################################ Reverse Hit
		    if ($NuclStarts{$seqID}) { 
			push(@{$Revcomp{$seqID}},1);
			push(@{$NuclStarts{$seqID}},$nuclend);
			push(@{$NuclEnds{$seqID}},$nuclstart);	
		    } else {
			@{$Revcomp{$seqID}}    = [1];
			@{$NuclStarts{$seqID}} = [$nuclend];
			@{$NuclEnds{$seqID}}   = [$nuclstart];
		    }
		    
		}
		
		# Record the start and end positions in the protein
		if ($ProtStarts{$seqID}) {
		    push(@{$ProtStarts{$seqID}},$protstart);
		    push(@{$ProtEnds{$seqID}},$protend);
		} else {
		    @{$ProtStarts{$seqID}} = [$protstart];
		    @{$ProtEnds{$seqID}}   = [$protend];
		}
		

		# Get the score
		if ($line =~ /(\S+)\s+(\S+)\s*$/) {
		    
		    my $e_value = $1;
		    my $score   = $2;

		    if ($Scores{$seqID}) { 
			push(@{$Scores{$seqID}},$score);
		    } else { 
			@{$Scores{$seqID}} = [$score];
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
		    } elsif ($e_value eq "0.0") {
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

    # Iterate over sequences, SPALNing 'em up and down
    #
    # UPDATE -- Increasingly wide ranges as hits fail to connect
    #
    foreach my $seqID (keys %BlatNameIndex) {

	my $seqname = $BlatNameIndex{$seqID};

	# Did we find any suitable hits?
	if ($ConfirmedHits{$seqname}) {

	    # We're going to need that sequence
	    my $eslsfetchCmd = 'esl-sfetch  '.$proteinfilename." '".$seqname."' > ".$protfilename;
	    if (system($eslsfetchCmd)) { die "\n  ERROR:  esl-sfetch command '$eslsfetchCmd' failed\n\n"; }
	    
	    $seqname =~ /GN:([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|(\S*\|)?([^\|\s]+)\s*$/;
	    my $orig_gene = $1;
	    my $isoform   = $2;
	    my $species   = $3;
	    my $iso_id    = $4;
	    my $group_id  = $6;

	    my $seqlength = $SeqLengths{$seqname};

	    # Record this sequence (and number of significant e-values) if
	    # it meets our threshold (100)
	    if ($TEtracker{$seqname} && $TEtracker{$seqname} >= 100) {
		print $TEfile "$seqname\n";
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
		    if ($revcomp) { $eslsfetchCmd = 'esl-sfetch -c '.$maxNucl.'..'.$minNucl; }
		    else          { $eslsfetchCmd = 'esl-sfetch -c '.$minNucl.'..'.$maxNucl; }	    
		    $eslsfetchCmd = $eslsfetchCmd.' -o '.$nuclfilename;
		    $eslsfetchCmd = $eslsfetchCmd.' '.$genomefilename.' '.$original_chr;
		    $eslsfetchCmd = $eslsfetchCmd." > /dev/null 2>&1";
		    if (system($eslsfetchCmd)) { die "\n\tERROR: Command '$eslsfetchCmd' failed (BAS)\n\n"; } 

		    # Assemble the sytem call to run SPALN
		    #my $spalnCmd = 'spaln -Q3 -O1 -S3 '.$nuclfilename.' '.$protfilename.' |';
		    my $spalnCmd = 'spaln -Q3 -O1 -S3 -ya3 '.$nuclfilename.' '.$protfilename;
		    $spalnCmd = $spalnCmd.' 2>/dev/null';# if (!$verbose);
		    $spalnCmd = $spalnCmd.' |';	
		    
		    print $CmdLog "> (".localtime()." $spalnCmd\n\n";
		    
		    # Make sure we're recording intelligently...
		    if ($revcomp) { $chrname = $chrname.'[revcomp]'; $revcomp = 3; } #DEBUGGING

		    # Parse SPALN's output
		    my ($nextScore,$nextLine) = ParseSPALNOutput($spalnCmd,$revcomp,$minNucl-1,$highscore,
								 $chrname,$seqname,$seqlength,$CmdLog,$spalner);

		    
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
		$hitline =~ s/Isoform ID \: GN\:/Isoform ID \: /;
		$hitline =~ s/Method     \: SPALN/Method     \: SPALN\+BLAT/;
		print $HitFile  "$hitline";
	    } else { 
		print $MissFile "$seqname\n"; 
	    }
	    
	} else {

	    # No hits from blat for this sequence
	    print $MissFile "$seqname\n";

	}
    }

    # Clean up after yourself, man!
    #
    if (-e $nuclfilename) { system("rm $nuclfilename"); }
    if (-e $protfilename) { system("rm $protfilename"); }

}






#########################################################################
#
#  Function Name: ParseSPALNOutput
#
#  About: The name says it all.  Run a SPALN command and make sense of
#         its wisdom.
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

    my $prot_len = shift;
    
    my $SpalnMisses = shift;

    # NOTE: When we're just looking for SPALN output we can't just
    #       run the command and jump ship, due to how the various
    #       functions depend on score comparison (well, BLAT).
    #       For this reason, we'll set the 'hitstring' to be the full
    #       output AFTER we get the score.
    #
    my $spalner = shift;
    
    my ($line,$hitstring);

    # DEBUGGING
    #my $printspalnCmd = $spalnCmd;
    #$printspalnCmd =~ s/\|$//;
    #system($printspalnCmd);
    # DEBUGGING
	
    open(my $stdout,$spalnCmd) || die "\n  ERROR: Spaln command '$spalnCmd' failed\n\n";


    # Look for where SPALN has called the start and end of the region.
    # Note that these will always be formatted as "low-high", even
    # if we're using the reverse complement
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
    # Note that these are (still) always low-high
    $line =~ /\S+\/(\d+)\-(\d+)/;
    my $start_pos = $1;
    my $end_pos   = $2;


    # Go to the section that lays out the score and identity
    my ($hitscore,$percent,$true_num_chars);
    while ($line = <$stdout>) {

	if ($line =~ /Score \= (\d+)/) {
	    
	    # SPALN's score
	    #
	    $hitscore = $1;

	    # Grab the percent identity
	    #
	    $line =~ /\((\d+)\.\d+ \%\)[\s\n\r]+$/;
	    $percent  = $1;

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

    
    # If we've hit the end of the file or have worse than 97% identity
    # jump ship.
    if (eof($stdout)) {
	close($stdout); 
	return(0,0); 
    } elsif ($percent < 97 || $true_num_chars < $prot_len) {
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
    my @FullNuclSeq;
    my @FullProtSeq;
    my @AAPositions;
    my @NuclPositions;    
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
	# translation of the nucleotide sequence (and we'll
	# want to skip it).
	#
	next if ($line !~ /\| \S+\/\d+\-\d+/);

	my $nn_line = $line;

	# EXCELLENT!  This is the next line of DNA characters
	#
	my @NextNucls = split('',$line);


	# The final check we need to do is to add/subtract
	# the initial 'jump' into the hit (but only if this
	# is the very beginning)
	#
	if ($first_jump == 0) {

	    $line =~ /\s+(\d+)/;
	    $first_jump = $1;

	    # Special catch for BLAT parsing
	    if ($revcomp==3) { $current_pos  = ($start_pos+1)-$first_jump; }
	    else             { $current_pos += $first_jump;                }

	    $first_jump = 1;

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


    # Check if we failed to get the right number
    # of amino acids.
    #
    return(0,0) if ($num_aas != $prot_len);



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
	my @MicroRun = split('-',$micro_run_str);
	$micro_start = $MicroRun[0];
	$micro_end   = $MicroRun[1]; # One past actual end, so we can do strict '<'


	# Identify the amino acid sequence that constitutes
	# this micro exon (or chain of micro exons)
	#
	my @MicroAAs;
	for ($j = $micro_start; $j < $micro_end; $j++) {
	    push(@MicroAAs,$FullProtSeq[$j]);
	}


	# Figure out how long the micro exon is
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
	    my $rear_codon_center = $NuclPositions[$AAPositions[$micro_start-1]];	    
	    my $nuclseq_pos = $AAPositions[$micro_start-1]+3;	    
	    while ($rear_ext_len < $micro_len) {
		
		# Make sure we don't overstep
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
		if (uc($FullProtSeq[$micro_start+($rear_ext_len+1)]) eq uc(TranslateCodon(\@NextCodon))) {
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
	    for ($i = 0; $i < $micro_len; $i++) { $AAPositions[$micro_start+$i] = $RearExt[$i]; }
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
	    
	    for ($i = 0; $i < $micro_len; $i++) { $AAPositions[$micro_start+$i] = $FwdExt[$i]; }
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
	$j = $micro_end - $fwd_ext_len;
	while ($i < $micro_len) {
	    $AAPositions[$j] = $FwdExt[($micro_len-1)-$i];
	    $i++;
	    $j++;
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




sub GenSortIndex
{
    my ($i,$j,$k);

    my $targetRef = shift;
    my @Target    = @{$targetRef};

    my @Index;
    my $length = 0;
    while ($length < @Target) {
	$Index[$length] = $length;
	$length++;
    }

    # Mergesort (pulled over from DiagonalSets.pm)
    # Run mergesort on the set, using ProtStarts to determine order.
    my ($start,$end,$placer);
    my @temp;
    my $groupSize = 1;
    while ($groupSize < $length) {
        foreach my $groupID (0..POSIX::ceil(($length/(2*$groupSize)))-1) {
            $start = $groupID * (2 * $groupSize);
            $end   = $start + (2 * $groupSize) - 1;
            $end   = $length-1 if ($end >= $length);
            foreach $i (0..$end-$start) {
                $temp[$i] = $Index[$start+$i];
            }
            $i = 0;
            $j = $groupSize;
            $placer = $start;
            while ($i < $groupSize && $j <= $end-$start) {
                if ($Target[$temp[$i]] < $Target[$temp[$j]]) {
                    $Index[$placer] = $temp[$i];
                    $i++;
                } else {
                    $Index[$placer] = $temp[$j];
                    $j++;
                }
                $placer++;
            }
            if ($i < $groupSize) {
                while ($i < $groupSize) {
                    $Index[$placer] = $temp[$i];
                    $i++;
                    $placer++;
                }
            } else {
                while ($j <= $end-$start) {
                    $Index[$placer] = $temp[$j];
                    $j++;
                    $placer++;
                }
            }
        }
        $groupSize *= 2;
    }

    # Great! Now pass that bad-boy back!
    return \@Index;

}




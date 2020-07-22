#!/usr/bin/env perl
#
# MIRAGE.pl - Multiple-sequence Isoform Alignment Tool Guided by Exon Boundaries
#           - Alex Nord
#           - 2016
#
# >>> M2 (2020) : Look forward to seeing big changes 'round here!
#
# ABOUT: This is a top-level script used to perform isoform alignment across
#        species.  It requires a fasta-formatted database of proteins and a
#        file containing information about where genomes and .gtf indices for each
#        species that is being considered can be found.  The program outputs
#        gene-wise MSAs to a directory (by default: MirageResults/FinalMSAs).
#
#
use warnings;
use strict;
use Getopt::Long;
use Time::HiRes;
use POSIX;
use Cwd;

# YUCK
sub GetThisDir { my $lib = $0; $lib =~ s/\/Mirage2.pl$//; return $lib; }
use lib GetThisDir();
use BureaucracyMirage;
use DisplayProgress;

sub PrintUsage;
sub DetailedUsage;
sub PrintVersion;
sub CheckInstall;
sub CheckSourceFiles;
sub ParseArgs;
sub VerifiedClean;
sub ParseSpeciesGuide;
sub ConfirmSSI;
sub GenerateSpeciesDBs;
sub CoverMinorSpecies;
sub ParseSeqNameAsMirage;
sub ParseSeqNameAsUniProt;
sub AttachSeqToMSA;


# VERSION INFO
my $mirage_version = 'Version 2.0.0a';


##################
##################
####          ####
####  SCRIPT  ####
####          ####
##################
##################



# In case the minimum number of arguments hasn't been provided,
# assist the user
if (@ARGV == 0) { PrintUsage(); }

# Generics
my ($i,$j,$k);


# Figure out what the location of the Mirage src directory is
my $location = $0;
$location =~ s/Mirage2\.pl$//;


# We're going to need these friends
my $eslsfetch  = $location.'../inc/easel/miniapps/esl-sfetch';
my $eslseqstat = $location.'../inc/easel/miniapps/esl-seqstat';


# Where we'll be storing timing information
my $StartTime = [Time::HiRes::gettimeofday()];
my $IntervalStart;
my $IntervalEnd;
my $FinalTime;
my @QuilterTimeStats;
my @MultiMSATimeStats;
my $MultiSeqNWTime;
my $AvgNWTime;
my $TotalRuntime;
my $defaultcpus = 2;


# Read in user arguments
my $optsRef = ParseArgs();
my %Options = %{$optsRef};


# Change options to more intelligible variables
my $ProteinDB    = ConfirmFile($ARGV[0]);
my $SpeciesGuide = ConfirmFile($ARGV[1]);
my $ResultsDir   = $Options{outdirname};
my $verbose      = $Options{verbose};
my $num_cpus     = $Options{cpus};
my $timed        = $Options{time};
my $stack_arfs   = $Options{stackarfs};
my $forcecompile = $Options{forcecompile}; # Hidden
my $cleanMSA     = $Options{cleanmsa};     # Hidden
my $just_spaln   = $Options{justspaln};    # Hidden
my $track_spaln  = $Options{trackspaln};   # Hidden


# We've officially started Mirage-ery!
DispProgMirage('init');

# Verify that we have all the files we need on-hand
#CheckSourceFiles($forcecompile);


# Convert species guide into three arrays.
# We'll ignore any species that don't show up
# in our protein database.
my ($speciesref,$gtfsref,$genomesref) = ParseSpeciesGuide($SpeciesGuide,$ProteinDB);
my @Species    = @{$speciesref};
my $numSpecies = scalar(@Species);
my @GTFs       = @{$gtfsref};
my @Genomes    = @{$genomesref};


# Make a directory to store results
$ResultsDir = CreateDirectory($ResultsDir);
print "\n  Results Directory:  $ResultsDir\n\n" if ($verbose);


# Stick in a directory to record progress information as we're running the
# various tools that we're so happy to have this opportunity to be running!
my $progress_dirname = CreateDirectory($ResultsDir.'.progress');

# We'll also create a directory to house any sequences belonging to species that
# we don't have genomes for.  'Misc' is safe to use because all actual species are
# lower-cased in our list.
push(@Species,'Misc');

# Make species-specific results folders
my %SpeciesDir;
my %SpeciesSeqDir;
my $AllSpeciesDir = CreateDirectory($ResultsDir.'SpeciesMSAs');
foreach my $species (@Species) {
    $SpeciesDir{$species}    = CreateDirectory($AllSpeciesDir.$species);
    $SpeciesSeqDir{$species} = CreateDirectory($SpeciesDir{$species}.'seqs');
}

# Create a temp directory in the results directory where we'll hide all of our secrets.
#my $tempdirname = $ResultsDir.'temp/';
#if (-d($tempdirname)) { system("rm -rf $tempdirname"); }
#if (system("mkdir $tempdirname")) { die "\n  ERROR:  Failed to generate temporary directory '$tempdirname'\n\n"; }
#
# EDIT: Each process is liable to have its own tempdir location, so we accomodate them

# Make species-specific databases, organized by gene family
#
# M2: We're going to want to be able to parse straight-from-SwissProt/TrEMBL seqnames.
#
#     We're also switching all names to a simple numeric value.  Everything else is
#     already implicit (i.e., species and gene names).
#
# TODO, MAYBE: Allow for users to provide a list of genes to try searching a family
#              against (e.g., "nord_pseudo_fam1: obscn titin fgfr2") -- synonyms

DispProgMirage('db-speciation');
my ($originalseqnames_ref) = GenerateSpeciesDBs($ProteinDB,$num_cpus,\%SpeciesSeqDir);
my @OriginalSeqNames = @{$originalseqnames_ref};

# DEBUGGING
# I'm going to print off all the original sequence names to a file so
# if things go wrong during testing I can look up the particular sequences
# and zero in on them.
my $debug_f = OpenOutputFile($ResultsDir.'seq-name-guide');
for (my $i=0; $i<scalar(@OriginalSeqNames); $i++) {
    print $debug_f "$i: $OriginalSeqNames[$i]\n";
}
close($debug_f);
# DEBUGGING

my $MinorSpeciesDBName;
my @SpeciesDBNames;

# Run Quilter
my %MultiMSADir;
my %GroupsBySpecies;
for ($i=0; $i<$numSpecies; $i++) {

    # DEBUGGING
    next if ($Species[$i] ne 'human');
    # DEBUGGING

    # Start a timer for Quilter
    $IntervalStart = [Time::HiRes::gettimeofday()];

    # Start assembling that command!
    my $QuilterCmd = 'perl '.$location.'Quilter2.pl';

    # KEY 1: Protein database (implicit in species directory, under 'seqs/')
    $QuilterCmd = $QuilterCmd.' '.$SpeciesDir{$Species[$i]};

    # KEY 2: Genome
    $QuilterCmd = $QuilterCmd.' '.$Genomes[$i];

    # KEY 3: GTF index (special case: - for "Just use BLAT")
    $QuilterCmd = $QuilterCmd.' '.$GTFs[$i];
	
    # Do we want a whole bunch of stuff spat at us?
    $QuilterCmd = $QuilterCmd.' --v' if ($verbose);

    # Are we timing?
    $QuilterCmd = $QuilterCmd.' --time' if ($timed);
    
    # Display the call if we're being loud
    if ($verbose) {
	print "\n  Aligning to the $Species[$i] genome using Quilter.pl\n";
	print "  $QuilterCmd\n";
    }

    # Make that daaaaaang call.
    RunSystemCommand($QuilterCmd);

    # DEBUGGING
    die;
    # DEBUGGING

    # Figure out how much time that took and record it
    $QuilterTimeStats[$i] = Time::HiRes::tv_interval($IntervalStart);


    #################################################################### COUNTING DENSITY OF MICRO-EXONS
    #
    # We're interested in knowing how many of the alignments to
    # the genome that quilter came up with have dense pockets
    # of really small exons
    #
    my $max_splice_sites   = 3;
    my $observed_max       = 0;
    my $splice_window_size = 7;
    my $num_reported_seqs  = 0;
    my $total_num_seqs     = 0;
    my $current_seq_name;
    
    my $WarningFileName = $SpeciesDir{$Species[$i]}.'/MicroExonWarnings.out';
    
    open(my $WarningFile,'>',$WarningFileName) || die "\n  ERROR:  Failed to open '$WarningFileName'\n\n";
    open(my $CodonFile,'<',$SpeciesDir{$Species[$i]}.'/Hits.Quilter.out') || die "\n  NASTEEEEY\n\n";
    
    while (my $line = <$CodonFile>) {

	$line =~ s/\n|\r//g;
	if ($line =~ /^Match Pos\.s\:\s+(\S+)\s?$/) {

	    my @Codons = split(',',$1);
	    my @SpliceDensity = ();

	    for (my $check = 0; $check < @Codons; $check++) {

		if ($Codons[$check] eq '*') {

		    for (my $fwd = 0; $check+$fwd < @Codons && $fwd <= $splice_window_size; $fwd++) {

			if ($SpliceDensity[$check+$fwd]) { $SpliceDensity[$check+$fwd]++; }
			else                             { $SpliceDensity[$check+$fwd]=1; }
			
			if ($SpliceDensity[$check+$fwd]==$max_splice_sites) {
			    print $WarningFile "$current_seq_name\n";
			    $observed_max = 1;
			    last;
			}
		    }
		}

		if ($observed_max) {
		    $num_reported_seqs++;
		    $observed_max = 0;
		    last;
		}
		
	    }

	    $total_num_seqs++;
	    
	} elsif ($line =~ /^Isoform\s+ID\s+\:\s+(\S+)\s?$/) {
	    
	    $current_seq_name = $1;

	}
	
    }

    if ($num_reported_seqs) { print $WarningFile "\n  $num_reported_seqs of $total_num_seqs\n\n"; }
    
    close($CodonFile);
    close($WarningFile);

    if (!$num_reported_seqs && -e $WarningFileName) { system("rm $WarningFileName"); }
    #
    #
    #################################################################### END COUNTING STUFF

        
    # Now we want to check how long it takes to run MultiMSA
    $IntervalStart = [Time::HiRes::gettimeofday()];

    # Grab a hold of your toys (confirming that they're where they should be)
    my $HitFileName  = $SpeciesDir{$Species[$i]}.'/Hits.Quilter.out';
    my $MissFileName = $SpeciesDir{$Species[$i]}.'/Misses.Quilter.out';
    my $NearFileName = $SpeciesDir{$Species[$i]}.'/NearHits.Quilter.out';
    #if (!(-e $HitFile)) { die "\n  Failed to locate Quilter output '$HitFile'\n\n"; }

    # Name the directory where we'll be putting the output from MultiMSA.pl
    $MultiMSADir{$Species[$i]} = $SpeciesDir{$Species[$i]}.'/Alignments';

    # Prepare to generate a file containing tallies for each gene with
    # disagreeing positions.
    my $DisFilename = $SpeciesDir{$Species[$i]}.'/CandidateARFs.txt';
    
    # Assembling the MultiMSA command.
    my $MultiMSACmd;
    $MultiMSACmd = 'perl '.$location.'MultiMSA.pl '; # base call
    $MultiMSACmd = $MultiMSACmd.$HitFileName.' ';    # Quilter output
    $MultiMSACmd = $MultiMSACmd.$SpeciesDBNames[$i]; # isoform file
    
    # Guiding output towards our desired directory
    $MultiMSACmd = $MultiMSACmd.' -folder '.$MultiMSADir{$Species[$i]};

    # Tell MultiMSA how many CPUs we want it to work with
    $MultiMSACmd = $MultiMSACmd.' -n '.$num_cpus;

    # Give it reference to the quilter output folder, so we can scoop
    # up any of those tasty multi-chromosomal gene families.
    $MultiMSACmd = $MultiMSACmd.' -misses '.$MissFileName;

    # If we want to stack ARFs in the output alignment files, now is
    # the time to communicate our desires
    $MultiMSACmd = $MultiMSACmd.' -stack-arfs' if ($stack_arfs);
    
    # Display the call and/or make it
    if ($verbose) {
	print "  Generating MSAs for all identified $Species[$i] genes\n";
	print "  $MultiMSACmd\n";
    }

    if (system($MultiMSACmd)) {
	# HAMSTERS!
	die "\n  *  ERROR:  MultiMSA.pl failed during execution  *\n\n";
    }


    # Grab all group filenames.  Note that these should be identical across
    # species.  We also count how many isoforms each group had within the species.
    # Because we might be editing the contents of the directory, we begin by
    # storing the filenames in an array, rather than doing it on-the-fly.
    opendir(my $dir,$MultiMSADir{$Species[$i]}) || die "\n  Failed to open directory '$MultiMSADir{$Species[$i]}'\n\n";
    open(my $DisFile,'>',$DisFilename);
    my @DirContents;
    while (my $filename = readdir($dir)) {

	# Don't need to mess around with non-MSA files
	if ($filename =~ /\.afa$/) {

	    # Add to our list of directory contents
	    push(@DirContents,$filename);
	    
	    # Check if there's a disagreement file
	    $filename =~ s/\.afa$/\_ARFs/;
	    $filename =  $MultiMSADir{$Species[$i]}.'/'.$filename;

	    # Open it up and record the first line to our big disagreement file
	    if (-e $filename) {
		open(my $geneDis,'<',$filename);
		my $discountline = <$geneDis>;
		close($geneDis);
		$filename =~ /([^\/]+)\_ARFs$/;
		print $DisFile "$1: $discountline";
	    }

	}
    }
    closedir($dir);
    close($DisFile);

    
    # If we didn't have any disagreements, kill the disagreement file
    system("rm $DisFilename") if (-e $DisFilename && !(-s $DisFilename));

    # We're also now splitting all of the 'Hit' data off into gene-family-
    # specific files, so presumably we should also be free to kill this file
    system("rm \"$HitFileName\"");


    # Make a hash of all the names of sequences that we missed.
    # This will include any sequences that hit to a different 
    # chromosome from the group's "leader" (arbitrary -- first seen)
    my %QuilterMisses;
    if (-e $MissFileName) {
	open(my $MissFile,'<',$MissFileName);
	while (my $line = <$MissFile>) {
	    $line =~ s/\n|\r//g;
	    next if (!$line);
	    $line =~ /\|([^\|\s]+)\s*$/;
	    my $fam_name = lc($1);
	    if ($QuilterMisses{$fam_name}) {
		# Some sequences have commas in the name,
		# so we use ampersands (sp?) to break things up, instead
		$QuilterMisses{$fam_name} = $QuilterMisses{$fam_name}.'&'.$line;
	    } else { 
		$QuilterMisses{$fam_name} = $line;
	    }
	}
	close($MissFile);
    }
    
    # We'll collect per-family miss information after attempting to identify partial
    # mappings, so we can get rid of this file.
    system("rm \"$MissFileName\"");
    
    
    # You'll also need to know who had decent (but presumably not good enough)
    # Spaln mapping.  These will be recorded as both Misses and as their own thing
    #
    # NOTE: The chromosomes labeled in this file will have '[revcomp]' if revcomp
    #
    my %QuilterNearHits;
    if (-e $NearFileName) {

	open(my $NearFile,'<',$NearFileName);
	while (my $line = <$NearFile>) {

	    $line =~ s/\n|\r//g;
	    next if (!$line);

	    $line =~ /^(\S+) (\S+) (\S+)$/;
	    my $seqname = $1;
	    my $chromosome = $2;
	    my $maprange = $3;

	    $seqname =~ /\|([^\|]+)$/;
	    my $fam_name = lc($1);

	    # First off, this is technically a miss
	    if ($QuilterMisses{$fam_name}) { 
		$QuilterMisses{$fam_name} = $QuilterMisses{$fam_name}.'&'.$seqname;
	    } else {
		$QuilterMisses{$fam_name} = $seqname;
	    }

	    # Second off, this was a close hit -- keep the three fields single-space separated
	    # NOTE: Since these may end up being provided as commandline arguments,
	    #       we'll need them to be all quoted-up
	    $seqname    = '"'.$seqname.'"';
	    $chromosome = '"'.$chromosome.'"';
	    $maprange   = '"'.$maprange.'"';
	    my $triple  = $seqname.' '.$chromosome.' '.$maprange;
	    if ($QuilterNearHits{$fam_name}) {
		$QuilterNearHits{$fam_name} = $QuilterNearHits{$fam_name}.'&'.$triple;
	    } else {
		$QuilterNearHits{$fam_name} = $triple;
	    }

	}

	close($NearFile);

    }


    # Figure out which groups have MSAs in this species
    foreach my $filename (@DirContents) {

	if (-s $MultiMSADir{$Species[$i]}.'/'.$filename) {

	    # Count '>'s (number of sequences)
	    open(my $MultiMSAFile,'<',$MultiMSADir{$Species[$i]}.'/'.$filename);
	    my $countedIsos = 0;
	    while (my $line = <$MultiMSAFile>) {
		$countedIsos++ if ($line =~ /^\>/);
	    }
	    close($MultiMSAFile);

	    # Do we know this sequence?
	    if ($GroupsBySpecies{$filename}) {
		push(@{$GroupsBySpecies{$filename}},$i+1);
		push(@{$GroupsBySpecies{$filename}},$countedIsos);
	    } else {
		${$GroupsBySpecies{$filename}}[0] = $i+1;
		${$GroupsBySpecies{$filename}}[1] = $countedIsos;
	    }

	} elsif (-e $MultiMSADir{$Species[$i]}.'/'.$filename) {
	    # Wipe out any files that are only nominally present
	    my $clearMSACmd = 'rm '.$MultiMSADir{$Species[$i]}.'/'.$filename;
	    system($clearMSACmd);
	}

    }

    # Now that we have MSAs for all the groups that we had hits for,
    # align any that we missed using Needleman-Wunsch.
    print "  Aligning missed isoforms to gene MSAs\n" if ($verbose && keys %QuilterMisses);
    my $singleseqfilename = $SpeciesDir{$Species[$i]}.'Mirage.Single.Temp.fa';
    my $singleseqresults  = $SpeciesDir{$Species[$i]}.'Mirage.Multi.Temp.afa';

    foreach my $missedgroup (keys %QuilterMisses) {
	my @MissedSeqs = split('&',$QuilterMisses{$missedgroup});
	foreach my $seq (@MissedSeqs) {
	    my $GBSref = AttachSeqToMSA($seq,$MultiMSADir{$Species[$i]},$missedgroup,
					$SpeciesDBNames[$i],\%GroupsBySpecies,$i+1,
					$singleseqfilename,$singleseqresults);
	    %GroupsBySpecies = %{$GBSref};
	}
    }

    system("rm $singleseqfilename") if (-e $singleseqfilename);
    system("rm $singleseqresults")  if (-e $singleseqresults);

    # Knock it off with that darn timing!
    $MultiMSATimeStats[$i] = Time::HiRes::tv_interval($IntervalStart);

    # Species[i] over and out!
    ClearProgress();
    print "Intra-species alignment complete for $Species[$i]\n";

}


# Cool beans, daddy-o!  Now we're going to make sure that all of
# the species that didn't make it into our guide get a little bit
# of attention (and gather in their names).
my $GBSref;
my $MMSAref;
($speciesref,$numSpecies,$GBSref,$MMSAref)
    = CoverMinorSpecies(\@Species,$ProteinDB,$numSpecies,
			\%GroupsBySpecies,$AllSpeciesDir,\%MultiMSADir);
@Species         = @{$speciesref};
%GroupsBySpecies = %{$GBSref};
%MultiMSADir     = %{$MMSAref};


# Generate first progress message
my $progmsg = '  Mirage.pl: Preparing to generate final MSAs';
while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
print "$progmsg\r" if (!$verbose);

# We're now going to track how long the whole final-MSA-generating
# part of the program takes
$IntervalStart = [Time::HiRes::gettimeofday()];


# Make a directory to place our final MSAs in
my $FinalDir = $ResultsDir.'/FinalMSAs';
if (system("mkdir $FinalDir")) { die "\n  Failed to create directory '$FinalDir'\n\n"; }


# Convert the hash entries to an array, so we can divvy it up
# among the threads.  Did I mention there were going to be threads?
# Oh, yeah, you'd better believe there are going to be threads.
my @GroupFileNames   = sort keys %GroupsBySpecies;
my $TotNumGroupFiles = @GroupFileNames;
if ($TotNumGroupFiles < $num_cpus) { $num_cpus = $TotNumGroupFiles; }


# Carve out a portion of the files for each process to work on,
# or else die.  This would be really strange to see happen, but
# it's good to cover your bases.
my $portion;
if ($num_cpus) {
    $portion = $TotNumGroupFiles/$num_cpus;
} else {
    print "\n";
    print "  No genes were identified for alignment.\n";
    print "  Program Terminating Peacefully\n";
    print "\n\n";
    exit(0);
}


# Let the user know that we're about to hit the (long) home stretch
print "\n  Computing final cross-species MSAs (large genes can take up to 30 seconds)\n\n"
    if ($verbose);


# Split up the work across the desired number of processes
my $tempdirname = $ResultsDir.'/temp/';
if (system("mkdir \"$tempdirname\"")) { die "\n  ERROR:  Failed to create Mirage tempdir '$tempdirname'\n\n"; }
my $progressbase = $tempdirname.'mirage.thread_progress.';
my $processes = 1;
my $ThreadID  = 0; # In case of only one process
my $pid;
while ($processes < $num_cpus) {

    # Create a new thread
    if ( $pid = fork ) {
	# Not liking this pid? Fuhget about it!
	if (not defined $pid) { die "\n\tFork failed\n\n"; }	
	# Master thread maintains ThreadID of 0
	$ThreadID = 0;
    } else {
        $ThreadID = $processes;
	last;
    }

    $processes++;

}


# Each thread writes its progress (number of completed families)
# to a file containing 
my $progtime = time();
my $progfilename = $progressbase.$ThreadID;
my $num_complete = 0;


# Figure out where you'll be starting and stopping
my $startpoint = $ThreadID * $portion;
my $endpoint   = ($ThreadID+1) * $portion;
if ($ThreadID == $num_cpus-1) { $endpoint = $TotNumGroupFiles; }


#
# M2: At some point we're going to have to transfer back to the original names.
#     Good thing we've moved DCE info to come after a '#' in the seqnames...
#


# Now we go through each sequence group and align it across species
my $tempfileA = $tempdirname.'tempfileA.'.$ThreadID.'.MIRAGE.afa';
my $tempfileB = $tempdirname.'tempfileB.'.$ThreadID.'.MIRAGE.afa';
for (my $groupindex = $startpoint; $groupindex < $endpoint; $groupindex++) {

    # First, at your leisure (once every couple seconds) write
    # out how far you've gotten to a file, unless you're thread 0
    # who gets a very special job
    if (!$verbose) {
	if (time()-$progtime > 2) {
	    
	    if ($ThreadID) {
		open(my $progfile,'>',$progfilename);
		print $progfile "$num_complete\n";
		close($progfile);
		
	    } else {
		
		# Generate call to ProgressTimer
		my $progressCmd = $location.'ProgressTimer.pl '.$progressbase.' '.$num_cpus;
		$progressCmd    = $progressCmd.' '.$num_complete.' |';
		if ($location !~ /^\//) { $progressCmd = './'.$progressCmd; }
		
		# Make the call and read the results
		open(my $progfile,$progressCmd);
		my $overall_progress = readline($progfile);
		close($progfile);
		$overall_progress =~ s/\n|\r//g;
		
		# Pull together the message
		$progmsg = '  Mirage.pl:   '.$overall_progress.' groups aligned';
		while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
		$progmsg = $progmsg."\r";
		
		# Print it out!
		print "$progmsg";
		
	    }
	    
	    # The time for progress is NOW! People, rise up!
	    $progtime = time();
	    
	}
    }

    
    #
    # Now let's get to work! Time is money, people!
    #
    
    my $groupfile = $GroupFileNames[$groupindex];
    my $groupID = $groupfile;
    $groupID =~ s/\.afa$//;

    my @GroupSpecies;
    my @GroupSeqs;
    my $speciesPerGroup = 0;

    $i = 0;
    while ($i < @{$GroupsBySpecies{$groupfile}}) {
	$GroupSpecies[$i/2] = $Species[${$GroupsBySpecies{$groupfile}}[$i]-1];
	$GroupSeqs[$i/2]    = ${$GroupsBySpecies{$groupfile}}[$i+1];
	$speciesPerGroup++;
	$i += 2;
    }

    #print "  -- $gene: $speciesPerGene (@GeneSpecies)\n" if ($verbose);

    my $ResultFile = $FinalDir.'/'.$groupfile;

    if ($speciesPerGroup == 1) {

	if ($cleanMSA) {

	    # still want to remove any splice site characters
	    my $FinalMSACmd = 'perl '.$location.'FinalMSA.pl "'.$MultiMSADir{$GroupSpecies[0]}.'/'.$groupfile.'"';
	    $FinalMSACmd    = $FinalMSACmd." '".$ResultFile."'";
	    
	    # Print the command and execute
	    #print "  $FinalMSACmd\n" if ($verbose);
	    if (system($FinalMSACmd)) { 
		system("rm -rf $tempdirname");
		die "\n  *  ERROR: FinalMSA failed during execution  *\n\n"; 
	    }
	    
	} else {
	
	    # Just copy it over
	    if (system("cp '$MultiMSADir{$GroupSpecies[0]}/$groupfile' '$ResultFile'")) {		
		die "\n  Failed to move '$MultiMSADir{$GroupSpecies[0]}/$groupfile' to '$ResultFile'\n\n";
	    }

	}

	# Don't be shy, now!
	next;

    }

    my $fileA     = $MultiMSADir{$GroupSpecies[0]}.'/'.$groupfile;
    my $fileAsize = $GroupSeqs[0];
    my $fileB     = $MultiMSADir{$GroupSpecies[1]}.'/'.$groupfile;
    my $fileBsize = $GroupSeqs[1];

    my $MultiSeqNWCmd = $location."MultiSeqNW '$fileA' $fileAsize"; # First species
    $MultiSeqNWCmd    = $MultiSeqNWCmd." '$fileB' $fileBsize";      # Second species
    $MultiSeqNWCmd    = $MultiSeqNWCmd.' > '.$tempfileA;            # Send to temp file    
    if ($location !~ /^\//) { $MultiSeqNWCmd = './'.$MultiSeqNWCmd; }

    print "  $MultiSeqNWCmd\n" if ($verbose);
    if (system($MultiSeqNWCmd)) { die "\n  *  ERROR: MultiSeqNW failed during execution of command '$MultiSeqNWCmd'  *\n\n"; }
    
    $i = 2;
    while ($i < $speciesPerGroup) {

	if ($i % 2) { $fileA = $tempfileB; }
	else        { $fileA = $tempfileA; }

	$fileAsize += $fileBsize;
	$fileB      = $MultiMSADir{$GroupSpecies[$i]}.'/'.$groupfile;
	$fileBsize  = $GroupSeqs[$i];
	$MultiSeqNWCmd = $location."MultiSeqNW '$fileA' $fileAsize"; # First species
	$MultiSeqNWCmd = $MultiSeqNWCmd." '$fileB' $fileBsize";      # Second species
	if ($location !~ /^\//) { $MultiSeqNWCmd = './'.$MultiSeqNWCmd; }

	# The tempfile we're writing to depends on how many iterations 
	# we've been through.
	if ($i % 2) { $MultiSeqNWCmd = $MultiSeqNWCmd.' > '.$tempfileA; }
	else        { $MultiSeqNWCmd = $MultiSeqNWCmd.' > '.$tempfileB; }

	#print "  $MultiSeqNWCmd\n" if ($verbose);
	if (system($MultiSeqNWCmd)) { die "\n  *  ERROR: MultiSeqNW failed during execution of command '$MultiSeqNWCmd'  *\n\n"; }

	$i++;

    }

    # Perform the final cleanup! (unless we're CRAZY!)
    if ($cleanMSA) {

	my $FinalMSACmd = 'perl '.$location.'FinalMSA.pl ';
	if ($i % 2) { $FinalMSACmd = $FinalMSACmd."'$tempfileB'"; }
	else        { $FinalMSACmd = $FinalMSACmd."'$tempfileA'"; }
	$FinalMSACmd = $FinalMSACmd." '$ResultFile'";
	
	# Print the FinalMSA command and execute
	#print "  $FinalMSACmd\n" if ($verbose);
	if (system($FinalMSACmd)) { die "\n  *  ERROR: FinalMSA failed during execution  *\n\n"; }    

    } else { # Ew! That MSA is gonna be riddled with intron markers!

	my $FinalMoveCmd = 'mv ';
	if ($i % 2) { $FinalMoveCmd = $FinalMoveCmd."'$tempfileB'"; }
	else        { $FinalMoveCmd = $FinalMoveCmd."'$tempfileA'"; }
	$FinalMoveCmd = $FinalMoveCmd." '$ResultFile'";

	# Print the move command and execute
	#print "  $FinalMoveCmd\n" if ($verbose);
	if (system($FinalMoveCmd)) { die "\n  *  ERROR: FinalMoveCmd failed during execution  *\n\n"; }

    }

    print "  - $groupID complete\n" if ($verbose);
    $num_complete++;

}


# Get rid of the temporary files
if (-e $tempfileA && system("rm $tempfileA")) { die "\n  Failed to remove file $tempfileA\n\n"; }
if (-e $tempfileB && system("rm $tempfileB")) { die "\n  Failed to remove file $tempfileB\n\n"; }


# De-fork
if ($ThreadID) { exit(0); }
while (wait() != -1) {}


# Well, why NOT put this here?!
for ($i=0; $i<$numSpecies; $i++) {
    system("rm \"$SpeciesDBNames[$i]\"")      if (-e $SpeciesDBNames[$i]);
    system("rm \"$SpeciesDBNames[$i]\.ssi\"") if (-e $SpeciesDBNames[$i].'.ssi');
}
system("rm \"$MinorSpeciesDBName\"")      if (-e $MinorSpeciesDBName);
system("rm \"$MinorSpeciesDBName\.ssi\"") if (-e $MinorSpeciesDBName.'.ssi');


# Get rid of progress files (won't apply to thread 0)
for ($i = 1; $i < $num_cpus; $i++) {
    $progfilename = $progressbase.$i;
    if (-e $progfilename && system("rm $progfilename")) { die "\n  Failed to remove file $progfilename\n\n"; }
}
$progmsg = ' ';
while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
print "$progmsg\r" if (!$verbose);


# When did we finish generating our final MSAs?
$MultiSeqNWTime = Time::HiRes::tv_interval($IntervalStart);
$AvgNWTime = $MultiSeqNWTime / $TotNumGroupFiles;


# One final thing: Move any suspected transposable elements to
# a special final folder
if (grep -f, glob $tempdirname.'*.TEs.Quilter.out') {
    system("mkdir $ResultsDir/SuspectedTEs");
    system("mv ".$tempdirname."*.TEs.Quilter.out $ResultsDir/SuspectedTEs");
}


# Remove the temp directory
if (system("rm -rf $tempdirname")) { die "\n  ERROR:  Failed to delete temp directory '$tempdirname'\n\n"; }


# Slap that stop-watch!
$TotalRuntime = Time::HiRes::tv_interval($StartTime);

# Print out runtime statistitcs
if ($verbose || $timed) {
    my $formattedTime;
    print "\n\n\n";
    print "  +---------------------------------------------------+\n";
    print "                   Runtime Statistics \n";
    print "  +---------------------------------------------------+\n";
    foreach $i (0..$numSpecies-1) {
	
	print "\n";
	print "    Species        : $Species[$i]\n";
	
	$formattedTime = sprintf("%.3f",$QuilterTimeStats[$i]);
	print "    Hit Stitching  : $formattedTime seconds\n";
	
	$formattedTime = sprintf("%.3f",$MultiMSATimeStats[$i]);
	print "    MSA Generation : $formattedTime seconds\n";
	print "\n";
	
    }
    
    print "\n";
    $formattedTime = sprintf("%.3f",$MultiSeqNWTime);
    print "    Final MSA Generation  : $formattedTime seconds\n";
    
    $formattedTime = sprintf("%.3f",$AvgNWTime);
    print "    Avg MSA Time per Gene : $formattedTime seconds\n";
    
    print "\n\n";
    my $TotalRunMins = int($TotalRuntime / 60);
    my $TotalRunSecs = $TotalRuntime - (60 * $TotalRunMins);
    $formattedTime = sprintf("%d minutes %.3f seconds",$TotalRunMins,$TotalRunSecs);
    print "    TOTAL RUNTIME : $formattedTime \n";
    print "\n\n";
    
    print "  +---------------------------------------------------+\n";
    print "\n\n\n";

}

# Let the user know where the final directory is
print "\n\n  All tasks complete.\n  Results in $FinalDir\n\n\n";

# No more progress to be made!
system("rm -rf \"\$progress_dirname\"");


# What, that's all you got? Pssssh, shoulda known it was gonna be EZ PZ :p
1;








#######################
#######################
####               ####
####  SUBROUTINES  ####
####               ####
#######################
#######################







########################################################################
#
# Function Name: PrintUsage
#
# About: This function gives summary information about the program and
#        its use.  Printed to user if incorrect number of commandline
#        arguments are provided.
#
sub PrintUsage
{
    print "\n\n";
    print " MIRAGE : Multiple-sequence Isoform Alignment Tool Guided by Exon Boundaries ($mirage_version)\n\n";
    print " USAGE  : mirage [OPT.s] <Isoform DB>  <Species Guide>                                  \n\n";
    print " ARG.s  : <IsoformDB>     : A FASTA-formatted protein database, with the following      \n";
    print "                            naming convention:                                          \n\n";
    print "                            >gene_name|protein_name|species|seqID|groupID               \n\n";
    print "          <Species Guide> : A file indicating, for each species being searched on,      \n";
    print "                            the location of a FASTA-formatted genome for that species   \n";
    print "                            and a .gtf index file corresponding to that species and     \n";
    print "                            genome.  These fields should be whitespace-separated and    \n";
    print "                            ordered as follows:  species, genome, .gtf index            \n";
    print "                            It is recommended that similar species are grouped together \n";
    print "                            and positioned near the top of the list.                    \n\n";
    print " OPT.s  : --help      : More detailed help.                                             \n";
    print "          --verbose   : Verbose output.                                                 \n";
    print "          --time      : Print timing data to stdout at end of program                   \n";
    print "          -outdirname : Specify output directory name.                                  \n";
    print "          -cpus       : Specify number of CPU cores (default: 2)                        \n";
    die "\n\n";
}





########################################################################
#
# Function Name: Detailed Usage
#
# About: This function provides more details about the program than
#        the standard "input-assist"
#
sub DetailedUsage
{
    print "\n\n";
    print " .-======-+--=-----------------=----=--=--------------=---------=---------------. \n";
    print " | MIRAGE :  Multiple-sequence Isoform Alignment Tool Guided by Exon Boundaries | \n";
    print " '-======-+--=-----------------=----=--=--------------=---------=---------------' \n";
    print "  $mirage_version\n";
    print "                                                                                  \n";
    print "    USAGE :  mirage  [OPT.s]  <Isoform DB>  <Species Guide>                       \n";
    print "                                                                                  \n";
    print "                                                                                  \n";
    print "    ARG.s :  Isoform DB    : A FASTA-formatted protein database using the         \n";
    print "                             following sequence naming convention:                \n";
    print "                                                                                  \n";
    print "                             >gene_name|protein_name|species|seqID|groupID        \n";
    print "                                                                                  \n";
    print "             Species Guide : A simple file indicating, for each species being     \n";
    print "                             searched on, the location of a FASTA-formatted       \n";
    print "                             genome for that species and a .gtf index file        \n";
    print "                             corresponding to that species and genome.            \n";
    print "                                                                                  \n";
    print "                             The filenames should be paths to the files from      \n";
    print "                             the directory containing MIRAGE.pl and must be       \n";
    print "                             separated by whitespace.                             \n";
    print "                                                                                  \n";
    print "                             Required order:   species, genome, index file        \n";
    print "                                                                                  \n";
    print "                               EXAMPLE:                                           \n";
    print "                             .------------------------------------------------.   \n";
    print "                             | human  ~/data/HumanGenome.fa  ~/data/human.gtf |   \n";
    print "                             | mouse  ~/data/MouseGenome.fa  ~/data/mouse.gtf |   \n";
    print "                             | ...                                            |   \n";
    print "                             '------------------------------------------------'   \n";
    print "                                                                                  \n";
    print "                             Species are introduced into the final genewise MSAs  \n";
    print "                             in the order that they are listed in the species     \n";
    print "                             guide.  It is recommended that similar species are   \n";
    print "                             grouped together for optimal alignments.             \n";
    print "                                                                                  \n";
    print "                                                                                  \n";
    print "    OPT.s :  --verbose            : Verbose output                                \n";
    print "             --time               : Print timing data to stdout at end of program \n";
    print "             -outdirname <string> : Specify ouptut directory name                 \n";
    print "             -cpus <int>          : Specify number of CPU cores (default: 2)      \n";
    die "\n\n\n";
}





########################################################################
#
# Function Name: PrintVersion
#
# About:  Print out the version of Mirage being used
#
sub PrintVersion { die "\n  Mirage $mirage_version\n\n"; }





########################################################################
#
# Function Name: CheckInstall
#
# About:  Verify that spaln is running and that everything is where
#         we hope it would be.
#
sub CheckInstall
{

    # Figure out what the location of the Mirage src directory is
    my $location = $0;
    $location =~ s/Mirage2\.pl$//;

    my @RequiredFiles;
    push(@RequiredFiles,$location.'Quilter2.pl');
    push(@RequiredFiles,$location.'FastMap2.c');
    push(@RequiredFiles,$location.'MultiMSA.pl');
    push(@RequiredFiles,$location.'MultiSeqNW.c');
    push(@RequiredFiles,$location.'MultiSeqNW.h');
    push(@RequiredFiles,$location.'FinalMSA.pl');
    push(@RequiredFiles,$location.'makefile');
    push(@RequiredFiles,$location.'../inc/easel/miniapps/esl-sfetch');
    push(@RequiredFiles,$location.'../inc/easel/miniapps/esl-seqstat');
    push(@RequiredFiles,$location.'../inc/spaln2.3.3/src/spaln');
    push(@RequiredFiles,$location.'../inc/blat/blat.linux.x86_64');
    push(@RequiredFiles,$location.'../inc/blat/blat.macOSX.x86_64');
    push(@RequiredFiles,$location.'../inc/blat/blat.macOSX.i386');

    foreach my $file (@RequiredFiles) {
	if (!(-e $file)) { die "\n  Failure: Could not locate required file '$file'\n\n"; }
    }

    die "\n  Lookin' good!\n\n";

}





########################################################################
#
# Function Name: CheckSourceFiles
#
# About: This function verifies that all required files are present
#        in the src directory and, if anything needs to be compiled,
#        runs some compilation commands.
#
sub CheckSourceFiles
{
    my @RequiredFiles;
    push(@RequiredFiles,$location.'Quilter2.pl');
    push(@RequiredFiles,$location.'MultiMSA.pl');
    push(@RequiredFiles,$location.'FinalMSA.pl');
    push(@RequiredFiles,$location.'FastMap2');
    push(@RequiredFiles,$location.'HitWeaver');
    push(@RequiredFiles,$location.'MultiSeqNW');

    foreach my $srcfile (@RequiredFiles) {
	if (!(-e $srcfile)) {
	    die "\n  Failed to locate required file '$srcfile'\n\n";
	}
    }
}







########################################################################
#
# Function Name: ParseArgs
#
# About: This function takes the command line arguments provided by
#        the user and makes sense of them.
#
sub ParseArgs
{

    my %Options = ( 
	cpus => $defaultcpus, 
	outdirname => 'MirageResults',
	cleanmsa => 1,
	);

    &GetOptions( 
	\%Options,
	"help",
	"check",
	"version",
	"outdirname=s",
	"verbose",
	"cpus=i",
	"time",
	"stackarfs",
	"forcecompile", # Hidden
	"cleanmsa=i",   # Hidden
	"justspaln",    # Hidden
	"trackspaln",   # Hidden
	)
	|| die "\n  ERROR:  Failed to parse command line arguments\n\n";

    # If the user just wants a little help, I don't think it's too
    # difficult to give them a hand
    if ($Options{help})    { DetailedUsage(); }
    if ($Options{check})   { CheckInstall();  }
    if ($Options{version}) { PrintVersion();  }

    # If we don't have the required files, give usage
    my $proteindb = $ARGV[0];
    my $reference = $ARGV[1];
    if (!$proteindb) {
	if (!$reference) {
	    print "\n  ERROR : Protein database and reference file not provided\n";
	    print "  --------------------------------------------------------";
	    PrintUsage();
	} 
	print "\n  ERROR : Protein database not provided\n";
	print "  -------------------------------------";
	PrintUsage();
    }
    if (!$reference) {
	print "\n  ERROR : Reference file not provided\n";
	print "  -----------------------------------";
	PrintUsage();
    }

    # While we're here, we make sure there's a .ssi index for the
    # protein database
    ConfirmSSI($proteindb);

    return \%Options;
}








########################################################################
#
# Function Name: VerifiedClean
#
# About: Make sure that our database meets the standards for mirage.
#        Return 1 if clean, otherwise 0.
#
sub VerifiedClean
{
    # Get a hold of database name
    my $dbname = shift;

    # Open database
    my $DB = OpenInputFile($dbname);
    
    # If there are any comments starting it off, we leave them well enough alone
    my $line = <$DB>;
    while ($line && $line !~ /^>/) {
	$line = <$DB>;
    }

    # If this is the end of the file complain
    if (eof($DB)) {
	close($DB);
	die "\n  ERROR: Database '$dbname' does not appear to be FASTA-formatted\n\n";
    }

    # Check sequences
    while (!eof($DB)) {
	
	# Check the header
	$line =~ s/\n|\r//g;
	my $header = $line;
	if ($header =~ /\s|\#|\&|\//) {
	    close($DB);
	    print "\n  ERROR: Sequence name '$header' uses illegal characters\n";
	    return 0;
	}

	# Press on (verifying that this is an actual sequence)
	my $hascontent = 0;
	$line = <$DB>;
	while ($line && $line !~ /^\>/) {
	    $line =~ s/\n|\r//g;
	    $hascontent++ if ($line);
	    $line = <$DB>;
	}

	# If there isn't any substance to a sequence report it and quit
	if ($hascontent == 0) {
	    close($DB);
	    print "\n  ERROR: Sequence name '$header' does not belong to a sequence\n";
	    return 0;
	}

    }

    # Close the database
    close($DB);

    # Yes, we are happy with this database.
    return 1;
}







########################################################################
#
# Function Name: ParseSpeciesGuide
#
# About: This function reads the necessary "Species Guide" file, making
#        sure everything has the right (easel) indexing structures.
#
sub ParseSpeciesGuide
{
    my $guide_fname  = shift;
    my $protdb_fname = shift;

    my $GuideFile = OpenInputFile($guide_fname);

    # It's possible that we'll need to know what '~' means
    open(my $homecheck,"echo ~ |");
    my $home = <$homecheck>;
    close($homecheck);
    $home =~ s/\n|\r//g;
    $home =~ s/\/$//g;
    
    my $i = 0;
    my (@Species,@GTFs,@Genomes);

    while (my $line = <$GuideFile>) {

	$line =~ s/\n|\r//g;
	next if (!$line);
	
	if ($line =~ /^(\S+)\s+(\S+)/) {
		
	    $Species[$i] =  lc($1);
	    $Genomes[$i] =  $2;
	    $Genomes[$i] =~ s/^\~/$home/;
	    
	    if ($line =~ /\S+\s+\S+\s+(\S+)/) {
		$GTFs[$i] = $1;
		$GTFs[$i] =~ s/^\~/$home/;
	    } else {
		$GTFs[$i] = '-';
	    }
	    
	    # Verify that the genome exists
	    if (!(-e $Genomes[$i])) {
		print "\n  Failed to locate genome '$Genomes[$i]'\n\n"; 
		die   "  Recommendation: Verify path from directory containing MIRAGE.pl.\n\n";
	    }
	    
	    # Make sure the genome has a .ssi index
	    ConfirmSSI($Genomes[$i]);

	    # If we don't have a '.chrom.sizes' file (std name for UCSC) then
	    # generate one
	    my $chrsizes = $Genomes[$i];
	    $chrsizes =~ s/\.fa$/\.chrom.sizes/;
	    if (!(-e $chrsizes)) {
		my $ChrSizes = OpenOutputFile($chrsizes);
		my $SeqStats = OpenSystemCommand($eslseqstat." -a \"$Genomes[$i]\"");
		while (my $statline = <$SeqStats>) {
		    if ($statline =~ /^\= (\S+)\s+(\d+)/) {
			my $chr = $1;
			my $len = $2;
			print $ChrSizes "$chr\t$len\n";
		    }
		}
		close($SeqStats);
		close($ChrSizes);
	    }
	    
	    # Verify that the gtf exists (unless it's '-')
	    if ($GTFs[$i] ne '-' && !(-e $GTFs[$i])) {
		print "\n  Failed to locate GTF index '$GTFs[$i]'\n"; 
		die   "  Recommendation: Verify path from directory containing MIRAGE.pl.\n\n";
	    }
	    
	    $i++;

	} else {
	    
	    # Something isn't right
	    die "\n  Species guide line '$line' does not match expected format\n\n";
	    	    
	}

    }
    
    close($GuideFile);

    return (\@Species,\@GTFs,\@Genomes);
}








########################################################################
#
# Function Name: ConfirmSSI
#
# About:  This function takes a fasta-formatted database and either
#         generates an SSI file (if one doesn't already exist) or 
#         uses a bit of goofiness to test the accuracy of the existing
#         one.  If the existing one is inaccurate, it's replaced.
#
sub ConfirmSSI
{
    my $dbname = shift;

    # If there's no index, just make one already
    if (!(-e "$dbname\.ssi")) {
	if (system($eslsfetch." \-\-index $dbname 1>/dev/null")) {
	    die "\n  ERROR:  Failed to create easel index for \"$dbname\"\n\n";
	}
	return;
    }

    # Aw jeez, looks like we're going to be double-checking the current
    # index...
 
    # Slightly change the name of the database and make the link
    my $symlinkdb = $dbname;
    $symlinkdb =~ s/\.[^\.]+$/\.MirageSymlink\.fa/;
    my $symlinkcmd = "ln -s \"$dbname\" \"$symlinkdb\"";
    if (system($symlinkcmd)) {
	RunSystemCommand("rm \"$symlinkdb\"") if (-e $symlinkdb);
	die "\n  ERROR:  Symlink command '$symlinkcmd' failed\n\n";
    }

    # Build an ssi for the symlink
    if (system($eslsfetch." \-\-index $symlinkdb 1>/dev/null")) {
	RunSystemCommand("rm \"$symlinkdb\"");
	RunSystemCommand("rm \"$symlinkdb\.ssi\"") if (-e $symlinkdb.'.ssi');
	die"\n  ERROR:  Failed to create easel index for \"$symlinkdb\"\n\n";
    }

    # Check for differences between the two databases
    my $diffcmd = "diff \"$dbname\.ssi\" \"$symlinkdb\.ssi\" \|";
    open(my $diff,$diffcmd) || die "\n  ERROR:  Failed to diff the SSIs for '$dbname' and '$symlinkdb'\n\n";
    my $diffline = <$diff>;
    close($diff);

    # If we have a substantive difference, the existing SSI needs replacement!
    if ($diffline) { system("mv \"$symlinkdb\.ssi\" \"$dbname\.ssi\""); }
    else           { system("rm \"$symlinkdb\.ssi\"");                  }
    system("rm \"$symlinkdb\"");

    return;

}






########################################################################
#
# Function Name: GenerateSpeciesDBs
#
# About:  This function breaks up the isoform database into individual
#         species-specific databases, organized by gene families.
#
#         It also identifies the specific line ranges for each of the
#         requested processes, so that each process can easily operate
#         on the level of gene families, rather than individual sequences.
#
#         Because the naming conventions might be a bit of hot rubbish,
#         here's a super fun little guide to what this function returns:
#
#             + SpeciesDBNames : The filenames for the species-specific
#                 databases, organized in the same order as the Species
#                 array.
#
#             + MinorSpeciesDBName : The filename for the database of all
#                 isoforms belonging to species that don't have an
#                 associated genome.
#
#             + SpecThreadCounts : The number of threads to use for each
#                 species.  This will typically be num_cpus, but in
#                 the event that we have some species with fewer gene
#                 families than requested processes we'll reduce this.
#
#             + SpecFamBreaks : The line number that each thread should
#                 count up to (but not hit).
#
# M2: To be a little bit more flexible, we'll accept sequence names to
#     formatted either as:
#
# [1] >species|gene|iso_num 
#
#     or using the UniProt format:
#
# [2] >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier
#      [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
#
#     Where:
#     1. db is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
#     2. UniqueIdentifier is the primary accession number of the UniProtKB entry.
#     3. EntryName is the entry name of the UniProtKB entry.
#     4. ProteinName is the recommended name of the UniProtKB entry as annotated
#        in the RecName field. 
#        For UniProtKB/TrEMBL entries without a RecName field, the SubName field
#        is used. 
#        In case of multiple SubNames, the first one is used. 
#        The 'precursor' attribute is excluded, 'Fragment' is included with the
#        name if applicable.
#     5. OrganismName is the scientific name of the organism of the UniProtKB entry.
#     6. OrganismIdentifier is the unique identifier of the source organism,
#        assigned by the NCBI.
#     7. GeneName is the first gene name of the UniProtKB entry. 
#        If there is no gene name, OrderedLocusName or ORFname, the GN field is
#        not listed.
#     8. ProteinExistence is the numerical value describing the evidence for the
#        existence of the protein.
#     9. SequenceVersion is the version number of the sequence.
#
#     NOTE that a similar format (which we'll also accept) is the 'Alternative Isoform'
#     format that UniProt uses, where we have the following:
#
# [3] >sp|IsoID|EntryName Isoform IsoformName of ProteinName OS=OrganismName
#      OX=OrganismIdentifier[ GN=GeneName]
#
#     Where:
#     1. IsoID is the isoform identifier as assigned in the ALTERNATIVE PRODUCTS
#        section of the UniProtKB entry.
#     2. IsoformName is the isoform name as annotated in the ALTERNATIVE PRODUCTS
#        Name field of the UniProtKB entry.
#     3. ProteinExistence and SequenceVersion do not apply to alternative isoforms
#
#     NOTE: Cases like [2] and [3] will be converted into the [1] format for search,
#           so we'll also assign each sequence a unique number to index into a list of
#           original names.
#
#           Also, [3] is just a special variant of [2], so we really only need to have
#           two acceptable formats for parsing.
#
sub GenerateSpeciesDBs
{
    my $ProteinDB_name = shift;
    my $num_cpus       = shift;
    my $specseqdir_ref = shift;

    my %SpeciesSeqDir = %{$specseqdir_ref};

    # The first thing we'll do is scan through our database
    # making a set of species-AND-gene-family mini-databases.

    # We'll have a hash that counts the number of sequence characters,
    # both for species as wholes and specific genes within species
    my %CharCounts;

    #ConfirmSSI($ProteinDB_name); # We do this in 'ParseArgs'
    my $ProteinDB = OpenInputFile($ProteinDB_name);
    my $seq_num = 0;
    my @OriginalSeqNames;
    my $spec_gene_filename;
    my $SpecGeneFile;
    my $species;
    my $genefam;
    while (my $line = <$ProteinDB>) {
	
	$line =~ s/\n|\r//g;
	next if (!$line);
	next if ($line =~ /^\s+\#/); # No comment lines in the search dbs, okay folks?

	if ($line =~ /^\>/) {

	    # The full name may include comment content, which we'll assume starts with
	    # a '#' character
	    my $seqname = $line;
	    $seqname =~ s/^\>//;
	    $seqname =~ /^([^\#]+)/;
	    my $parsename = $1;
	    
	    # Figure out what format this sequence name is using and play in the space
	    # with it.
	    if (ParseSeqNameAsMirage($parsename)) {

		$parsename =~ /^([^\|]+)\|([^\|]+)\|/;
		$species = lc($1);
		$genefam = lc($2);

	    } elsif (ParseSeqNameAsUniProt($parsename)) {

		$parsename =~ /OS\=(\S+)/;
		$species = lc($1);

		$parsename =~ /GN\=(\S+)/;
		$genefam = lc($1);

	    } else {
		die "\n  SHIT!\n\n";
	    }

	    # Record the original name of this sequence, and then switch over to our
	    # nice, minimal name format for the duration.
	    push(@OriginalSeqNames,$seqname);
	    $seqname = $seq_num;

	    # If we already have a file open, close it
	    if ($spec_gene_filename) { close($SpecGeneFile); }

	    # Determine the specific file we want to write to
	    $species = 'Misc' if (!$SpeciesDir{$species});
	    $spec_gene_filename = $SpeciesSeqDir{$species}.$genefam.'.fa';

	    # Print that beautiful nameline!
	    open($SpecGeneFile,'>>',$spec_gene_filename) || die "\n  ERROR:  Failed to open output species-gene database '$spec_gene_filename' ($species)\n\n";
	    print $SpecGeneFile ">$seq_num\n";

	    $seq_num++;

	} elsif ($spec_gene_filename) {

	    print $SpecGeneFile "$line\n";

	    # Wait! Before you move on, how many characters are we lookin' at here?
	    my $linelen = length($line);
	    if ($CharCounts{$species}) { $CharCounts{$species} += $linelen; }
	    else                       { $CharCounts{$species}  = $linelen; }

	    if ($CharCounts{$species.'&'.$genefam}) {
		$CharCounts{$species.'&'.$genefam} = $CharCounts{$species.'&'.$genefam} += $linelen;
	    } else {
		$CharCounts{$species.'&'.$genefam} = $CharCounts{$species.'&'.$genefam}  = $linelen;
	    }
	    
	}

    }
    if ($spec_gene_filename) { close($SpecGeneFile); }
    close($ProteinDB);

    # Next, we'll run through each of our species and determine which gene
    # families belong to which threads.
    foreach $species (keys %SpeciesSeqDir) {

	# You gotta catch this ol' trickster before it fouls things up!
	next if ($species eq 'Misc');

	# How many sequence characters belong to this family?
	# How many characters (on average) should go to each thread?
	my $tot_chars = $CharCounts{$species};
	my $thread_portion = int($tot_chars/$num_cpus);

	# Great! Now we can run through and assign each thread its share of sequences
	my $ThreadGuide = OpenOutputFile($SpeciesSeqDir{$species}.'Thread-Guide');
	my $threadnum   = 0;
	my $threadchars = 0;

	# We'll also print the number of threads at the top of the file,
	# to spare the option parsing.
	print $ThreadGuide "Num CPUs: $num_cpus\n";

	# This feels a tad inefficient, but in the grand scheme of things it should be
	# negligible -- esp. compared to the memory cost of hash-of-hashing
	foreach my $species_gene (keys %CharCounts) {

	    next if ($species_gene !~ /\&/);
	    $species_gene =~ /^(\S+)\&(\S+)$/;
	    my $sp   = $1; # 'species' is already in use :(
	    $genefam = $2;

	    next if ($sp ne $species);
	    
	    print $ThreadGuide "$threadnum $genefam\n";

	    # If we've hit our portion, jump ship!
	    $threadchars += $CharCounts{$species_gene};
	    if ($threadchars >= $thread_portion && $threadnum < $num_cpus-1) {
		$threadnum++;
		$threadchars = 0;
	    }

	}
	close($ThreadGuide);

    }
    
    # This is it, baby!
    return (\@OriginalSeqNames);

}





########################################################################
#
# Function Name: CoverMinorSpecies
#
# About: In order to handle the presence of species that aren't listed
#        in the species guide file, we do a little bit of extra work
#        making sure they get worked into the fold
#
sub CoverMinorSpecies
{
    # Grab the original (major) species
    my $speciesRef   = shift;
    my @Species      = @{$speciesRef};

    # Convert to a hash, for fun and games (and work and seriousness)
    my %MajorSpecies;
    foreach my $i (0..@Species-1) {
	$MajorSpecies{$Species[$i]} = 1;
    }
    
    # Grab the name of the database we're interested in looking through
    my $DBname = shift;

    # How many species were there to start with?
    my $num_species = shift;
    my %SpeciesIDs;

    # Grab a reference to the group hash.  Keep in mind that the
    # values stored in this hash are (species_index,number_of_isoforms)
    # pairs
    my $grouphashRef    = shift;
    my %GroupsBySpecies = %{$grouphashRef};

    # The name of the results directory
    my $speciesDir = shift;
    $speciesDir    =~ s/\/$//;

    # The hash of intra-species MSA file locations within species
    my $MultiMSAref = shift;
    my %MultiMSADir = %{$MultiMSAref};

    # A temporary "GroupsBySpecies"-like hash
    my %MinorGroupsBySpecies;
    
    # Scan along the protein database and record any sequences that don't
    # belong to one of the major species.  It would be nice if this could
    # be parallelized, but we'd run the risk of processes trampling over
    # directories and files, so for the time being we'll just do a one-thread
    # scan.
    my $species;
    my %MinorSpecies;
    my $tempresults  = 'CoverMinorSpecies.results.fa';
    my $tempfilename = 'CoverMinorSpecies.temp.fa';
    open(my $DB,'<',$ARGV[0]);
    my $line = readline($DB); # Priming the regex
    while (!eof($DB)) {
	
	# Go to the start of the next sequence
	while (!eof($DB) && $line !~ /^\>/) { # "GN:" EXCISION
	    $line = readline($DB);
	    $line =~ s/\r|\n//g;
	}
	
	# Are we done?
	last if (eof($DB));
	
	# Not yet we ain't! Grab the species name
	my $header = $line;
	$line =~ /^\>([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|(\S*\|)?([^\|\s]+)\s*$/; # "GN:" EXCISION

	if (!$6) { print ">>> $line\n"; }
	
	my $gene    = uc($1);
	my $species = lc($3);
	my $groupID = lc($6);

	if (!$MajorSpecies{$species}) {
	    
	    # Hot diggity! We got ourselves a minor species!
	    $MinorSpecies{$species} = 1;

	    # Oo la la! A brand new one, at that!
	    if (!$SpeciesIDs{$species}) {
		$MultiMSADir{$species} = $speciesDir.'/'.$species;
		$SpeciesIDs{$species}  = $num_species;
		$Species[$num_species] = $species;
		$num_species++;
		if (system("mkdir $speciesDir/$species")) {
		    die "\n  ERROR: Failed to make directory '$speciesDir/$species'\n\n";
		}
	    }
	    
	    # Recording the sequence.
	    open(my $tempfile,'>',$tempfilename);
	    print $tempfile "$line\n";
	    $line = readline($DB);
	    while ($line !~ /\>/) {
		print $tempfile "$line";
		last if (eof($DB));
		$line = readline($DB);
	    }
	    print $tempfile "\n";
	    close($tempfile);
	    
	    # Have we hit on any other groups for this species?
	    if ($MinorGroupsBySpecies{$groupID.'#'.$species}) { 

		# Need to get all Needleman-Wunschy about it
		my $NWcmd = $location.'MultiSeqNW '.$tempfilename.' 1 '."'$speciesDir/$species/$groupID.afa' ";
		$NWcmd    = $NWcmd.$MinorGroupsBySpecies{$groupID.'#'.$species};
		$NWcmd    = $NWcmd.' > '.$tempresults;
		if ($location !~ /^\//) { $NWcmd = './'.$NWcmd; }
		if (system($NWcmd)) { die "\n  Failed to align non-indexed sequence '$header'\n\n"; }

		# Transfer the results
		my $mvCmd = "mv $tempresults '$speciesDir/$species/$groupID.afa'";
		if (system($mvCmd)) { die "\n  Failed to move '$tempfilename' to '$speciesDir/$species/$groupID.afa'\n\n"; }
		$MinorGroupsBySpecies{$groupID.'#'.$species}++; 

	    } else { 

		# Free pass!
		my $mvCmd = "mv $tempfilename '$speciesDir/$species/$groupID.afa'";
		if (system($mvCmd)) { die "\n  Failed to move '$tempfilename' to '$speciesDir/$species/$groupID.afa'\n\n"; }
		$MinorGroupsBySpecies{$groupID.'#'.$species} = 1;

	    }

	} else {
	    
	    # Prime the next search
	    $line = readline($DB) if (!eof($DB));
	    
	}
	
    }
    close($DB);
    
    
    # If we ever got around to writing into that tempfile, destroy it
    system("rm $tempfilename") if (-e $tempfilename);

    
    # NOW that we've logged all these fun new species, let's make sure they
    # make it onto our big datastructures
    foreach my $pair (keys %MinorGroupsBySpecies) {
	my @asArray  = split('#',$pair);
	my $filename = $asArray[0].'.afa';
	if ($GroupsBySpecies{$filename}) {
	    push(@{$GroupsBySpecies{$filename}},$SpeciesIDs{$asArray[1]}+1);
	    push(@{$GroupsBySpecies{$filename}},$MinorGroupsBySpecies{$pair});
	} else {
	    ${$GroupsBySpecies{$filename}}[0] = $SpeciesIDs{$asArray[1]}+1;
	    ${$GroupsBySpecies{$filename}}[1] = $MinorGroupsBySpecies{$pair};
	}
    }


    # Return the minor species
    return (\@Species,$numSpecies,\%GroupsBySpecies,\%MultiMSADir);
    
    
}





########################################################################
#
# Function Name: ParseSeqNameAsMirage
#
# About: species|gene|iso_num
#
sub ParseSeqNameAsMirage
{
    my $seqname = shift;
    return 1 if ($seqname =~ /^[^\|]+\|[^\|]+\|[^\|]+$/);
    return 0;
}


########################################################################
#
# Function Name: ParseSeqNameAsUniProt
#
# About: db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
#
#     Where:
#     1. db is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
#     2. UniqueIdentifier is the primary accession number of the UniProtKB entry.
#     3. EntryName is the entry name of the UniProtKB entry.
#     4. ProteinName is the recommended name of the UniProtKB entry as annotated
#        in the RecName field. 
#        For UniProtKB/TrEMBL entries without a RecName field, the SubName field
#        is used. 
#        In case of multiple SubNames, the first one is used. 
#        The 'precursor' attribute is excluded, 'Fragment' is included with the
#        name if applicable.
#     5. OrganismName is the scientific name of the organism of the UniProtKB entry.
#     6. OrganismIdentifier is the unique identifier of the source organism,
#        assigned by the NCBI.
#     7. GeneName is the first gene name of the UniProtKB entry. 
#        If there is no gene name, OrderedLocusName or ORFname, the GN field is
#        not listed.
#     8. ProteinExistence is the numerical value describing the evidence for the
#        existence of the protein.
#     9. SequenceVersion is the version number of the sequence.
#
#
sub ParseSeqNameAsUniProt
{
    my $seqname = shift;
    if ($seqname =~ /^([^\|]+)\|[^\|]+\|[^\|]+ /) {
	my $origin_db = lc($1);
	if ($origin_db ne 'sp' && $origin_db ne 'tr') { return 0; }
	# Make sure we can get everything we need to do our Mirage-y work with this
	# sequence.
	if ($seqname =~ /GN\=\S+/ && $seqname =~ /OS\=\S+/) { return 1; }
	return 0;
    }
    return 0;
}






########################################################################
#
# Function Name: AttachSeqToMSA
#
# About: This function attaches a single sequence to an MSA, and is
#        used to handle cases where neither FastMap nor SPALN
#        are able to identify a correct alignment of the sequence to
#        the genome.
#
sub AttachSeqToMSA
{

    my ($i,$j,$k);

    my $seqname   = shift;
    my $MSADir    = shift;
    my $groupname = shift;
    my $ProteinDB = shift;

    # As a reminder, GroupsBySpecies is a hash of arrays of
    # pairs where the first item is a numeric value (indicating
    # index into Species array) and the second item is the number
    # of isoforms for the given group name (the hash key).
    my $GBSref          = shift;
    my %GroupsBySpecies = %{$GBSref};
    my $species_index   = shift;

    # First, let's grab a hold of the missed sequence
    my $seqfilename  = shift;
    my $resultsname  = shift;
    my $origseq = $seqname;
    $origseq =~ s/^\>|\s//g;
    my $eslsfetchCmd = $eslsfetch." -o \"$seqfilename\" \"$ProteinDB\" \"$origseq\" > /dev/null 2>\&1";
    
    
    # run it!
    if (system($eslsfetchCmd)) { die "\n  ERROR: esl-sfetch command '$eslsfetchCmd' failed\n\n"; }


    # Now we're ready to join this sequence to the MSA (if there is one!!!)
    my $baseMSAfile = $MSADir.'/'.$groupname.'.afa';
    if (-s $baseMSAfile) {

	# How many seq.s are in the existing MSA?
	$i = 0;
	while (${$GroupsBySpecies{$groupname.'.afa'}}[$i] 
	       && ${$GroupsBySpecies{$groupname.'.afa'}}[$i] != $species_index) {
	    $i += 2;
	}

	# If we haven't shot over the edge (i.e., exhausted known species with this group), 
	# then we're appending
	if (${$GroupsBySpecies{$groupname.'.afa'}}[$i]) {
	    
	    my $num_isos = ${$GroupsBySpecies{$groupname.'.afa'}}[$i+1];
	    
	    if ($num_isos == 0) { # Final special check, I promise

		# We just call this our MSA, and make a note of it in our hash
		system("mv '$seqfilename' '$baseMSAfile'");
		${$GroupsBySpecies{$groupname.'.afa'}}[$i+1] = 1;

	    } else {

		# Compose and run that daaaaang command!
		my $alignCmd = $location."MultiSeqNW '".$seqfilename."' 1 '".$baseMSAfile."' ".$num_isos;
		$alignCmd    = $alignCmd.' -igBase 0';
		$alignCmd    = $alignCmd." > '".$resultsname."'";
		if ($location !~ /^\//) { $alignCmd = './'.$alignCmd; }
		if (system($alignCmd)) { die "  *  ERROR: Alignment of missing sequence '$origseq' to MSA failed ($alignCmd)\n\n"; }
		
		# Move the results and count this sequence
		system("mv $resultsname '$baseMSAfile'");
		${$GroupsBySpecies{$groupname.'.afa'}}[$i+1] = $num_isos+1;

	    }		

	} else {

	    # We just call this our MSA, and make a note of it in our hash
	    system("mv '$seqfilename' '$baseMSAfile'");
	    push(@{$GroupsBySpecies{$groupname.'.afa'}},$species_index);
	    push(@{$GroupsBySpecies{$groupname.'.afa'}},1);

	}

    } else {
	
	# We just call this our MSA, and make a note of it in our hash
	system("mv '$seqfilename' '$baseMSAfile'");
	push(@{$GroupsBySpecies{$groupname.'.afa'}},$species_index);
	push(@{$GroupsBySpecies{$groupname.'.afa'}},1);

    }
    
    # Consider your hash updated!
    return \%GroupsBySpecies;
    
}





# EOF

#!/usr/bin/env perl
#
# MIRAGE.pl - Multiple-sequence Isoform Alignment Tool Guided by Exon Boundaries
#           - Alex Nord
#           - 2016
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


sub PrintUsage;
sub DetailedUsage;
sub CheckInstall;
sub CheckSourceFiles;
sub ParseArgs;
sub VerifiedClean;
sub ParseSpeciesGuide;
sub GenerateSpeciesDBs;
sub CoverMinorSpecies;
sub AttachSeqToMSA;



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


# Figure out what the location of the Mirage src directory is
my $location = $0;
$location =~ s/Mirage\.pl$//;


# Generics
my ($i,$j,$k);


# We're going to need these friends
my $eslsfetch  = $location.'../inc/easel/miniapps/esl-sfetch';


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
my $ProteinDB    = $ARGV[0];
my $SpeciesGuide = $ARGV[1];
my $ResultsDir   = $Options{outdirname};
my $verbose      = $Options{verbose};
my $numProcesses = $Options{cpus};
my $timed        = $Options{time};
my $stack_arfs   = $Options{stackarfs};
my $forcecompile = $Options{forcecompile}; # Hidden
my $cleanMSA     = $Options{cleanmsa};     # Hidden
my $just_spaln   = $Options{justspaln};    # Hidden


# Verify that we have all the files we need on-hand
#CheckSourceFiles($forcecompile);


# Make sure that the two key files are for real
if (!(-e $ProteinDB))    { die "\n  Failed to locate file '$ARGV[0]'\n\n"; } 
if (!(-e $SpeciesGuide)) { die "\n  Failed to locate file '$ARGV[1]'\n\n"; } 


# Check if the protein database satisfies the requirements
# we've set (no names with #, &, or /, no whitespace in names)
if (!VerifiedClean($ProteinDB)) {
    print "\n";
    print "  Protein database '$ProteinDB' needs to be cleaned.\n";
    print "  Try command './CleanMirageDB.pl $ProteinDB'\n";
    die "\n";
}


# Convert species guide into three arrays.
# We'll ignore any species that don't show up
# in our protein database.
my ($speciesref,$gtfsref,$genomesref) = ParseSpeciesGuide($SpeciesGuide,$ProteinDB);
my @Species    = @{$speciesref};
my $numSpecies = @Species;
my @GTFs       = @{$gtfsref};
my @Genomes    = @{$genomesref};


# Make a directory to store results
if (-e $ResultsDir) {
    my $testDir = $ResultsDir;
    $i = 0;
    while (-e $testDir) {
	$i++;
	$testDir = $ResultsDir."_$i";
    }
    $ResultsDir = $testDir;
}
if (system("mkdir $ResultsDir")) { die "\n  Failed to create directory '$ResultsDir'\n\n"; }
print "\n  Results Directory:  $ResultsDir\n\n";


# Make species-specific results folders
my $specname;
my %SpeciesDir;
my $AllSpeciesDir = $ResultsDir.'/SpeciesMSAs';
if (system("mkdir $AllSpeciesDir")) { die "\n  Failed to create directory '$AllSpeciesDir'\n\n"; }
foreach $specname (@Species) {
    $SpeciesDir{$specname} = $AllSpeciesDir.'/'.$specname;
    if (system("mkdir $SpeciesDir{$specname}")) {
	die "\n  Failed to create directory '$SpeciesDir{$specname}'\n\n";
    }
}


# Create a temp directory in src where we'll hide all of our secrets.
my $tempdirname = $location.'temp/';
if (-d($tempdirname)) { system("rm -rf $tempdirname"); }
if (system("mkdir $tempdirname")) { die "\n  ERROR:  Failed to generate temporary directory '$tempdirname'\n\n"; }


# Make species-specific databases, organized by gene family
my ($speciesdbnames_ref,$MinorSpeciesDBName,$specprocesses_ref,$fambreaks_ref) 
    = GenerateSpeciesDBs(\@Species,$ProteinDB,$numProcesses);
my @SpeciesDBNames   = @{$speciesdbnames_ref};
my @SpeciesProcesses = @{$specprocesses_ref};
my @SpeciesFamBreaks = @{$fambreaks_ref};


# Run Quilter
my %MultiMSADir;
my %GroupsBySpecies;
for ($i=0; $i<$numSpecies; $i++) {

    # Start a timer for Quilter
    $IntervalStart = [Time::HiRes::gettimeofday()];

    # Assembling the Quilter command.  Because we're a directory back
    # from src, we need to append a '../' to each of the file names.
    my $QuilterCmd;
    $QuilterCmd = 'perl '.$location.'Quilter.pl '; # base call

    # Depending on whether we're using a relative path (w.r.t. the current
    # directory) or absolute path, what we pass to Quilter will vary for
    # each of the files

    # Protein database
    $QuilterCmd = $QuilterCmd.$SpeciesDBNames[$i].' ';

    # Genome
    $QuilterCmd = $QuilterCmd.$Genomes[$i].' ';

    # GTF index (special case: - for "Just use BLAT")
    $QuilterCmd = $QuilterCmd.$GTFs[$i].' ';
	
    # Species name
    $QuilterCmd = $QuilterCmd.$Species[$i];
    
    # Guide output towards results directory
    $QuilterCmd = $QuilterCmd.' -o '.$SpeciesDir{$Species[$i]};

    # How many CPUs do we want to work with?
    my $speciesProcesses = $SpeciesProcesses[$i];

    # Tell Quilter how many CPUs we want it to work with (-nplus is -n plus stop points)
    $QuilterCmd = $QuilterCmd.' -nplus '.$speciesProcesses;

    # Tell Quilter which lines each process stops at
    for ($j=0; $j<$speciesProcesses; $j++) { $QuilterCmd = $QuilterCmd.' '.$SpeciesFamBreaks[$i][$j]; }

    # Do we want a whole bunch of stuff spat at us?
    $QuilterCmd = $QuilterCmd.' -v' if ($verbose);

    # Are we timing?
    $QuilterCmd = $QuilterCmd.' -time' if ($timed);

    # Display the call if we're being loud
    if ($verbose) {
	print "\n  Aligning to the $Species[$i] genome using Quilter.pl\n";
	print "  $QuilterCmd\n";
    }

    # Make that daaaaaang call.
    if (system($QuilterCmd)) {
	system("rm -rf $tempdirname");
	die "\n  *  ERROR:  Quilter.pl failed during execution  *\n\n";
    }

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
    #if (!(-e $HitFile)) { die "\n  Failed to locate Quilter output '$HitFile'\n\n"; }

    # Name the directory where we'll be putting the output from MultiMSA.pl
    $MultiMSADir{$Species[$i]} = $SpeciesDir{$Species[$i]}.'/MultiMSA_output';

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
    $MultiMSACmd = $MultiMSACmd.' -n '.$numProcesses;

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
	# Hamsters!
	system("rm -rf $tempdirname");
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
    my $singleseqfilename = $tempdirname.'Mirage.Single.Temp.fa';
    my $singleseqresults  = $tempdirname.'Mirage.Multi.Temp.afa';

    foreach my $missedgroup (keys %QuilterMisses) {
	my @seqnames = split('&',$QuilterMisses{$missedgroup});
	foreach my $seq (@seqnames) {
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
if (system("mkdir $FinalDir")) { 
    system("rm -rf $tempdirname");
    die "\n  Failed to create directory '$FinalDir'\n\n"; 
}


# Convert the hash entries to an array, so we can divvy it up
# among the threads.  Did I mention there were going to be threads?
# Oh, yeah, you'd better believe there are going to be threads.
my @GroupFileNames   = sort keys %GroupsBySpecies;
my $TotNumGroupFiles = @GroupFileNames;
if ($TotNumGroupFiles < $numProcesses) { $numProcesses = $TotNumGroupFiles; }


# Carve out a portion of the files for each process to work on,
# or else die.  This would be really strange to see happen, but
# it's good to cover your bases.
my $portion;
if ($numProcesses) {
    $portion = $TotNumGroupFiles/$numProcesses;
} else {
    system("rm -rf $tempdirname");
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
my $progressbase = $tempdirname.'mirage.thread_progress.';
my $processes = 1;
my $ThreadID  = 0; # In case of only one process
my $pid;
while ($processes < $numProcesses) {

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
if ($ThreadID == $numProcesses-1) { $endpoint = $TotNumGroupFiles; }


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
		my $progressCmd = $location.'ProgressTimer.pl '.$progressbase.' '.$numProcesses;
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
system("rm \"$MinorSpeciesDBName\"")      if (-e $MinorSpeciesDBName);
system("rm \"$MinorSpeciesDBName\.ssi\"") if (-e $MinorSpeciesDBName.'.ssi');


# Get rid of progress files (won't apply to thread 0)
for ($i = 1; $i < $numProcesses; $i++) {
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
    print " MIRAGE : Multiple-sequence Isoform Alignment Tool Guided by Exon Boundaries            \n\n";
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
# Function Name: CheckInstall
#
# About:  Verify that spaln is running and that everything is where
#         we hope it would be.
#
sub CheckInstall
{

    # Figure out what the location of the Mirage src directory is
    my $location = $0;
    $location =~ s/Mirage\.pl$//;

    my @RequiredFiles;
    push(@RequiredFiles,$location.'Quilter.pl');
    push(@RequiredFiles,$location.'DiagonalSets.pm');
    push(@RequiredFiles,$location.'FastDiagonals.c');
    push(@RequiredFiles,$location.'Diagonals.c');
    push(@RequiredFiles,$location.'Diagonals.h');
    push(@RequiredFiles,$location.'TransSW.c');
    push(@RequiredFiles,$location.'MultiMSA.pl');
    push(@RequiredFiles,$location.'MultiSeqNW.c');
    push(@RequiredFiles,$location.'MultiSeqNW.h');
    push(@RequiredFiles,$location.'FinalMSA.pl');
    push(@RequiredFiles,$location.'makefile');
    push(@RequiredFiles,$location.'../inc/easel/miniapps/esl-sfetch');
    push(@RequiredFiles,$location.'../inc/easel/miniapps/esl-seqstat');
    push(@RequiredFiles,$location.'../inc/spaln2.2.2/src/spaln');
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
    push(@RequiredFiles,$location.'Quilter.pl');
    push(@RequiredFiles,$location.'DiagonalSets.pm');
    push(@RequiredFiles,$location.'MultiMSA.pl');
    push(@RequiredFiles,$location.'FinalMSA.pl');
    push(@RequiredFiles,$location.'FastDiagonals');
    push(@RequiredFiles,$location.'TransSW');
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
	);

    &GetOptions( 
	\%Options,
	"help",
	"check",
	"outdirname=s",
	"verbose",
	"cpus=i",
	"time",
	"stackarfs",
	"forcecompile", # Hidden
	"cleanmsa",     # Hidden
	"justspaln",    # Hidden
	)
	|| die "\n  ERROR:  Failed to parse command line arguments\n\n";

    # If the user just wants a little help, I don't think it's too
    # difficult to give them a hand
    if ($Options{help})  { DetailedUsage(); }
    if ($Options{check}) { CheckInstall();  }

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
    if (!(-e $proteindb.'.ssi')) {
	if (system($eslsfetch." \-\-index $proteindb")) { 
	    die "\n  Failed to create easel index for $proteindb\n\n"; 
	}
    }

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
    open(my $db,'<',$dbname) || die "\n  Failed to open protein database '$dbname'\n\n";
    
    # If there are any comments starting it off, we leave them well enough alone
    my $line = <$db>;
    while ($line && $line !~ /^>/) {
	$line = <$db>;
    }

    # If this is the end of the file complain
    if (eof($db)) {
	close($db);
	die "\n  ERROR: Database '$dbname' does not appear to be FASTA-formatted\n\n";
    }

    # Check sequences
    while (!eof($db)) {
	
	# Check the header
	$line =~ s/\n|\r//g;
	my $header = $line;
	if ($header =~ /\s|\#|\&|\//) {
	    close($db);
	    print "\n  ERROR: Sequence name '$header' uses illegal characters\n";
	    return 0;
	}

	# Press on (verifying that this is an actual sequence)
	my $hascontent = 0;
	$line = <$db>;
	while ($line && $line !~ /^\>/) {
	    $line =~ s/\n|\r//g;
	    $hascontent++ if ($line);
	    $line = <$db>;
	}

	# If there isn't any substance to a sequence report it and quit
	if ($hascontent == 0) {
	    close($db);
	    print "\n  ERROR: Sequence name '$header' does not belong to a sequence\n";
	    return 0;
	}

    }

    # Close the database
    close($db);

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
    my $GuideFilename = shift;
    my $ProtFilename  = shift;

    open(my $GuideFile,'<',$GuideFilename) || die "\n  Failed to open species guide file '$GuideFilename'\n\n";

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

	    # Guarantee that this species is actually represented in
	    # the protein database.
	    my $spec = lc($1);
	    my $SpecGrepCmd = "grep -m 1 -i '|$spec|' $ProtFilename |";
	    open(my $SG,$SpecGrepCmd) || die "\n  ERROR:  Species grep command '$SpecGrepCmd' failed\n\n";
	    my $specgrepline = <$SG>;
	    close($SG);

	    # Is there anything to be stoked on?
	    if ($specgrepline) { 

		# Recover the old line regex
		$line =~ /^(\S+)\s+(\S+)/;
		
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
		if (!(-e $Genomes[$i].'.ssi')) {
		    if (system($eslsfetch." --index $Genomes[$i] 1>/dev/null")) { 
			die "\n  Failed to create easel index for $Genomes[$i]\n\n"; 
		    }
		}
		
		# Verify that the gtf exists (unless it's '-')
		if ($GTFs[$i] ne '-' && !(-e $GTFs[$i])) {
		    print "\n  Failed to locate GTF index '$GTFs[$i]'\n"; 
		    die   "  Recommendation: Verify path from directory containing MIRAGE.pl.\n\n";
		}

		$i++;

	    }
		
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
#                 species.  This will typically be numProcesses, but in
#                 the event that we have some species with fewer gene
#                 families than requested processes we'll reduce this.
#
#             + SpecFamBreaks : The line number that each thread should
#                 count up to (but not hit).
#
sub GenerateSpeciesDBs
{
    my $species_ref    = shift;
    my $ProteinDB_name = shift;
    my $numProcesses   = shift;

    my $DB_base_name = $ProteinDB_name;
    $DB_base_name =~ s/\.[^\/\.]+$//; # remove the file extension.

    my @Species = @{$species_ref};
    my %SpeciesToDBNames;
    foreach my $species (@Species) { $SpeciesToDBNames{$species} = $tempdirname.$species.'.prot.fa'; }
    my $MinorSpeciesDBName = $tempdirname.'.minor-species.prot.fa';
    
    # The very first thing we'll do is scan through our database
    # making a set of species-AND-gene-family mini-databases, and
    # then we can concatenate within species.

    open(my $ProteinDB,'<',$ProteinDB_name) || die "\n  ERROR:  Failed to open input protein database '$ProteinDB_name'\n\n";
    my $spec_gene_filename = 0;
    my %SpeciesGeneFileNames;
    my $SpecGeneFile;
    my $next_entry = 0;
    while (my $line = <$ProteinDB>) {
	
	$line =~ s/\n|\r//g;
	next if (!$line);

	if ($line =~ /\>(\S+)/) {

	    if ($spec_gene_filename) { close($SpecGeneFile); }

	    my $seqname = $1;
	    $seqname =~ /^[^\|]+\|[^\|]+\|([^\|]+)\|\S+\|([^\|]+)$/;
	    my $species = lc($1);
	    my $genefam = lc($2);

	    $spec_gene_filename = $tempdirname.$species.'-'.$genefam.'.fa';
	    $SpeciesGeneFileNames{$species.';'.$genefam} = $spec_gene_filename;
	    
	    open($SpecGeneFile,'>>',$spec_gene_filename) || die "\n  ERROR:  Failed to open output species-gene database '$spec_gene_filename'\n\n";
	    
	    print $SpecGeneFile "$line\n";

	} elsif ($spec_gene_filename) {
	    print $SpecGeneFile "$line\n";
	}

    }
    if ($spec_gene_filename) { close($SpecGeneFile); }
    close($ProteinDB);

    #
    #  Now that we've segrated species and gene families into TONS of little databases,
    #  we can run through our little databases and concatenate them into species databases.
    #
    #  Observe that the primary sorting on species name and secondary sorting on gene family
    #  guarantees that we'll have a really clean re-organization into gene families.
    #

    my %GeneFamsPerSpecies;
    foreach my $spec_gene (sort keys %SpeciesGeneFileNames) {

	$spec_gene =~ /^(\S+)\;/;
	my $species = $1;

	if ($GeneFamsPerSpecies{$species}) { $GeneFamsPerSpecies{$species}++; }
	else                               { $GeneFamsPerSpecies{$species}=1; }

	my $spec_gene_filename = $SpeciesGeneFileNames{$spec_gene};

	my $species_db_name;
	if ($SpeciesToDBNames{$species}) { $species_db_name = $SpeciesToDBNames{$species}; }
	else                             { $species_db_name = $MinorSpeciesDBName;         }

	my $cat_cmd = "cat \"$spec_gene_filename\" >> \"$species_db_name\"";
	if (system($cat_cmd)) { die "\n  ERROR:  Concatenation command '$cat_cmd' failed during execution\n\n"; }

	my $rm_cmd = "rm \"$spec_gene_filename\"";
	if (system($rm_cmd)) { die "\n  ERROR:  Failed to remove species-gene file '$spec_gene_filename'\n\n"; }


    }


    # Next, we can figure out the right order for our species DB names
    my @SpeciesDBNames;
    my @SpeciesFamCounts;
    foreach my $species (@Species) {
	push(@SpeciesDBNames,$SpeciesToDBNames{$species}); 
	push(@SpeciesFamCounts,$GeneFamsPerSpecies{$species});
    }


    # Were there any unindexed species?
    if (-e $MinorSpeciesDBName) { 
	
	# Let's just slip this nasty boy in right quick
	
    }



    #
    #  We can actually give Quilter some cheat codes by identifying the specific lines
    #  at which we want each process to start up.
    #
    #  Note that these are one above the final line, so it's an [i-1,i) sort of thing
    #



    my @SpecFamBreaks;
    my @SpecThreadCounts;
    for (my $i=0; $i<scalar(@SpeciesDBNames); $i++) {

	my $species_db_name = $SpeciesDBNames[$i];
	my $wc_cmd = "wc -l \"$species_db_name\" \|";
	open(my $WC,$wc_cmd) || die "\n  ERROR:  Failed to run linecount command '$wc_cmd'\n\n";
	my $wc_line = <$WC>;
	close($WC);

	my $tot_line_count;
	if ($wc_line =~ /^\s*(\d+)/) {
	    $tot_line_count = $1;
	} else {
	    die "\n  ERROR:  Failed to get a line count for file '$species_db_name'\n\n";
	}

	# Just for good measure, we'll slip in some esl-indexing here
	if (system($eslsfetch." --index \"$species_db_name\" 1>/dev/null")) { 
	    die "\n  ERROR:  Failed to generate esl index for species database '$species_db_name'\n\n"; 
	}
	
	open(my $SpeciesDB,'<',$species_db_name) || die "\n  ERROR:  Failed to open species database '$species_db_name'\n\n";

	# We need some sort of mechanism for adjusting the number of requested
	# cores in the event that there's a species with very few entries.
	#
	if ($SpeciesFamCounts[$i] <= $numProcesses) {

	    #
	    # We'll just assign each core a single gene family
	    #
	    
	    $SpecThreadCounts[$i] = $SpeciesFamCounts[$i];

	    my $linenum = 0;
	    my $lastfam;
	    my $j=0;

	    my $line;
		
	    while ($linenum < $tot_line_count) {
		$line = <$SpeciesDB>;
		$linenum++;
		if ($line =~ /^\>\S+\|([^\|]+)$/) {
		    my $nextfam = lc($1);
		    if ($lastfam && $lastfam ne $nextfam) {
			$SpecFamBreaks[$i][$j] = $linenum;
			$j++;
		    }
		    $lastfam = $nextfam;			
		}
	    }
	    
	    $SpecFamBreaks[$i][$j] = $linenum;

	} else {

	    my $avg_frac = $tot_line_count / $numProcesses;

	    $SpecThreadCounts[$i] = $numProcesses;

	    my $linenum = 0;
	    my $lastfam;
	    my $j=0;
	    while ($j<$SpecThreadCounts[$i]) {
		
		my $line;
		my $fams_in_block = 0;
	    		
		while ($linenum < $tot_line_count && $linenum < ($j+1)*$avg_frac) {
		    $line = <$SpeciesDB>;
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
			$line = <$SpeciesDB>;
			$linenum++;
			if ($line && $line =~ /^\>\S+\|([^\|]+)$/) {
			    my $nextfam = lc($1);
			    if ($nextfam eq $lastfam) {
				$line = <$SpeciesDB>;
				$linenum++;
			    }
			}
		    }
		    
		    $SpecFamBreaks[$i][$j] = $linenum;
		    $j++;

		} else {

		    $SpecThreadCounts[$i]--;

		}
		
	    }

	    if ($SpecThreadCounts[$i]) {
		$SpecFamBreaks[$i][$SpecThreadCounts[$i]-1] = $tot_line_count;
	    } else {
		die "\n  Hey, dude, what's going on with '$species_db_name' and division into threadable blocks?\n\n";
	    }

	}

	close($SpeciesDB);
    
    }


    # This is it, baby!
    return (\@SpeciesDBNames,$MinorSpeciesDBName,\@SpecThreadCounts,\@SpecFamBreaks);

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
# Function Name: AttachSeqToMSA
#
# About: This function attaches a single sequence to an MSA, and is
#        used to handle cases where neither FastDiagonals nor SPALN
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

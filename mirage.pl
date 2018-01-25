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
if (@ARGV == 1) {
    DetailedUsage() if (lc($ARGV[0]) eq '-h' || lc($ARGV[0]) =~ /help/);
    CheckInstall()  if (lc($ARGV[0]) =~ /\-check$/);
}

# Generics
my ($i,$j,$k);


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


# Read in user arguments
my $optsRef = ParseArgs(\@ARGV);
my @Options = @{$optsRef};


# Change options to more intelligible variables
my $ProteinDB    = $Options[0];
my $SpeciesGuide = $Options[1];
my $ResultsDir   = $Options[2];
my $verbose      = $Options[3];
my $numProcesses = $Options[4];
my $timed        = $Options[5];
my $stack_arfs   = $Options[6];
my $forcecompile = $Options[7]; # Hidden
my $cleanMSA     = $Options[8]; # Hidden
my $just_spaln   = $Options[9]; # Hidden


# Verify that we have all the files we need on-hand
CheckSourceFiles($forcecompile);


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


# Make a bunch of species-specific databases
my @SpeciesDBNames;
my @SpeciesDBs;
my %SpeciesToID;
for (my $i=0; $i<$numSpecies; $i++) {
    $SpeciesDBNames[$i] = $ProteinDB;
    $SpeciesDBNames[$i] =~ s/\.fa$/\.Mirage.$Species[$i]\.fa/;
    open($SpeciesDBs[$i],'>',$SpeciesDBNames[$i]);
    $SpeciesToID{lc($Species[$i])} = $i+1;
}
open(my $db_file,'<',$ProteinDB);
my $db_line = <$db_file>;
while (!eof($db_file)) {
    if ($db_line =~ /\>[^\|]+\|[^\|]+\|([^\|]+)\|/) {
	my $seq_spec = lc($1);
	if ($SpeciesToID{$seq_spec}) {
	    my $species_db = $SpeciesDBs[$SpeciesToID{$seq_spec}-1];
	    print $species_db "$db_line";
	    $db_line = <$db_file>;
	    while ($db_line !~ /\>/) {
		print $species_db "$db_line";
		last if (eof($db_file));
		$db_line = <$db_file>;
	    }
	} else {
	    $db_line = <$db_file>;
	}
    } else {
	$db_line = <$db_file>;
    }
}
for (my $i=0; $i<$numSpecies; $i++) { 
    close($SpeciesDBs[$i]);
    system("esl-sfetch --index \"$SpeciesDBNames[$i]\" 1>/dev/null 2>\&1");
}


# Run Quilter
my %MultiMSADir;
my %GroupsBySpecies;
foreach $i (0..$numSpecies-1) {

    # Start a timer for Quilter
    $IntervalStart = [Time::HiRes::gettimeofday()];

    # Assembling the Quilter command.  Because we're a directory back
    # from src, we need to append a '../' to each of the file names.
    my $QuilterCmd;
    $QuilterCmd = 'perl src/Quilter.pl '; # base call

    # Depending on whether we're using a relative path (w.r.t. the current
    # directory) or absolute path, what we pass to Quilter will vary for
    # each of the files

    # Protein database
    my $SpeciesDB = $SpeciesDBNames[$i];
    if ($SpeciesDB =~ /^\//) { $QuilterCmd = $QuilterCmd.$SpeciesDB.' ';       }
    else                     { $QuilterCmd = $QuilterCmd.'../'.$SpeciesDB.' '; }

    # Genome
    if ($Genomes[$i] =~ /^\//) { $QuilterCmd = $QuilterCmd.$Genomes[$i].' ';       }
    else                       { $QuilterCmd = $QuilterCmd.'../'.$Genomes[$i].' '; }

    # GTF index (special case: - for "Just use BLAST")
    if ($GTFs[$i] =~ /^\// || $GTFs[$i] eq '-') { $QuilterCmd = $QuilterCmd.$GTFs[$i].' ';       }
    else                                        { $QuilterCmd = $QuilterCmd.'../'.$GTFs[$i].' '; }
	
    # Species name
    $QuilterCmd = $QuilterCmd.$Species[$i];
    
    # Guide output towards results directory
    $QuilterCmd = $QuilterCmd.' -o ../'.$SpeciesDir{$Species[$i]};

    # Set Quilter to work in the src directory, since it needs to interact with
    # other tools there.
    $QuilterCmd = $QuilterCmd.' -setcwd '.cwd().'/src';
    
    # Tell Quilter how many CPUs we want it to work with
    $QuilterCmd = $QuilterCmd.' -n '.$numProcesses;

    # Do we want a whole bunch of stuff spat at us?
    $QuilterCmd = $QuilterCmd.' -v' if ($verbose);

    # Are we timing?
    $QuilterCmd = $QuilterCmd.' -time' if ($timed);

    # Are we only using SPALN?
    $QuilterCmd = $QuilterCmd.' -fast' if ($just_spaln);

    # Display the call if we're being loud
    if ($verbose) {
	print "\n  Aligning to the $Species[$i] genome using Quilter.pl\n";
	print "  $QuilterCmd\n";
    }

    # Make that daaaaaang call.
    if (system($QuilterCmd)) {
	# Rats!
	foreach my $species_db_name (@SpeciesDBNames) {
	    if (-e $species_db_name) {
		system("rm \"$species_db_name\"");
		system("rm \"$species_db_name\.ssi\"");
	    }
	}
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
    $MultiMSACmd = 'perl src/MultiMSA.pl ';       # base call
    $MultiMSACmd = $MultiMSACmd.$HitFileName.' '; # Quilter output
    $MultiMSACmd = $MultiMSACmd.$SpeciesDB;       # isoform file
    
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
	foreach my $species_db_name (@SpeciesDBNames) {
	    if (-e $species_db_name) {
		system("rm \"$species_db_name\"");
		system("rm \"$species_db_name\.ssi\"");
	    }
	}
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
    my $singleseqfilename = 'Mirage.Single.Temp.fa';
    my $singleseqresults  = 'Mirage.Multi.Temp.afa';

    foreach my $missedgroup (keys %QuilterMisses) {
	my @seqnames = split('&',$QuilterMisses{$missedgroup});
	foreach my $seq (@seqnames) {
	    my $GBSref = AttachSeqToMSA($seq,$MultiMSADir{$Species[$i]},$missedgroup,
					$SpeciesDB,\%GroupsBySpecies,$i+1,
					$singleseqfilename,$singleseqresults);
	    %GroupsBySpecies = %{$GBSref};
	}
    }
    system("rm $singleseqfilename") if (-e $singleseqfilename);
    system("rm $singleseqresults")  if (-e $singleseqresults);

    # Knock it off with that darn timing!
    $MultiMSATimeStats[$i] = Time::HiRes::tv_interval($IntervalStart);

    # Get rid of that stinky old pair of files!
    system("rm \"$SpeciesDB\"");
    system("rm \"$SpeciesDB\.ssi\"");
    
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
if ($TotNumGroupFiles < $numProcesses) { $numProcesses = $TotNumGroupFiles; }


# Carve out a portion of the files for each process to work on,
# or else die.  This would be really strange to see happen, but
# it's good to cover your bases.
my $portion;
if ($numProcesses) {
    $portion = $TotNumGroupFiles/$numProcesses;
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
my $progressbase = 'mirage.thread_progress.';
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
my $tempfileA = 'tempfileA.'.$ThreadID.'.MIRAGE.afa';
my $tempfileB = 'tempfileB.'.$ThreadID.'.MIRAGE.afa';
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
		my $progressCmd = './src/ProgressTimer.pl '.$progressbase.' '.$numProcesses;
		$progressCmd = $progressCmd.' '.$num_complete.' |';
		
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
	    my $FinalMSACmd = "perl src/FinalMSA.pl '".$MultiMSADir{$GroupSpecies[0]}."/".$groupfile."'";
	    $FinalMSACmd    = $FinalMSACmd." '".$ResultFile."'";
	    
	    # Print the command and execute
	    #print "  $FinalMSACmd\n" if ($verbose);
	    if (system($FinalMSACmd)) { die "\n  *  ERROR: FinalMSA failed during execution  *\n\n"; }
	    
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

    my $MultiSeqNWCmd;
    $MultiSeqNWCmd = "./src/MultiSeqNW '$fileA' $fileAsize"; # First species
    $MultiSeqNWCmd = $MultiSeqNWCmd." '$fileB' $fileBsize";  # Second species
    $MultiSeqNWCmd = $MultiSeqNWCmd.' > '.$tempfileA;        # Send to temp file    

    print "  $MultiSeqNWCmd\n" if ($verbose);
    if (system($MultiSeqNWCmd)) { die "\n  *  ERROR: MultiSeqNW failed during execution of command '$MultiSeqNWCmd'  *\n\n"; }
    
    $i = 2;
    while ($i < $speciesPerGroup) {

	if ($i % 2) { $fileA = $tempfileB; }
	else        { $fileA = $tempfileA; }

	$fileAsize += $fileBsize;
	$fileB      = $MultiMSADir{$GroupSpecies[$i]}.'/'.$groupfile;
	$fileBsize  = $GroupSeqs[$i];
	$MultiSeqNWCmd = "./src/MultiSeqNW '$fileA' $fileAsize"; # First species
	$MultiSeqNWCmd = $MultiSeqNWCmd." '$fileB' $fileBsize";  # Second species

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

	my $FinalMSACmd = 'perl src/FinalMSA.pl ';
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
if (grep -f, glob 'src/*.TEs.Quilter.out') {
    system("mkdir $ResultsDir/SuspectedTEs");
    system("mv src/*.TEs.Quilter.out $ResultsDir/SuspectedTEs");
}


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
    print " USAGE  : perl mirage.pl  <Isoform DB>  <Species Guide>  [OPT.s]                        \n\n";
    print " ARG.s  : <IsoformDB>     : A FASTA-formatted protein database, with the following      \n";
    print "                            naming convention: > gene_name|prot_name|species|ID         \n\n";
    print "          <Species Guide> : A file indicating, for each species being searched on,      \n";
    print "                            the location of a FASTA-formatted genome for that species   \n";
    print "                            and a .gtf index file corresponding to that species and     \n";
    print "                            genome.  These fields should be comma-separated and         \n";
    print "                            ordered as follows:  species, genome, .gtf index            \n";
    print "                            It is recommended that similar species are grouped together \n";
    print "                            and positioned near the top of the list.                    \n\n";
    print " OPT.s  : -h     : More detailed help.                                                  \n";
    print "          -o     : Specify output directory name.                                       \n";
    print "          -v     : Verbose output.                                                      \n";
    print "          -n     : Specify number of CPU cores (default: 2)                             \n";
    print "          --time : Print timing data to stdout at end of program                        \n";
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
    print "    USAGE :  perl  mirage.pl  <Isoform DB>  <Species Guide>  [OPT.s]              \n";
    print "                                                                                  \n";
    print "                                                                                  \n";
    print "    ARG.s :  Isoform DB    : A FASTA-formatted protein database using the         \n";
    print "                             following sequence naming convention:                \n";
    print "                                                                                  \n";
    print "                             >GN:gene_name|protein_name|species|seqID|groupID     \n";
    print "                                                                                  \n";
    print "             Species Guide : A simple file indicating, for each species being     \n";
    print "                             searched on, the location of a FASTA-formatted       \n";
    print "                             genome for that species and a .gtf index file        \n";
    print "                             corresponding to that species and genome.            \n";
    print "                                                                                  \n";
    print "                             The filenames should be paths to the files from      \n";
    print "                             the directory containing MIRAGE.pl and must be       \n";
    print "                             separated by commas.                                 \n";
    print "                                                                                  \n";
    print "                             Required order:   species, genome, index file        \n";
    print "                                                                                  \n";
    print "                               EXAMPLE:                                           \n";
    print "                             .------------------------------------------------.   \n";
    print "                             | human, ~/data/HumanGenome.fa, ~/data/human.gtf |   \n";
    print "                             | mouse, ~/data/MouseGenome.fa, ~/data/mouse.gtf |   \n";
    print "                             | ...                                            |   \n";
    print "                             '------------------------------------------------'   \n";
    print "                                                                                  \n";
    print "                             Species are introduced into the final genewise MSAs  \n";
    print "                             in the order that they are listed in the species     \n";
    print "                             guide.  It is recommended that similar species are   \n";
    print "                             grouped together for optimal alignments.             \n";
    print "                                                                                  \n";
    print "                                                                                  \n";
    print "    OPT.s :  -o <string> : Specify ouptut directory name                          \n";
    print "             -v          : Verbose output                                         \n";
    print "             -n <int>    : Specify number of CPU cores (default: 2)               \n";
    print "             --time      : Print timing data to stdout at end of program          \n";
    print "             -h          : Print this help page                                   \n";
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
    open(my $cmdcheck,'which spaln |');
    my $line = <$cmdcheck>;
    close($cmdcheck);
    if (!$line) { die "\n  Failure: spaln does not appear to be on the PATH\n\n"; }    

    open($cmdcheck,'which esl-sfetch |');
    $line = <$cmdcheck>;
    close($cmdcheck);
    if (!$line) { die "\n  Failure: esl-sfetch does not appear to be on the PATH\n\n"; }

    open($cmdcheck,'which esl-seqstat |');
    $line = <$cmdcheck>;
    close($cmdcheck);
    if (!$line) { die "\n  Failure: esl-seqstat does not appear to be on the PATH\n\n"; }

    my @RequiredFiles;
    push(@RequiredFiles,'src/Quilter.pl');
    push(@RequiredFiles,'src/DiagonalSets.pm');
    push(@RequiredFiles,'src/FastDiagonals.c');
    push(@RequiredFiles,'src/Diagonals.c');
    push(@RequiredFiles,'src/Diagonals.h');
    push(@RequiredFiles,'src/TransSW.c');
    push(@RequiredFiles,'src/MultiMSA.pl');
    push(@RequiredFiles,'src/MultiSeqNW.c');
    push(@RequiredFiles,'src/MultiSeqNW.h');
    push(@RequiredFiles,'src/FinalMSA.pl');
    push(@RequiredFiles,'src/makefile');

    foreach my $file (@RequiredFiles) {
	if (!(-e $file)) { die "\n  Failure: Could not locate required file '$file'\n\n"; }
    }

    open($cmdcheck,'which blat |');
    $line = <$cmdcheck>;
    close($cmdcheck);
    if (!$line) {
	die "\n  Success! (Although installing blat is strongly recommended)\n\n";
    } else {
	die "\n  Success!\n\n";
    }
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
    my $forcecompile = shift;

    my @RequiredFiles;
    push(@RequiredFiles,'src/Quilter.pl');
    push(@RequiredFiles,'src/DiagonalSets.pm');
    push(@RequiredFiles,'src/FastDiagonals.c');
    push(@RequiredFiles,'src/Diagonals.c');
    push(@RequiredFiles,'src/Diagonals.h');
    push(@RequiredFiles,'src/TransSW.c');
    push(@RequiredFiles,'src/MultiMSA.pl');
    push(@RequiredFiles,'src/MultiSeqNW.c');
    push(@RequiredFiles,'src/MultiSeqNW.h');
    push(@RequiredFiles,'src/FinalMSA.pl');
    push(@RequiredFiles,'src/makefile');

    foreach my $srcfile (@RequiredFiles) {
	if (!(-e $srcfile)) {
	    die "\n  Failed to locate required file '$srcfile'\n\n";
	}
    }

    my @RequiredPrograms;
    push(@RequiredPrograms,'src/FastDiagonals');
    push(@RequiredPrograms,'src/TransSW');
    push(@RequiredPrograms,'src/MultiSeqNW');

    my $missingProgram = 0;
    foreach my $srcprogram (@RequiredPrograms) {
	if (!(-e $srcprogram) || $forcecompile) {
	    $missingProgram = 1;
	}
    }

    if ($missingProgram) {
	chdir('./src');
	print "\n  Compiling from src/makefile\n";
	if (system("make --silent")) {
	    die "  ERROR: Compilation failed.\n\n";
	}
	chdir('..');
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
    my $argsRef = shift;
    my @ARGS    = @{$argsRef};

    my @Options;
    $Options[0] = $ARGS[0];        # Isoform file
    $Options[1] = $ARGS[1];        # Species guide
    $Options[2] = 'MirageResults'; # Output folder name
    $Options[3] = 0;               # Verbose output
    $Options[4] = 2;               # Number of CPUs
    $Options[5] = 0;               # Report timing information with final output
    $Options[6] = 0;               # Stack ARFs with the rest of the sequence
    $Options[7] = 0;               # (hidden) Force Compile
    $Options[8] = 1;               # (hidden) perform final cleanup
    $Options[9] = 0;               # (hidden) Only use SPALN (no FastDiagonals)

    my $i = 2;
    while ($i < @ARGS) {
	if ($ARGS[$i] eq '-h') {
	    DetailedUsage(); # Kills program
	} elsif ($ARGS[$i] eq '-o') {
	    $i++;
	    $Options[2] = $ARGS[$i];
	    $Options[2] =~ s/\/$//; # Clip terminal '/'
	} elsif ($ARGS[$i] eq '-v') {
	    $Options[3] = 1;
	} elsif ($ARGS[$i] eq '-n') {
	    $i++;
	    $Options[4] = int($ARGS[$i]);
	    if ($Options[4] <= 0) {
		print "  Unsupported number of CPUs requested ($Options[4])\n";
		print "  Reverting to default (2)\n.";
		$Options[4] = 2;
	    }
	} elsif ($ARGS[$i] eq '--time') {
	    $Options[5] = 1;
	} elsif ($ARGS[$i] eq '--stack-arfs') {
	    $Options[6] = 1;
	} elsif ($ARGS[$i] eq '--forcecompile') {
	    $Options[7] = 1;
	} elsif ($ARGS[$i] eq '--showSpliceSites') {
	    $Options[8] = 0;
	} elsif ($ARGS[$i] eq '--fast') {
	    $Options[9] = 1;
	} else {
	    print "  Unrecognized option '$ARGS[$i]' ignored\n";
	}
	$i++;
    }

    # While we're here, we make sure there's a .ssi index for the
    # protein database
    if (!(-e $ARGS[0].'.ssi')) {
	if (system("esl\-sfetch \-\-index $ARGS[0]")) { 
	    die "\n  Failed to create easel index for $ARGS[0]\n\n"; 
	}
    }

    return \@Options;
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
		    if (system("esl-sfetch --index $Genomes[$i]")) { 
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
    
    # Grab the name of the protein database
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
		my $NWcmd = './src/MultiSeqNW '.$tempfilename.' 1 '."'$speciesDir/$species/$groupID.afa' ";
		$NWcmd    = $NWcmd.$MinorGroupsBySpecies{$groupID.'#'.$species};
		$NWcmd    = $NWcmd.' > '.$tempresults;
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
    my $eslsfetchCmd = "esl-sfetch -o \"$seqfilename\" \"$ProteinDB\" \"$origseq\" > /dev/null 2>\&1";
    
    
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
		my $alignCmd = "./src/MultiSeqNW '".$seqfilename."' 1 '".$baseMSAfile."' ".$num_isos;
		$alignCmd    = $alignCmd.' -igBase 0';
		$alignCmd    = $alignCmd." > '".$resultsname."'";
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

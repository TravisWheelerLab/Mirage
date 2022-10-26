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
sub GetScriptDir { return '.' if ($0 !~ /\//); $0 =~ /^(.+)\/[^\/]+$/; return $1; }
use lib GetScriptDir();
use BureaucracyMirage;
use DisplayProgress;

sub PrintUsage;
sub DetailedUsage;
sub PrintVersion;
sub CheckSourceFiles;
sub ParseArgs;
sub VerifiedClean;
sub ParseSpeciesGuide;
sub SetMergeOrder;
sub ConfirmValidTree;
sub GenerateSpeciesDBs;
sub ParseSeqNameAsMirage;
sub ParseSeqNameAsUniProt;
sub SaveSpeciesMappings;
sub AggregateMappingMisses;
sub AlignUnmappedSeqs;
sub AlignMiscSeqs;
sub MergeAlignments;
sub FinalizeIntraSpeciesMSA;
sub ReorganizeResultsForMapping;
sub EvaluateMissDir;
sub ReportGranularQuilterTiming;
sub PrintTimingInfo;
sub FormatTimeString;


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


# Figure out what the location of the Mirage src directory is
my $location = $0;
$location =~ s/Mirage2\.pl$//;

# Where are all our dependencies?
my $dependencies_ref = FindDependencies();
my %Dependencies = %{$dependencies_ref};

# We're going to need these friends
my $sindex       = $Dependencies{'sindex'};
my $sfetch       = $Dependencies{'sfetch'};
my $sstat        = $Dependencies{'sstat'};
my $quilter      = $Dependencies{'quilter2'};
my $maps_to_msas = $Dependencies{'mapstomsas'};
my $multi_seq_nw = $Dependencies{'multiseqnw'};
my $final_msa    = $Dependencies{'finalmsa'};


# Read in user arguments
my $optsRef = ParseArgs();
my %Options = %{$optsRef};


# Change options to more intelligible variables
my $ProteinDB       = ConfirmFile($ARGV[0]);
my $SpeciesGuide    = ConfirmFile($ARGV[1]);
my $ResultsDir      = $Options{outdirname};
my $verbose         = $Options{verbose};
my $num_cpus        = $Options{cpus};
my $timing          = $Options{time};
my $only_map        = $Options{only_map};
my $blat_off        = $Options{blat_off};
my $stack_arfs      = $Options{stackarfs};     # Hidden
my $forcecompile    = $Options{forcecompile};  # Hidden
my $cleanMSA        = $Options{cleanmsa};      # Hidden
my $just_spaln      = $Options{justspaln};     # Hidden
my $track_spaln     = $Options{trackspaln};    # Hidden
my $gene_timing     = $Options{genetiming};    # Hidden
my $max_spaln_nucls = $Options{maxspalnnucls}; # Hidden


# Things have officially kicked off, so start timing
my $MirageTimer = StartTimer();

# Verify that we have all the files we need on-hand
#CheckSourceFiles($forcecompile);

# Make a directory to store results
$ResultsDir = CreateDirectory($ResultsDir);

# Stick in a directory to record progress information as we're running the
# various tools that we're so happy to have this opportunity to be running!
my $progress_dirname = CreateDirectory($ResultsDir.'.progress');
InitMirageProgressVars($progress_dirname,$num_cpus);

# We'll also need a folder to record any mapping errors
my $misses_dirname = CreateDirectory($ResultsDir.'Mapping-Misses');

# Convert species guide into three arrays.
# We'll ignore any species that don't show up
# in our protein database.
my ($speciesref,$gtfsref,$genomesref,$mergeorderref)
    = ParseSpeciesGuide($SpeciesGuide,$ProteinDB);
my @Species     = @{$speciesref};
my @GTFs        = @{$gtfsref};
my @Genomes     = @{$genomesref};
my @MergeOrder  = @{$mergeorderref};
my $num_species = scalar(@Species);

# Make species-specific results folders
my %SpeciesDir;
my $AllSpeciesDir = CreateDirectory($ResultsDir.'Species-MSAs');
foreach my $species (@Species) {
    $SpeciesDir{$species} = CreateDirectory($AllSpeciesDir.$species);
    CreateDirectory($SpeciesDir{$species}.'seqs');
    CreateDirectory($SpeciesDir{$species}.'alignments');
    CreateDirectory($SpeciesDir{$species}.'mappings');
    CreateDirectory($SpeciesDir{$species}.'timing') if ($gene_timing);
}

# Divide the database according to species
DispProgMirage('db-speciation');
my ($origseqnames_ref,$allgenes_ref,$misc_seqs,$speciesdir_ref)
    = GenerateSpeciesDBs($ProteinDB,$num_cpus,\%SpeciesDir);
my @OrigSeqNames = @{$origseqnames_ref};
my @AllGenes = @{$allgenes_ref};
%SpeciesDir = %{$speciesdir_ref};

# I'm going to print off all the original sequence names to a file so
# if things go wrong during a run we can actually figure out who's who
my $seqnamefname = $ResultsDir.'seq-name-guide';
my $seqnamef = OpenOutputFile($seqnamefname);
for (my $i=0; $i<scalar(@OrigSeqNames); $i++) {
    print $seqnamef "$i: $OrigSeqNames[$i]\n";
}
close($seqnamef);


# Do our intra-species magic!
for (my $i=0; $i<$num_species-1; $i++) {

    # Did we discover that this species wasn't actually present?
    next if (!$SpeciesDir{$Species[$i]});

    my $SpeciesTimer = StartTimer();

    my $species_dirname = $SpeciesDir{$Species[$i]};
    my $species_seqdir  = $species_dirname.'seqs/';

    #
    #  Q U I L T E R
    #

    # Start assembling that command!
    my $QuilterCmd = $quilter;

    # OPT: Collect timing data, with varying granularity
    $QuilterCmd = $QuilterCmd.' --time' if ($timing);
    $QuilterCmd = $QuilterCmd.' --genetiming' if ($gene_timing);

    # OPT: Turn off BLAT for faster runtime
    $QuilterCmd = $QuilterCmd.' --blatoff' if ($blat_off);

    # OPT: Manually set cap on nucleotide sequence length for spaln search
    $QuilterCmd = $QuilterCmd.' -maxspalnnucls='.$max_spaln_nucls
	if ($max_spaln_nucls);
	
    # KEY 1: Protein database (implicit in species directory, under 'seqs/')
    $QuilterCmd = $QuilterCmd.' '.$species_dirname;

    # KEY 2: Genome
    $QuilterCmd = $QuilterCmd.' '.$Genomes[$i];

    # KEY 3: GTF index (special case: - for "Just use BLAT")
    $QuilterCmd = $QuilterCmd.' '.$GTFs[$i];

    # Start timing Quilter specifically
    my $QuilterTimer = StartTimer();

    # Make that daaaaaang call.
    RunSystemCommand($QuilterCmd);

    # If we're timing, let's see how long Quilter took
    if ($timing) {

	my $quilter_timing_str
	    = SecondsToSMHDString(GetElapsedTime($QuilterTimer));
	ClearProgress();
	print "\n";
	print "  $Species[$i]\n";
	print "    - Trans. Mapping   : $quilter_timing_str\n";

	my $quilter_timer_fname = $species_dirname.'quilter-timing.out';
	ReportGranularQuilterTiming($quilter_timer_fname);
	RunSystemCommand("rm $quilter_timer_fname");

    }


    #
    #  M A P s   t o   M S A s
    #

    # Now we want to check how long it takes to run MapsToMSAs
    my $MapsToMSAsTimer = StartTimer();

    # Construct the command (typically just the executable and the sequence directory).
    my $MapsToMSAsCmd = $maps_to_msas;
    $MapsToMSAsCmd = $MapsToMSAsCmd.' --genetiming' if ($gene_timing);
    $MapsToMSAsCmd = $MapsToMSAsCmd.' '.$species_seqdir;

    # Rock 'n' roll 'n' align!
    RunSystemCommand($MapsToMSAsCmd);


    # NOTE: Because some key data management the takes place in
    # MapsToMSAs, we'll run it even in the case that we're just
    # generating our mappings.  Luckily, this is a really fast step
    # in the program, so we shouldn't be hurting on account of the
    # unnecessary extra work...
    if ($only_map) {
	ClearProgress();
	if ($timing) {
	    my $species_timing_str
		= SecondsToSMHDString(GetElapsedTime($SpeciesTimer));
	    print "\n  + Mapping Time for $Species[$i] : $species_timing_str\n";
	} else {
	    print "  Mapping complete for $Species[$i]\n";
	}
	next;
    }
    

    # Knock it off with that darn timing!
    if ($timing) {
	my $maps2msas_timing_str
	    = SecondsToSMHDString(GetElapsedTime($MapsToMSAsTimer));
	ClearProgress();
	print "    - Mappings -> MSAs : $maps2msas_timing_str\n";
    }

    
    #  M u l t i S e q N W

    # Time NW-style alignment (intra-species)
    my $NWTimer = StartTimer();

    # Align any sequences that didn't fully map by getting all Needleman-Wunsch-y
    my $missesbygene_ref
	= AggregateMappingMisses($species_dirname,$misses_dirname,\@OrigSeqNames);
    AlignUnmappedSeqs($missesbygene_ref,$species_seqdir);

    # #SaveTheMappings #MappingsArePeopleToo!
    SaveSpeciesMappings($species_dirname,\@OrigSeqNames);

    # How'd we do?
    if ($timing) {
	my $nw_timing_str
	    = SecondsToSMHDString(GetElapsedTime($NWTimer));
	ClearProgress();
	print "    - Unmapped Seq Ali : $nw_timing_str\n";
    }


    # Species[i] over and out!
    ClearProgress();
    print "\n" if ($timing);
    print "  Intra-species alignment complete for $Species[$i]\n";

    if ($timing) {
	my $species_timing_str
	    = SecondsToSMHDString(GetElapsedTime($SpeciesTimer));
	print "  + Total Intraspecies Alignment Time : $species_timing_str\n";
    }

}



# Are we only mapping?
if ($only_map) {

    # We'll need to be able to associate sequence IDs with names
    my $SeqNameFile = OpenInputFile($seqnamefname);
    my %SeqIDsToNames;
    while (my $line = <$SeqNameFile>) {
	if ($line =~ /^(\d+)\:\s*(\S.+\S)\s*$/) {
	    my $seq_id = $1;
	    my $seq_name = $2;
	    $SeqIDsToNames{$seq_id} = $seq_name;
	}
    }
    close($SeqNameFile);

    # Make our mapping data output directory
    my $mapdirname = CreateDirectory($ResultsDir.'Mappings-by-Gene');

    # Move things over by species
    for (my $species_id=0; $species_id < $num_species-1; $species_id++) {

	my $species = $Species[$species_id];

	next if (!$SpeciesDir{$species});
	next if (!(-d $SpeciesDir{$species}.'seqs'));

	ReorganizeResultsForMapping($species,$SpeciesDir{$species}.'seqs/',$mapdirname,
				    \%SeqIDsToNames);

    }

    # No more need for 'Species-MSAs'!
    RunSystemCommand("rm -rf \"$AllSpeciesDir\" \&");

} else {

    # NOPE! We're aligning, dude!

    # But don't you dare forget about the 'Misc' sequences...
    if ($misc_seqs) {

	my $MiscTimer = StartTimer();

	AlignMiscSeqs($SpeciesDir{'Misc'});
	
	if ($timing) {
	    my $misc_timing_str
		= SecondsToSMHDString(GetElapsedTime($MiscTimer));
	    print "  Alignment of Species without Genomes: $misc_timing_str\n";
	}

    }

    
    # We're now going to track how long the whole final-MSA-generating
    # part of the program takes
    my $AliMergeTimer = StartTimer();

    # Perform a most unnatural merging of alignments! (interspecies -- scandalous!)
    MergeAlignments(\@Species,\%SpeciesDir,\@MergeOrder,\@AllGenes,\@OrigSeqNames);

    # Slap that stop-watch!
    if ($timing) {
	my $ali_merge_timing_str = SecondsToSMHDString(GetElapsedTime($AliMergeTimer));
	ClearProgress();
	print "\n  Interspecies Alignment: $ali_merge_timing_str\n";
    }

    
    # Even though there'll be a few moments for cleanup, who says we can't pop
    # some champagne bottles?
    print "  ";
    print "+ " if ($timing);
    print "Inter-species alignment complete\n";
    print "\n" if ($timing);

    
    # If there weren't any genomeless sequences, we can clear out the Misc dir
    if (!$misc_seqs) {
	my $miscdir = $SpeciesDir{'Misc'};
	RunSystemCommand("rm -rf \"$miscdir\" \&");
    }

}



# Time for the end-of-execution mop-up!


# For a touch of cleanup, get rid of our species-specific protein databases
for (my $i=0; $i<$num_species; $i++) {
    my $seqdirname = $SpeciesDir{$Species[$i]}.'seqs/';
    RunSystemCommand("rm -rf \"$seqdirname\"");
}


# We'll also want to clear out the alias and name guide files
RunSystemCommand("rm \"$seqnamefname\"");
my $aliasfname = $ResultsDir.'gene-aliases';
RunSystemCommand("rm \"$aliasfname\"") if (-e $aliasfname);


# Check whether we actually had any mapping misses to report -- if not,
# celebrate by burning the miss directory to the ground
EvaluateMissDir($misses_dirname);


# No more progress to be made!
system("rm -rf \"$progress_dirname\" \&");


# Last chance to get timing data
if ($timing) {
    my $total_runtime_str = SecondsToSMHDString(GetElapsedTime($MirageTimer));
    ClearProgress();
    print "  Total Runtime  : $total_runtime_str\n";
}


# WE DID IT!
ClearProgress();
print "\n  Mirage complete; results in $ResultsDir\n\n";


####  NOICE!
1; #  NOICE!
####  NOICE!







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
    print " Mirage2: Multiple-sequence Isoform Alignment Tool Guided by Exon Boundaries ($mirage_version)\n";
    print "\n";
    print " USAGE  : mirage [OPT.s] <Isoform DB>  <Species Guide>                                  \n";
    print "\n";
    print " ARG.s  : <IsoformDB>     : A FASTA-formatted protein database, with the following      \n";
    print "                            naming convention:                                          \n";
    print "\n";
    print "                            >gene_name|protein_name|species|seqID|groupID               \n";
    print "\n";
    print "          <Species Guide> : A file indicating, for each species being searched on,      \n";
    print "                            the location of a FASTA-formatted genome for that species   \n";
    print "                            and a .gtf index file corresponding to that species and     \n";
    print "                            genome.  These fields should be whitespace-separated and    \n";
    print "                            ordered as follows:  species, genome, .gtf index            \n";
    print "                            It is recommended that similar species are grouped together \n";
    print "                            and positioned near the top of the list.                    \n";
    print "\n";
    print " OPT.s  : --help      : More detailed help.                                             \n";
    print "          --verbose   : Verbose output.                                                 \n";
    print "          --time      : Print timing data to stdout at end of program                   \n";
    print "          --only_map  : Stop after producing protein-to-genome mappings                 \n";
    print "          --blat_off  : Prevents BLAT search (faster, may miss some mappings)           \n";
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
    print " .-=======-+--=-----------------=----=--=--------------=---------=-------------------.  \n";
    print " | Mirage2 :  Multiple-sequence Isoform Alignment Tool Guided by Exon Boundaries (2) |  \n";
    print " '-=======-+--=-----------------=----=--=--------------=---------=-------------------'  \n";
    print "  $mirage_version\n";
    print "                                                                                   \n";
    print "    USAGE :  mirage  [OPT.s]  <Isoform DB>  <Species Guide>                        \n";
    print "                                                                                   \n";
    print "                                                                                   \n";
    print "    ARG.s :  Isoform DB    : A FASTA-formatted protein database using the          \n";
    print "                             following sequence naming convention:                 \n";
    print "                                                                                   \n";
    print "                             >gene_name|protein_name|species|seqID|groupID         \n";
    print "                                                                                   \n";
    print "             Species Guide : A simple file indicating, for each species being      \n";
    print "                             searched on, the location of a FASTA-formatted        \n";
    print "                             genome for that species and a .gtf index file         \n";
    print "                             corresponding to that species and genome.             \n";
    print "                                                                                   \n";
    print "                             The filenames should be paths to the files from       \n";
    print "                             the directory containing MIRAGE.pl and must be        \n";
    print "                             separated by whitespace.                              \n";
    print "                                                                                   \n";
    print "                             Required order:   species, genome, index file         \n";
    print "                                                                                   \n";
    print "                               EXAMPLE:                                            \n";
    print "                             .------------------------------------------------.    \n";
    print "                             | human  ~/data/HumanGenome.fa  ~/data/human.gtf |    \n";
    print "                             | mouse  ~/data/MouseGenome.fa  ~/data/mouse.gtf |    \n";
    print "                             | ...                                            |    \n";
    print "                             '------------------------------------------------'    \n";
    print "                                                                                   \n";
    print "                             Species are introduced into the final genewise MSAs   \n";
    print "                             in the order that they are listed in the species      \n";
    print "                             guide.  It is recommended that similar species are    \n";
    print "                             grouped together for optimal alignments.              \n";
    print "                                                                                   \n";
    print "                                                                                   \n";
    print "    OPT.s :  --verbose            : Verbose output                                 \n";
    print "             --time               : Print timing data to stdout at end of program  \n";
    print "             --only_map           : Stop after producing protein-to-genome mappings\n";
    print "             --blat_off           : Prevents BLAT search (faster, may miss some mappings)\n";
    print "             -outdirname <string> : Specify ouptut directory name                  \n";
    print "             -cpus <int>          : Specify number of CPU cores (default: 2)       \n";
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
# Function Name: ParseArgs
#
# About: This function takes the command line arguments provided by
#        the user and makes sense of them.
#
sub ParseArgs
{

    my %Options = ( 
	cpus => 2,
	outdirname => 'Mirage-Results',
	cleanmsa => 1,
	maxspalnnucls => 0,
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
	"only_map",
	"blat_off",
	"stackarfs",       # Hidden
	"forcecompile",    # Hidden
	"cleanmsa=i",      # Hidden
	"justspaln",       # Hidden
	"trackspaln",      # Hidden
	"genetiming",      # Hidden
	"maxspalnnucls=i", # Hidden
	)
	|| die "\n  ERROR:  Failed to parse command line arguments\n\n";

    # If the user just wants a little help, I don't think it's too
    # difficult to give them a hand
    if ($Options{help})    { DetailedUsage();             }
    if ($Options{check})   { die "\n  Looking good!\n\n"; } # FindDependencies already passed
    if ($Options{version}) { PrintVersion();              }

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
    # UPDATE: We've switched to .hsi, which will automatically check if this
    #         is a good index, so just run that
    RunSystemCommand($sindex." \"$proteindb\"");

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
	if ($header =~ /\s|\#|\&|\/|\$/) {
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
#        sure everything has the right (.hsi) indexing structures.
#
sub ParseSpeciesGuide
{
    
    my $guide_fname  = shift;
    my $protdb_fname = shift;

    my $GuideFile = OpenInputFile($guide_fname);

    # It's possible that we'll need to know what '~' means
    my $homecheck = OpenSystemCommand("echo ~");
    my $home = <$homecheck>;
    close($homecheck);
    $home =~ s/\n|\r//g;
    $home =~ s/\/$//g;
    
    my $i = 0;
    my (@Species,@GTFs,@Genomes);
    my $species_tree = 0;
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
	    
	    # Make sure the genome has an .hsi index
	    RunSystemCommand($sindex." \"$Genomes[$i]\"");

	    # If we don't have a '.chrom.sizes' file (std name for UCSC) then
	    # generate one
	    my $chrsizes = $Genomes[$i];
	    $chrsizes =~ s/\.fa$/\.chrom.sizes/;
	    if (!(-e $chrsizes)) {
		my $ChrSizes = OpenOutputFile($chrsizes);
		my $SeqStats = OpenSystemCommand($sstat." \"$Genomes[$i]\"");
		while (my $statline = <$SeqStats>) {
		    if ($statline =~ /^\:/) {
			$statline =~ /^\:\s+(\S+)/;
			my $chr = $1;
			$statline =~ /(\d+)$/;
			my $len = $1;
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

	} elsif ($line =~ /^\s*(\(.*\))\s*$/) {

	    # This looks like a species tree string to me!
	    $species_tree = lc($1);
	    $species_tree =~ s/\s//g;
	    
	} else {
	    
	    # Something isn't right
	    die "\n  Species guide line '$line' does not match expected format\n\n";
	    	    
	}

    }
    
    close($GuideFile);

    # If we didn't have any species, this is a prank species guide
    die "\n  No species entries found in '$guide_fname'\n\n" if (scalar(@Species) == 0);

    # We'll also create a directory to house any sequences belonging to species that
    # we don't have genomes for.  'Misc' is safe to use because all actual species are
    # lower-cased in our list.
    push(@Species,'Misc');

    # Now we'll generate a list of merges that will represent our species tree.
    # The list contents will be indices into Species or the Merge list, indicated
    # by 'S' or 'M'
    my $mergeorder_ref = SetMergeOrder($species_tree,\@Species);

    return (\@Species,\@GTFs,\@Genomes,$mergeorder_ref);

}






########################################################################
#
#  Function: SetMergeOrder
#
sub SetMergeOrder
{
    my $tree_str = shift;
    my $species_ref = shift;

    my @Species = @{$species_ref};
    my $num_species = scalar(@Species);

    # We'll benefit from being able to quickly check our species indices
    my %SpeciesToIndex;
    for (my $i=0; $i<$num_species; $i++) {
	$SpeciesToIndex{$Species[$i]} = $i+1; # Can't have '0' errors!
    }

    # If we don't have a species tree, we'll fake one.
    if (!$tree_str) {
	$tree_str = '('.$Species[0].','.$Species[1].')';
	for (my $i=2; $i<$num_species; $i++) {
	    $tree_str = '('.$tree_str.','.$Species[$i].')';
	}
    }

    # Before we dig in, let's confirm that all of our species are represented.
    # If any are missing, just toss 'em on.
    foreach my $species (@Species) {
	if ($tree_str !~ /\($species\,|\,$species\)/) {
	    $tree_str = '('.$tree_str.','.$species.')';
	}
    }

    # Let's also make sure that this looks like a valid tree
    $tree_str = ConfirmValidTree($tree_str);
    
    # Alrighty then, looks like it's time to build up a nice lil' merge order list!
    my @MergeOrder;
    my $num_merges = 0;
    while ($tree_str =~ /(\([\w\-\:]+\,[\w\-\:]+\))/) {
	
	my $merge_pair = $1;

	# Who's involved in this merge?
	$merge_pair =~ /\((\S+)\,(\S+)\)/;
	my $species1 = $1;
	my $species2 = $2;

	$merge_pair =~ s/\(/\\\(/;
	$merge_pair =~ s/\)/\\\)/;

	# Are these species names, or merge order coordinates?
	if ($species1 !~ /^M\d+/) {
	    if ($SpeciesToIndex{$species1}) {
		$species1 = $SpeciesToIndex{$species1}-1;
		$species1 = 'S'.$species1; 
	    } else {
		# UNKNOWN SPECIES!
		$species1 = 'S'.$num_species;
	    }
	}
	if ($species2 !~ /^M\d+/) {
	    if ($SpeciesToIndex{$species2}) {
		$species2 = $SpeciesToIndex{$species2}-1;
		$species2 = 'S'.$species2; 
	    } else {
		# UNKNOWN SPECIES!
		$species2 = 'S'.$num_species;
	    }
	}

	push(@MergeOrder,$species1.','.$species2);

	my $merge_id = 'M'.$num_merges;
	$num_merges++;
	
	# Clear that stinky ol' pair from the tree!
	$tree_str =~ s/$merge_pair/$merge_id/;
	
    }

    return(\@MergeOrder);
    
}




########################################################################
#
#  Function: ConfirmValidTree
#
sub ConfirmValidTree
{
    my $tree_str = shift;

    # First off, we'll do a raw count of parentheses, to make sure we have
    # equal numbers of open and close parentheses.  We'll also cover for
    # the user in case they accidentally left out any commas
    my $num_open = 0;
    my $num_close = 0;
    my $final_tree_str = '';
    my @TreeChrs = split(//,$tree_str);
    for (my $i=0; $i<scalar(@TreeChrs); $i++) {

	my $char = $TreeChrs[$i];
	$num_open++  if ($char eq '(');
	$num_close++ if ($char eq ')');

	if ($i && $char eq '(' && $TreeChrs[$i-1] eq ')') {
	    $final_tree_str = $final_tree_str.',';
	}
	
	if ($char =~ /[A-Za-z]/ && $TreeChrs[$i-1] eq ')') {
	    $final_tree_str = $final_tree_str.',';
	}
	
	$final_tree_str = $final_tree_str.$char;
	
	if ($char =~ /[A-Za-z]/ && $TreeChrs[$i+1] eq '(') {
	    $final_tree_str = $final_tree_str.',';
	}
	
    }

    if ($num_open != $num_close) {
	die "\n  ERROR:  Species tree has unequal numbers of open and close parens\n\n";
    }

    return $final_tree_str;
    
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
#     Where 'gene' can be a list of /-separated gene names (the FIRST will be
#     considered the "main" gene).  This allows the GTF to consider multiple
#     gene families' coordinates when mapping.
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

    my %SpeciesDir = %{$specseqdir_ref};

    # The first thing we'll do is scan through our database
    # making a set of species-AND-gene-family mini-databases.

    # We'll have a hash that counts the number of sequence characters,
    # both for species as wholes and specific genes within species
    my %CharCounts;

    # We'll also hash 'aliases' for genes, if the user has provided
    # a /-separated list in the 'gene' field
    my %GeneAliases;

    # We'll want to know whether we actually found any genomeless sequences
    my $misc_seqs = 0;

    my $ProteinDB = OpenInputFile($ProteinDB_name);
    my $seq_num = 0;
    my @OriginalSeqNames;
    my %GeneHash;
    my $spec_gene_filename;
    my $SpecGeneFile;
    my $species;
    my $gene;
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
		$gene = lc($2);

		# Any aliases? NOTE: We'll allow redundancy in this stage
		if ($gene =~ /\//) {
		    my @Aliases = split(/\//,$gene);
		    $gene = $Aliases[0];
		    for (my $i=0; $i<scalar(@Aliases); $i++) {
			if ($GeneAliases{$gene}) {
			    $GeneAliases{$gene} = $GeneAliases{$gene}.'/'.$Aliases[$i];
			} else {
			    $GeneAliases{$gene} = $Aliases[$i];
			}
		    }
		}

	    } elsif (ParseSeqNameAsUniProt($parsename)) {

		$parsename =~ /OS\=([^\=]+\=?)/;
		$species = lc($1);
		$species =~ s/\s+\S\S\=$//;

		# Replace illegal characters
		$species =~ s/\s|\||\&|\*|\$/\_/g;
		$species =~ s/\[/\(/g;
		$species =~ s/\]/\)/g;
		
		$parsename =~ /GN\=([^\=]+\=?)/;
		$gene = lc($1);
		$gene =~ s/\s+\S\S\=$//;
		
		# Replace illegal characters
		$gene =~ s/\s|\||\&|\*|\$/\_/g;
		$gene =~ s/\[/\(/g;
		$gene =~ s/\]/\)/g;

	    } else {
		die "\n  ERROR:  Unable to parse sequence name '$parsename'\n\n";
	    }

	    # Record the original name of this sequence, and then switch over to our
	    # nice, minimal name format for the duration.
	    push(@OriginalSeqNames,$seqname);
	    $seqname = $seq_num;

	    # If we already have a file open, close it
	    if ($spec_gene_filename) { close($SpecGeneFile); }

	    # Genomeless sequence ahoy!
	    if (!$SpeciesDir{$species}) {
		$species = 'Misc';
		$misc_seqs++;
	    }
	    
	    # Determine the specific file we want to write to
	    $spec_gene_filename = $SpeciesDir{$species}.'seqs/'.$gene.'.fa';
	    
	    # Print that beautiful nameline!
	    open($SpecGeneFile,'>>',$spec_gene_filename) || die "\n  ERROR:  Failed to open output species-gene database '$spec_gene_filename' ($species)\n\n";
	    print $SpecGeneFile ">$seq_num\n";

	    # Record that we've seen a gene
	    $GeneHash{$gene} = 1;

	    $seq_num++;

	} elsif ($spec_gene_filename) {

	    # If we the lengths of lines are large, we'll split them up
	    # (to make buffering easier for FastMap and friends)
	    my $max_chars_per_line = 240;
	    $line =~ s/\s//g; # Avoiding one pathological case
	    if (length($line) > $max_chars_per_line) {

		my $line_len = 0;
		foreach my $char (split(//,$line)) {
		    if ($line_len == $max_chars_per_line) {
			print $SpecGeneFile "\n";
			$line_len = 0;
		    }
		    print $SpecGeneFile "$char";
		    $line_len++;
		}
		print $SpecGeneFile "\n";
		
	    } else {
		print $SpecGeneFile "$line\n";
	    }

	    # Wait! Before you move on, how many characters are we lookin' at here?
	    my $linelen = length($line);
	    if ($CharCounts{$species}) { $CharCounts{$species} += $linelen; }
	    else                       { $CharCounts{$species}  = $linelen; }

	    if ($CharCounts{$species.'&'.$gene}) {
		$CharCounts{$species.'&'.$gene} = $CharCounts{$species.'&'.$gene} += $linelen;
	    } else {
		$CharCounts{$species.'&'.$gene} = $CharCounts{$species.'&'.$gene}  = $linelen;
	    }
	    
	}

    }
    if ($spec_gene_filename) { close($SpecGeneFile); }
    close($ProteinDB);

    # Next, we'll run through each of our species and determine which gene
    # families belong to which threads.
    foreach $species (keys %SpeciesDir) {

	# If this species didn't actually appear in the database, then we'll want
	# to remove it from circulation
	if (!$CharCounts{$species}) {
	    RunSystemCommand("rmdir \"$SpeciesDir{$species}\/alignments\"");
	    RunSystemCommand("rmdir \"$SpeciesDir{$species}\/mappings\"");
	    RunSystemCommand("rmdir \"$SpeciesDir{$species}\/seqs\"");
	    RunSystemCommand("rmdir \"$SpeciesDir{$species}\/timing\"") if ($gene_timing);
	    RunSystemCommand("rmdir \"$SpeciesDir{$species}\"");
	    $SpeciesDir{$species} = 0;
	    next;
	}

	# You gotta catch this ol' trickster before it fouls things up!
	next if ($species eq 'Misc');

	# How many sequence characters belong to this family?
	# How many characters (on average) should go to each thread?
	my $tot_chars = $CharCounts{$species};
	my $thread_portion = int($tot_chars/$num_cpus);

	# Great! Now we can run through and assign each thread its share of sequences
	my $ThreadGuide = OpenOutputFile($SpeciesDir{$species}.'seqs/Thread-Guide');
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
	    my $sp = $1; # 'species' is already in use :(
	    $gene  = $2;

	    next if ($sp ne $species);
	    
	    print $ThreadGuide "$threadnum $gene\n";

	    # If we've hit our portion, jump ship!
	    $threadchars += $CharCounts{$species_gene};
	    if ($threadchars >= $thread_portion && $threadnum < $num_cpus-1) {
		$threadnum++;
		$threadchars = 0;
	    }

	}
	close($ThreadGuide);

    }

    # If there were any aliases being used by proteins, be sure to record them!
    my @GenesWithAliases = sort keys %GeneAliases;
    if (scalar @GenesWithAliases) {

	# Invent a file to record alias information
	my $aliasf = OpenOutputFile($ResultsDir.'gene-aliases');
	foreach $gene (@GenesWithAliases) {

	    # We'll hash all of the aliases so we can avoid redundancy
	    my %UniqueAliases;
	    foreach my $alias (split(/\//,$GeneAliases{$gene})) {
		$UniqueAliases{$alias} = 1 unless ($alias eq $gene);
	    }

	    # Write 'em out!
	    print $aliasf "$gene";
	    foreach my $alias (sort keys %UniqueAliases) {
		print $aliasf " $alias";
	    }
	    print $aliasf "\n";
	    
	}
	close($aliasf);
	
    }

    # Do a quick switcheroo from a hash of seen genes to a list of genes
    my @AllGenes = sort keys %GeneHash;
    
    # This is it, baby!
    return (\@OriginalSeqNames,\@AllGenes,$misc_seqs,\%SpeciesDir);

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
    return 1 if ($seqname =~ /^[^\|]+\|[^\|]+\|[^\|]+\s*$/);
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
    # Make sure we can get everything we need to do our Mirage-y work with this
    # sequence.
    if ($seqname =~ /GN\=/ && $seqname =~ /OS\=/) {
	return 1;
    }
    return 0;
}





########################################################################
#
#  Function: SaveSpeciesMappings
#
sub SaveSpeciesMappings
{
    my $species_dirname = shift;
    my $origseqnames_ref = shift;

    my @OrigSeqNames = @{$origseqnames_ref};
    
    my $seqdirname = $species_dirname.'seqs/';
    my $mapdirname = $species_dirname.'mappings/';

    # First off, we'll make a list of files to check out, so we can get all
    # thready with it
    my $SeqDir = OpenDirectory($seqdirname);
    my @MappedGenes;
    while (my $fname = readdir($SeqDir)) {
	if ($fname =~ /^(\S+)\.quilter\.out/) {
	    push(@MappedGenes,$1);
	}
    }
    closedir($SeqDir);

    # Gettin' thready with it!
    my $num_threads = Min($num_cpus,scalar(@MappedGenes));
    return if ($num_threads == 0);

    my $threadID = SpawnProcesses($num_threads);
    my $start =  $threadID    * int(scalar(@MappedGenes)/$num_threads);
    my $end   = ($threadID+1) * int(scalar(@MappedGenes)/$num_threads);
    $end = scalar(@MappedGenes) if ($threadID == $num_threads-1);

    # Go through your assigned families and move them into their permanent residences
    for (my $i=$start; $i<$end; $i++) {
    
	my $gene = $MappedGenes[$i];
	my $infname = $seqdirname.$gene.'.quilter.out';
	my $inf = OpenInputFile($infname);
	my $outf = OpenOutputFile($mapdirname.$gene.'.out');
	
	while (my $line = <$inf>) {
	    $line =~ s/\n|\r//g;
	    if ($line =~ /Sequence ID\: (\d+)/) {
		my $seq_id = $1;
		my $seqname = $OrigSeqNames[$seq_id];
		print $outf "Sequence ID: $seqname\n";
	    } else {
		print $outf "$line\n";
	    }
	}
	
	close($inf);
	close($outf);
	
	RunSystemCommand("rm \"$infname\"");

	# If this family had any misses, move those over too
	my $missfname = $seqdirname.$gene.'.misses.out';
	if (-e $missfname) {

	    $inf = OpenInputFile($missfname);
	    $outf = OpenOutputFile($mapdirname.$gene.'.misses');

	    while (my $line = <$inf>) {
		$line =~ s/\n|\r//g;
		if ($line =~ /^(.*)\: (\d+)$/) {
		    my $details = $1;
		    my $seq_id = $2;
		    my $seqname = $OrigSeqNames[$seq_id];
		    print $outf "$details\: $seqname\n";
		}
	    }

	    close($inf);
	    close($outf);

	    RunSystemCommand("rm \"$missfname\"");
	    
	}
	
    }

    if ($threadID) { exit(0); }
    while (wait() != -1) {}
    
}





########################################################################
#
#  Function: AggregateMappingMisses
#
sub AggregateMappingMisses
{
    my $species_dirname  = shift;
    my $misses_dirname   = shift;
    my $origseqnames_ref = shift;

    my @OrigSeqNames = @{$origseqnames_ref};

    # What species is this?
    $species_dirname =~ /\/([^\/]+)\/$/;
    my $species = $1;

    # Make a hash of all the names of sequences that we missed.
    # This will include any sequences that hit to a 'noncanonical'
    # chromosome
    my %MappingMissesByGene;
    if (-e $species_dirname.'mapping-misses') {

	# Open up the list of mapping misses handed down by MapsToMSAs,
	# as well as an output file (where we'll write the final miss list)
	my $infname = $species_dirname.'mapping-misses';
	my $inf  = OpenInputFile($infname);
	my $outf = OpenOutputFile($misses_dirname.$species.'.misses');

	while (my $line = <$inf>) {
	    
	    $line =~ s/\n|\r//g;
	    next if (!$line);
	    
	    $line =~ /^(\S+) (\d+)\: (.*)$/;
	    my $gene = lc($1);
	    my $seq_id = $2;
	    my $reason = $3;

	    $seq_id++; # Need to handle possibility of sequence '0' :p
	    if ($MappingMissesByGene{$gene}) {
		$MappingMissesByGene{$gene} = $MappingMissesByGene{$gene}.','.$seq_id;
	    } else { 
		$MappingMissesByGene{$gene} = $seq_id;
	    }
	    $seq_id--;

	    my $seqname = $OrigSeqNames[$seq_id];
	    print $outf "$seqname -- $reason\n";
	    
	}

	close($outf);
	close($inf);

	# Since we now have the final record of misses in the designated folder,
	# let's clear out this file
	RunSystemCommand("rm \"$infname\"");
	
    }

    return \%MappingMissesByGene;
    
}





########################################################################
#
#  Function: AlignUnmappedSeqs
#
sub AlignUnmappedSeqs
{

    my $missesbygene_ref = shift;
    my $dirname = shift;

    my %MappingMissesByGene = %{$missesbygene_ref};

    # We'll need to do a quick reverse-engineer of the species name...
    $dirname =~ /\/([^\/]+)\/seqs\//;
    my $species = $1;

    # We'll make a list and determine how many threads to spawn
    my @GeneList = keys %MappingMissesByGene;
    my $num_threads = Min($num_cpus,scalar(@GeneList));
    return if ($num_threads == 0);

    # Let 'em know the truth
    DispProgMirage('msnw-init|0|'.$species);
    
    # Spawn them there threads!
    my $threadID = SpawnProcesses($num_threads);
    my $start =  $threadID    * int(scalar(@GeneList)/$num_threads);
    my $end   = ($threadID+1) * int(scalar(@GeneList)/$num_threads);
    $end = scalar(@GeneList) if ($threadID == $num_threads-1);

    # Go off!
    my $genes_completed = 0;
    for (my $gene_id=$start; $gene_id<$end; $gene_id++) {

	my $gene = $GeneList[$gene_id];

	my %SeqsToAlign;
	foreach my $seq_id (split(/\,/,$MappingMissesByGene{$gene})) {
	    # Recall that we incremented seq_id by 1 to be able to handle
	    # sequence '0' in a hash in 'AggregateMappingMisses'
	    $SeqsToAlign{$seq_id-1} = 1;
	}
	
	# We'll extract each of the unmapped sequences
	my @UnalignedSeqFiles;
	my $num_aligned = 0;
	my $num_unaligned = 0;
	my $unaligned = 0;
	my $outf;
	my $inf = OpenInputFile($dirname.$gene.'.fa');
	while (my $line = <$inf>) {

	    $line =~ s/\n|\r//g;
	    
	    if ($line =~ /\>(\d+)/) {

		# New seq! Who dis?
		my $seq_id = $1;

		# If we were writing to an output file, close it up
		close($outf) if ($unaligned);

		# Are we writing out this sequence, or has it been mapped?
		if ($SeqsToAlign{$seq_id}) {

		    $unaligned = 1;
		    $num_unaligned++;

		    my $outfname = $dirname.$seq_id.'.tmp.fa';
		    push(@UnalignedSeqFiles,$outfname);
		    $outf = OpenOutputFile($outfname);
		    print $outf ">$seq_id\n";
		    
		} else {
		    
		    $unaligned = 0;
		    $num_aligned++;
		    
		}

	    } elsif ($unaligned) {
		print $outf "$line\n";
	    }
	    
	}
	close($outf) if ($unaligned);
	close($inf);

	# If we have zero mapped sequences, surprise -- the first sequence
	# is our MSA! NOTE that we could more easily pick the last sequence,
	# but just to preserve order we do this sorta silly thing...
	my $msa_fname = $dirname.$gene.'.afa';
	my $i=0;
	if ($num_aligned == 0) {
	    RunSystemCommand("mv \"$UnalignedSeqFiles[0]\" \"$msa_fname\"");
	    $num_aligned++;
	    $i++;
	}

	# Where will we store the results of each alignment iteration?
	my $tmp_outfname = $dirname.$gene.'.tmp.afa';

	# Now we can just chug right on through, aligning seqs the ol'-fashioned way!
	while ($i<$num_unaligned) {
	    RunSystemCommand($multi_seq_nw." \"$UnalignedSeqFiles[$i]\" 1 \"$msa_fname\" $num_aligned -igBase 0 > \"$tmp_outfname\"");
	    RunSystemCommand("mv \"$tmp_outfname\" \"$msa_fname\"");
	    RunSystemCommand("rm \"$UnalignedSeqFiles[$i]\"");
	    $num_aligned++;
	    $i++;
	}

	# We don't want any fake '.afa' files clogging up the ol' pipey-pipes
	RunSystemCommand("rm \"$tmp_outfname\"") if (-e $tmp_outfname);

	# Nice work!
	$genes_completed++;
	DispProgMirage('msnw-loop|'.$threadID.'|'.$genes_completed);
	
    }

    # Make sure everybody's done before we head back home
    if ($threadID) { exit(0); }
    while (wait() != -1) {}

}





########################################################################
#
#  Function:  AlignMiscSeqs
#
sub AlignMiscSeqs
{
    my $dirname = shift;
    $dirname = $dirname.'seqs/';

    DispProgMirage('misc-ali');

    # We'll make a hash of all sequence names by gene, so we can take advantage
    # of the infrastructure of 'AlignUnmappedSeqs'
    my %SeqsByGene;
    my $Dir = OpenDirectory($dirname);
    while (my $fname = readdir($Dir)) {

	$fname = $dirname.$fname;
	if ($fname =~ /\/([^\/]+)\.fa/) {

	    my $gene = $1;

	    # We'll count the number of sequences, and if it's absurd (say, >30k)
	    # we don't generate an alignment (because it would take A BIT)
	    my $num_seqs = 0;
	    my $max_seqs = 30000;

	    my $grep = OpenSystemCommand("grep '>' \"$fname\"");
	    while (my $line = <$grep>) {
		if ($line =~ /\>(\d+)/) {

		    my $seq_id = $1+1; # Be careful with '0'!

		    # Too many unmapped sequences?
		    $num_seqs++;
		    if ($num_seqs > $max_seqs) {
			print "\n";
			print "  NOTE: Gene family '$gene' has a very large number of unmapped sequences (>30k).\n";
			print "        Final alignment will only incorporate mapped sequences for this family.\n";
			print "\n";
			$SeqsByGene{$gene} = 0;
			last;
		    }
		    
		    if ($SeqsByGene{$gene}) {
			$SeqsByGene{$gene} = $SeqsByGene{$gene}.','.$seq_id;
		    } else {
			$SeqsByGene{$gene} = $seq_id;
		    }

		}
	    }
	    close($grep);
	}
	
    }

    # Takin' advantage!
    AlignUnmappedSeqs(\%SeqsByGene,$dirname);

    ClearProgress();
    
}





########################################################################
#
#  Function: MergeAlignments
#
sub MergeAlignments
{
    my $species_ref = shift;
    my $speciesdir_ref = shift;
    my $mergeorder_ref = shift;
    my $allgenes_ref = shift;
    my $origseqnames_ref = shift;

    my @Species = @{$species_ref};
    my %SpeciesDir = %{$speciesdir_ref};
    my @MergeOrder = @{$mergeorder_ref};
    my @AllGenes = @{$allgenes_ref};
    my @OrigSeqNames = @{$origseqnames_ref};

    my $num_species = scalar(@Species);
    my $num_merges = scalar(@MergeOrder);

    # Announcement! We're the best!
    DispProgMirage('msnw-init|0|FINAL');

    # We'll need to add ARF information to the sequence names, so
    # let's get that data on-hand
    foreach my $species (@Species) {
	my $dirname = $SpeciesDir{$species};
	my $arfname = $dirname.'arfs';
	if (-e $arfname) {
	    my $inf = OpenInputFile($arfname);
	    while (my $line = <$inf>) {
		if ($line =~ /^(\d+) (\S+)/) {
		    my $seq_id = $1;
		    my $arf_data = $2;
		    if ($OrigSeqNames[$seq_id] =~ /\#/) {
			$OrigSeqNames[$seq_id] = $OrigSeqNames[$seq_id].' '.$arf_data;
		    } else {
			$OrigSeqNames[$seq_id] = $OrigSeqNames[$seq_id].' #'.$arf_data;
		    }
		}
	    }
	    close($inf);
	    RunSystemCommand("rm \"$arfname\"");
	}
    }

    # Where are we writing out our final MSAs?
    my $FinalDir = CreateDirectory($ResultsDir.'Final-MSAs');

    # Make a secret temporary directory to store intermediate results
    my $tmpdir = CreateDirectory($ResultsDir.'.tmp-msas');

    # Threadify and get 'em done!
    my $num_threads = Min($num_cpus,scalar(@AllGenes));
    my $threadID = SpawnProcesses($num_threads);

    my $start_gene_id =  $threadID    * int(scalar(@AllGenes)/$num_threads);
    my $end_gene_id   = ($threadID+1) * int(scalar(@AllGenes)/$num_threads);
    $end_gene_id = scalar(@AllGenes) if ($threadID == $num_threads-1);

    # Merge time!
    my $genes_completed = 0;
    for (my $gene_id=$start_gene_id; $gene_id<$end_gene_id; $gene_id++) {

	# For each species, see who has this gene in their repetoire
	my $gene = $AllGenes[$gene_id];
	my @SpeciesFiles;
	for (my $i=0; $i<$num_species; $i++) {
	    my $fname = $SpeciesDir{$Species[$i]}.'seqs/'.$gene.'.afa';
	    if (-e $fname) { $SpeciesFiles[$i] = $fname; }
	    else           { $SpeciesFiles[$i] = 0;      }
	}
	$SpeciesFiles[$num_species] = 0; # In case of unknown species

	# Merge 'em if you got 'em
	my @MergeFiles;
	for (my $merge=0; $merge<$num_merges; $merge++) {

	    $MergeOrder[$merge] =~ /^(\S+)\,(\S+)$/;
	    my $species1 = $1;
	    my $species2 = $2;

	    # Are we working with individual species, or pre-merged files?
	    if ($species1 =~ /S(\d+)/) {
		$species1 = $SpeciesFiles[$1];
		FinalizeIntraSpeciesMSA($species1,\@OrigSeqNames) if ($species1);
	    } elsif ($species1 =~ /M(\d+)/) {
		$species1 = $MergeFiles[$1];
	    }

	    if ($species2 =~ /S(\d+)/) {
		$species2 = $SpeciesFiles[$1];
		FinalizeIntraSpeciesMSA($species2,\@OrigSeqNames) if ($species2);
	    } elsif ($species2 =~ /M(\d+)/) {
		$species2 = $MergeFiles[$1];
	    }

	    # If these gene isn't attributable to either species, things are
	    # really easy!
	    if (!$species1 && !$species2) {
		$MergeFiles[$merge] = 0;
		next;
	    }

	    # Well, we at least know that there's a file we'll be creating
	    $MergeFiles[$merge] = $tmpdir.$threadID.'-'.$merge.'.afa';
	    
	    # If this gene is only present in one species, then we just move
	    # the file.
	    if (!$species1) {
		RunSystemCommand("mv \"$species2\" \"$MergeFiles[$merge]\"");
		next;
	    }
	    if (!$species2) {
		RunSystemCommand("mv \"$species1\" \"$MergeFiles[$merge]\"");
		next;
	    }

	    # UGH, we have to do WORK?!  UGGGGGHHHHHH

	    # How many sequences are in each file?
	    my $grep = OpenSystemCommand("grep '>' \"$species1\" | wc -l |");
	    my $numseqs1 = <$grep>;
	    $numseqs1 = int($numseqs1);
	    close($grep);

	    $grep = OpenSystemCommand("grep '>' \"$species2\" | wc -l |");
	    my $numseqs2 = <$grep>;
	    $numseqs2 = int($numseqs2);
	    close($grep);

	    # Gotta git down with that Needleman-Wunsch
	    my $nwcmd = $multi_seq_nw." \"$species1\" $numseqs1 ";
	    $nwcmd = $nwcmd."\"$species2\" $numseqs2 > \"$MergeFiles[$merge]\"";
	    RunSystemCommand($nwcmd);

	    # We'll remove the input files, just to avoid excessive storage use
	    RunSystemCommand("rm \"$species1\"");
	    RunSystemCommand("rm \"$species2\"");

	}

	# Files merged! Time to finalize the final sequence!
	my $infname = $MergeFiles[$num_merges-1];
	my $outfname = $FinalDir.$gene.'.afa';

	# (Sometimes finalization begins with a touch of cleanup)
	# (Otherwise, we'll need to confirm that splice site markers
	# extend through every row of the final MSA)
	if ($cleanMSA) {

	    my $tmpname = $outfname;
	    $tmpname =~ s/\.afa$/\.tmp/;

	    RunSystemCommand($final_msa." \"$infname\" \"$tmpname\"");
	    RunSystemCommand("rm \"$infname\"");

	    $infname = $tmpname;

	}

	# Big move's a-comin'!
	my $inf = OpenInputFile($infname);
	my $outf = OpenOutputFile($outfname);
	while (my $line = <$inf>) {
	    $line =~ s/\n|\r//g;
	    if ($line =~ /\>(\d+)/) {
		my $seq_id = $1;
		print $outf ">$OrigSeqNames[$seq_id]\n";
	    } else {
		print $outf "$line\n";
	    }
	}
	close($outf);
	close($inf);

	# Don't need that input file anymore
	RunSystemCommand("rm \"$infname\"");
	
	# Ummm, you okay, gorgeous? 'Cuz you're KILLING IT!
	$genes_completed++;
	DispProgMirage('msnw-loop|'.$threadID.'|'.$genes_completed);
	
    }
    
    # Amazing work, team!  Now bring it in to the chillout tent
    if ($threadID) { exit(0); }
    while (wait() != -1) {}

    # Clear out the secret temporary directory and bump on back to the top level!
    RunSystemCommand("rm -rf \"$tmpdir\" \&");
    
}





########################################################################
#
#  Function: FinalizeIntraSpeciesMSA
#
sub FinalizeIntraSpeciesMSA
{
    my $infname = shift;
    my $origseqnames_ref = shift;

    my @OrigSeqNames = @{$origseqnames_ref};

    $infname =~ /^(.*\/)seqs(\/[^\/]+)$/;
    my $outfname = $1.'alignments'.$2;

    my $inf  = OpenInputFile($infname);
    my $outf = OpenOutputFile($outfname);
    while (my $line = <$inf>) {
	$line =~ s/\n|\r//g;
	if ($line =~ /\>(\d+)/) {
	    my $seq_id = $1;
	    print $outf ">$OrigSeqNames[$seq_id]\n";
	} else {
	    print $outf "$line\n";
	}
    }
    close($outf);
    close($inf);
    
}






########################################################################
#
#  Function: ReorganizeResultsForMapping
#
sub ReorganizeResultsForMapping
{
    my $species = shift;
    my $species_data_dirname = shift;
    my $outdir_name = shift;
    my $seq_ids_to_names_ref = shift;

    my %SeqIDsToNames = %{$seq_ids_to_names_ref};

    my %GenesToMapFiles;
    my $SpeciesDataDir = OpenDirectory($species_data_dirname);
    while (my $fname = readdir($SpeciesDataDir)) {

	next if ($fname !~ /^(\S+)\.quilter\.out$/);
	my $gene = $1;

	$GenesToMapFiles{$gene} = $species_data_dirname.$fname;
	
    }
    my @SpeciesGeneList = keys %GenesToMapFiles;
    my $num_species_genes = scalar(@SpeciesGeneList);

    my $num_reorg_cpus = Min($num_cpus,$num_species_genes);
    my $thread_id = SpawnProcesses($num_reorg_cpus);

    my $AllThreadMissesFile = OpenOutputFile($species_data_dirname.$thread_id.'-misses');

    my $start_gene_id = $thread_id * int($num_species_genes / $num_reorg_cpus);
    my $end_gene_id = ($thread_id + 1) * int($num_species_genes / $num_reorg_cpus);
    $end_gene_id = $num_species_genes if ($thread_id == $num_reorg_cpus-1);

    for (my $gene_id = $start_gene_id; $gene_id < $end_gene_id; $gene_id++) {

	my $gene = $SpeciesGeneList[$gene_id];

	my $gene_dirname = $outdir_name.$gene.'/';
	if (!(-d $gene_dirname)) {
	    CreateDirectory($gene_dirname);
	}

	my $map_outfile_name = $gene_dirname.$species.'.mappings';
	my $MapOutFile = OpenOutputFile($map_outfile_name);
	my $num_maps = 0;

	my $miss_outfile_name = $gene_dirname.$species.'.misses';
	my $MissOutFile = OpenOutputFile($miss_outfile_name);
	my $num_misses = 0;
	
	my $MapInFile = OpenInputFile($GenesToMapFiles{$gene});
	while (my $line = <$MapInFile>) {

	    $line =~ s/\n|\r//g;

	    if ($line =~ /Sequence ID\: (\d+)/) {

		my $seq_name = $SeqIDsToNames{$1};

		$line = <$MapInFile>;
		$line =~ s/\n|\r//g;

		if ($line =~ /Unmapped/) {
		    print $MissOutFile "$seq_name\n";
		    print $AllThreadMissesFile "$seq_name\n";
		    $line = <$MapInFile>; # Eat the empty line
		    $num_misses++;
		} else {
		    print $MapOutFile "Sequence ID: $seq_name\n";
		    print $MapOutFile "$line\n";
		    $num_maps++;
		}
		
	    } else {

		print $MapOutFile "$line\n";

	    }
	    
	}
	close($MapInFile);

	close($MapOutFile);
	RunSystemCommand("rm -rf \"$map_outfile_name\"") if (!$num_maps);

	close($MissOutFile);
	RunSystemCommand("rm -rf \"$miss_outfile_name\"") if (!$num_misses);

    }

    close($AllThreadMissesFile);
    
    if ($thread_id) { exit(0); }
    while (wait() != -1) {}

    my $all_species_misses_fname = $misses_dirname.$species.'.misses';

    for (my $cpu_check_id = 0; $cpu_check_id < $num_cpus; $cpu_check_id++) {

	if (-e $species_data_dirname.$cpu_check_id.'-misses') {
	    my $cpu_misses_fname = $species_data_dirname.$cpu_check_id.'-misses';
	    RunSystemCommand("cat $cpu_misses_fname >> $all_species_misses_fname");
	}
	
	my $arf_filename = $species_data_dirname.$cpu_check_id.'-ARFs';
	next if (!(-e $arf_filename));

	my $ARFFile = OpenInputFile($arf_filename);
	while (my $line = <$ARFFile>) {

	    if ($line =~ /^(\d+)\s+ARFs?\:(\S+)/) {

		my $seq_name = $SeqIDsToNames{$1};
		my $arf_coords = $2;

		$seq_name =~ /^[^\|]+\|([^\|]+)\|/;
		my @SeqGeneList = split(/\//,$1);
		my $gene = $SeqGeneList[0];
		
		open(my $ARFOutFile,'>>',$outdir_name.$gene.'/'.$species.'.alt-reading-frames')
		    || die "\n  ERROR:  Failed to open file to write ARF data out for gene '$gene'\n\n";
		print $ARFOutFile "$seq_name\n";
		print $ARFOutFile "- Amino acids $arf_coords identified as ";
		if ($arf_coords =~ /\,/) {
		    print $ARFOutFile "alternative reading frames";
		} else {
		    print $ARFOutFile "an alternative reading frame";
		}
		print $ARFOutFile "\n\n";
		close($ARFOutFile);

	    }
	    
	}
	close($ARFFile);

    }

}






########################################################################
#
#  Function: EvaluateMissDir
#
sub EvaluateMissDir
{
    my $dirname = shift;

    my $dir = OpenDirectory($dirname);
    my $misses = 0;
    while (my $fname = readdir($dir)) {
	if ($fname =~ /\.misses/) {
	    $misses = 1;
	    last;
	}
    }
    closedir($dir);

    # Do we even need this directory?
    RunSystemCommand("rmdir \"$dirname\"") if ($misses == 0);
    
}





########################################################################
#
#  Function: ReportGranularQuilterTiming
#
sub ReportGranularQuilterTiming
{
    my $fname = shift;

    my $inf = OpenInputFile($fname);
    while (my $line = <$inf>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	$line =~ /^(.*\:) (\S+)$/;
	my $component = $1;
	my $seconds = $2;

	my $timing_str = SecondsToSMHDString($seconds);
	print "      - $component $timing_str\n";

    }
    close($inf);

}





# EOF

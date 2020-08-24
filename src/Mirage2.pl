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
sub ParseSeqNameAsMirage;
sub ParseSeqNameAsUniProt;
sub AggregateMappingMisses;
sub AlignUnmappedSeqs;
sub AlignMiscSeqs;


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


# We're going to need these friends
my $eslsfetch  = $location.'../inc/easel/miniapps/esl-sfetch';
my $eslseqstat = $location.'../inc/easel/miniapps/esl-seqstat';


# Where we'll be storing timing information
my $StartTime = [Time::HiRes::gettimeofday()];
my $IntervalStart;
my $IntervalEnd;
my $FinalTime;
my @QuilterTimeStats;
my @MapsToMSAsTimeStats;
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
my @Species = @{$speciesref};
my @GTFs    = @{$gtfsref};
my @Genomes = @{$genomesref};
my $num_species = scalar(@Species);

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
my %SpeciesMapDir;
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
my ($originalseqnames_ref,$genestospecies_ref)
    = GenerateSpeciesDBs($ProteinDB,$num_cpus,\%SpeciesSeqDir);
my @OriginalSeqNames = @{$originalseqnames_ref};
my %GenesToSpecies = %{$genestospecies_ref};

# I'm going to print off all the original sequence names to a file so
# if things go wrong during a run we can actually figure out who's who
my $seqnamef = OpenOutputFile($ResultsDir.'seq-name-guide');
for (my $i=0; $i<scalar(@OriginalSeqNames); $i++) {
    print $seqnamef "$i: $OriginalSeqNames[$i]\n";
}
close($seqnamef);


# Do our intra-species magic!
for (my $i=0; $i<$num_species; $i++) {

    my $species_dirname = $SpeciesDir{$Species[$i]};
    my $species_seqdir  = $SpeciesSeqDir{$Species[$i]};

    #
    #  Q U I L T E R
    #

    # Start a timer for Quilter
    $IntervalStart = [Time::HiRes::gettimeofday()];

    # Start assembling that command!
    my $QuilterCmd = 'perl '.$location.'Quilter2.pl';

    # KEY 1: Protein database (implicit in species directory, under 'seqs/')
    $QuilterCmd = $QuilterCmd.' '.$species_dirname;

    # KEY 2: Genome
    $QuilterCmd = $QuilterCmd.' '.$Genomes[$i];

    # KEY 3: GTF index (special case: - for "Just use BLAT")
    $QuilterCmd = $QuilterCmd.' '.$GTFs[$i];
	
    # Are we timing?
    $QuilterCmd = $QuilterCmd.' --time' if ($timed);
    
    # Make that daaaaaang call.
    RunSystemCommand($QuilterCmd);

    # Figure out how much time that took and record it
    $QuilterTimeStats[$i] = Time::HiRes::tv_interval($IntervalStart);

    #
    #  M A P s   T O   M S A s
    #

    # Now we want to check how long it takes to run MapsToMSAs
    $IntervalStart = [Time::HiRes::gettimeofday()];

    # As easy as it gets!
    my $MapsToMSAsCmd = 'perl '.$location.'MapsToMSAs.pl '.$species_seqdir;
    if (system($MapsToMSAsCmd)) {
	# HAMSTERS!
	die "\n  *  ERROR:  MapsToMSAs.pl failed during execution  *\n\n";
    }

    # Align any sequences that didn't fully map by getting all Needleman-Wunsch-y
    my $missesbygene_ref = AggregateMappingMisses($species_dirname);
    AlignUnmappedSeqs($missesbygene_ref,$species_seqdir);
    
    # Knock it off with that darn timing!
    $MapsToMSAsTimeStats[$i] = Time::HiRes::tv_interval($IntervalStart);

    # Species[i] over and out!
    ClearProgress();
    print "Intra-species alignment complete for $Species[$i]\n";

}


# HOLY COW! You just made a heckin' ton of MSAs!
# But don't forget about the 'Misc' sequences...
AlignMiscSeqs($SpeciesSeqDir{'Misc'});


# We're now going to track how long the whole final-MSA-generating
# part of the program takes
$IntervalStart = [Time::HiRes::gettimeofday()];


# Slap that stop-watch!
$TotalRuntime = Time::HiRes::tv_interval($StartTime);

# Print out runtime statistitcs
if ($verbose || $timed) {
    my $formattedTime;
    print "\n\n\n";
    print "  +---------------------------------------------------+\n";
    print "                   Runtime Statistics \n";
    print "  +---------------------------------------------------+\n";
    for (my $i=0; $i<$num_species; $i++) {
	
	print "\n";
	print "    Species        : $Species[$i]\n";
	
	$formattedTime = sprintf("%.3f",$QuilterTimeStats[$i]);
	print "    Hit Stitching  : $formattedTime seconds\n";
	
	$formattedTime = sprintf("%.3f",$MapsToMSAsTimeStats[$i]);
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
    push(@RequiredFiles,$location.'MapsToMSAs.pl');
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
    push(@RequiredFiles,$location.'MapsToMSAs.pl');
    push(@RequiredFiles,$location.'FinalMSA.pl');
    push(@RequiredFiles,$location.'FastMap2');
    push(@RequiredFiles,$location.'ExonWeaver');
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

    my %SpeciesSeqDir = %{$specseqdir_ref};

    # The first thing we'll do is scan through our database
    # making a set of species-AND-gene-family mini-databases.

    # We'll have a hash that counts the number of sequence characters,
    # both for species as wholes and specific genes within species
    my %CharCounts;

    # We'll also hash 'aliases' for genes, if the user has provided
    # a /-separated list in the 'gene' field
    my %GeneAliases;

    #ConfirmSSI($ProteinDB_name); # We do this in 'ParseArgs'
    my $ProteinDB = OpenInputFile($ProteinDB_name);
    my $seq_num = 0;
    my @OriginalSeqNames;
    my %SpeciesGeneHash;
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

		$parsename =~ /OS\=(\S+)/;
		$species = lc($1);

		$parsename =~ /GN\=(\S+)/;
		$gene = lc($1);

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
	    $spec_gene_filename = $SpeciesSeqDir{$species}.$gene.'.fa';

	    # Print that beautiful nameline!
	    open($SpecGeneFile,'>>',$spec_gene_filename) || die "\n  ERROR:  Failed to open output species-gene database '$spec_gene_filename' ($species)\n\n";
	    print $SpecGeneFile ">$seq_num\n";

	    # Record that this species has members of this family
	    $SpeciesGeneHash{$species.'|'.$gene} = 1;

	    $seq_num++;

	} elsif ($spec_gene_filename) {

	    print $SpecGeneFile "$line\n";

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

    # A final little piece of work we'll do is switch from our hash of gene/species
    # pairs to '1' to a hash from genes to a list of species that include them.
    my %GenesToSpecies;
    foreach my $gene_species_pair (keys %SpeciesGeneHash) {
	$SpeciesGeneHash{$gene_species_pair} =~ /^(\S+)\|(\S+)$/;
	$species = $1;
	$gene = $2;
	if ($GenesToSpecies{$gene}) {
	    $GenesToSpecies{$gene} = $GenesToSpecies{$gene}.','.$species;
	} else {
	    $GenesToSpecies{$gene} = $species;
	}
    }
    
    # This is it, baby!
    return (\@OriginalSeqNames,\%GenesToSpecies);

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
#  Function: AggregateMappingMisses
#
sub AggregateMappingMisses
{
    my $dirname = shift;

    # Make a hash of all the names of sequences that we missed.
    # This will include any sequences that hit to a 'noncanonical'
    # chromosome
    my %MappingMissesByGene;
    if (-e $dirname.'mapping-misses') {
	open(my $MissFile,'<',$dirname.'mapping-misses');
	while (my $line = <$MissFile>) {
	    $line =~ s/\n|\r//g;
	    next if (!$line);
	    $line =~ /^(\S+) (\d+)\:/;
	    my $gene = lc($1);
	    my $seq_id = $2;
	    if ($MappingMissesByGene{$gene}) {
		$MappingMissesByGene{$gene} = $MappingMissesByGene{$gene}.','.$seq_id;
	    } else { 
		$MappingMissesByGene{$gene} = $seq_id;
	    }
	}
	close($MissFile);
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

    # We'll make a list and determine how many threads to spawn
    my @GeneList = keys %MappingMissesByGene;
    my $num_threads = Min($num_cpus,scalar(@GeneList));
    return if ($num_threads == 0);

    # Spawn them there threads!
    my $threadID = SpawnProcesses($num_threads);
    my $start =  $threadID    * int(scalar(@GeneList)/$num_threads);
    my $end   = ($threadID+1) * int(scalar(@GeneList)/$num_threads);
    $end = scalar(@GeneList) if ($threadID == $num_threads-1);

    # Go off!
    for (my $gene_id=$start; $gene_id<$end; $gene_id++) {

	my $gene = $GeneList[$gene_id];

	my %SeqsToAlign;
	foreach my $seq_id (split(/\,/,$MappingMissesByGene{$gene})) {
	    $SeqsToAlign{$seq_id} = 1;
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
	    RunSystemCommand($location."MultiSeqNW \"$UnalignedSeqFiles[$i]\" 1 \"$msa_fname\" $num_aligned -igBase 0 > \"$tmp_outfname\"");
	    RunSystemCommand("mv \"$tmp_outfname\" \"$msa_fname\"");
	    $num_aligned++;
	    $i++;
	}

	# We don't want any fake '.afa' files clogging up the ol' pipey-pipes
	RunSystemCommand("rm \"$tmp_outfname\"") if (-e $tmp_outfname);
	
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

    # We'll make a hash of all sequence names by gene, so we can take advantage
    # of the infrastructure of 'AlignUnmappedSeqs'
    my %SeqsByGene;
    my $Dir = OpenDirectory($dirname);
    while (my $fname = readdir($Dir)) {

	$fname = $dirname.$fname;
	if ($fname =~ /\/([^\/]+)\.fa/) {
	    my $gene = $1;
	    my $grep = OpenSystemCommand("grep '>' \"$fname\"");
	    while (my $line = <$grep>) {
		if ($line =~ /\>(\d+)/) {
		    my $seq_id = $1;
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
    
}



# EOF

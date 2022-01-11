#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;
use Cwd;
use Getopt::Long;
use Time::HiRes;

# I AM UNAVOIDABLE
sub GetThisDir { my $lib = $0; $lib =~ s/\/Diviner.pl$//; return $lib; }
use lib GetThisDir();
use BureaucracyMirage;
use DisplayProgress;

# Subroutines
sub PrintUsage;
sub ParseArgs;
sub ParseGTF;
sub GetMappedSeqMSA;
sub ParseAFA;
sub RecordSplicedMSA;
sub ReduceMSAToSpecies;
sub FindGhostExons;
sub FindAliQualityDrops;
sub RecordGhostMSAs;
sub MatchMismatchScore;
sub LocalMatchMismatchAli;
sub GetB62Score;
sub MultiAminoSeqAli;
sub GetMapSummaryStats;
sub CollapseAndCountOverlaps;

# Added 'X' as ambiguity character (index 20)
my @Blosum62
    = ( 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,  0,
       -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3,  0,
       -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  0,
       -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  0,
	0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,  0,
       -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,
       -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  0,
	0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,  0,
       -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,
       -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,  0,
       -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,  0,
       -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,
       -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,  0,
       -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,  0,
       -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,  0,
	1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,
	0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,  0,
       -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,  0,
       -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,  0,
	0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0);
my %AminoIndex
    = ('A', 0,'C', 1,'D', 2,'E', 3,'F', 4,'G', 5,'H', 6,'I', 7,'K', 8,'L', 9,
       'M',10,'N',11,'P',12,'Q',13,'R',14,'S',15,'T',16,'V',17,'W',18,'Y',19,'X',20);
    

##############
#            #
#   SCRIPT   #
#            #
##############



if (@ARGV < 2) { PrintUsage(); }


# Figure out what the location of the Mirage src directory is
my $location = $0;
$location =~ s/Diviner\.pl$//;

# Find all the friends we're going to need inside Diviner
my $dependencies_ref = FindDependencies();
my %Dependencies = %{$dependencies_ref};

# We're going to need these friends
my $sindex  = $Dependencies{'sindex'};
my $sfetch  = $Dependencies{'sfetch'};
my $sstat   = $Dependencies{'sstat'};
my $tblastn = $Dependencies{'tblastn'};

# An astute observer will notice that these aren't the same settings as Quilter
# uses, which is because this isn't frickin' Quilter, geez.
# It's rad that our standardized filenames let us play like this!
$tblastn = $tblastn.' -outfmt 6 ';


# TODO: Make these options available as commandline arguments
my $options_ref = ParseArgs();
my %Options = %{$options_ref};
my $num_cpus = $Options{cpus};
my $outdirname = CreateDirectory($Options{outdirname});
my $outgenesdir = CreateDirectory($outdirname.'Results-by-Gene');
my $save_msas = $Options{savemsas}; # Do we want to write our spliced MSAs to files?
my $bad_ali_cutoff = $Options{alicutoff};

# Confirm that the input directory looks like the real deal
my $input_dirname = ConfirmDirectory($ARGV[0]);
my $final_results_dirname = ConfirmDirectory($input_dirname.'Final-MSAs');
my $all_species_dirname = ConfirmDirectory($input_dirname.'Species-MSAs');


# Start off the real work by parsing the species guide, which will give
# us genome locations and chromosome lengths.
my $tildedir_check = OpenSystemCommand('echo ~');
my $tildedir = <$tildedir_check>;
$tildedir =~ s/\n|\r//g;
$tildedir = ConfirmDirectory($tildedir);
close($tildedir_check);

my $SpeciesGuide = OpenInputFile($ARGV[1]);
my %SpeciesToGenomes;
my %SpeciesHasGTF;
my %GTFStartsToEnds;
my %GTFEndsToStarts;
my %ChrLensBySpecies;
while (my $line = <$SpeciesGuide>) {

    $line =~ s/\n|\r//g;
    next if ($line !~ /(\S+)\s+(\S+)\s+(\S+)/);
    my $species = lc($1);
    my $genome = $2;
    my $gtf = $3;
    $genome =~ s/\~\//$tildedir/;
    $gtf =~ s/\~\//$tildedir/;

    $SpeciesToGenomes{$species} = $genome;

    # I don't know why we'd be at this point, where we have
    # Mirage results without a genome index, but just in case...
    if (!(-e $genome.'.hsi')) {
	print "\n";
	print "  Warning: Couldn't find existing 'hsi' index for genome file '$genome'\n";
	print "           Index creation may take a couple of seconds\n\n";
	RunSystemCommand($sindex.' '.$genome);
    }

    # Get chromosome lengths for the genome
    my $SstatOut = OpenSystemCommand($sstat.' '.$genome);
    while (my $sstat_line = <$SstatOut>) {
	$sstat_line =~ s/\n|\r//g;
	if ($sstat_line =~ /^\: (\S+)\s+(\d+)/) {
	    my $chr = $1;
	    my $len = $2;
	    $ChrLensBySpecies{$species.'|'.$chr} = $len;
	}
    }
    close ($SstatOut);

    # We'll read in the species' GTF info so we can determine whether
    # any ghost exons we hit on are known coding regions (and not just
    # members of proteoforms that are missing from our database).
    if (-e $gtf) {
	$SpeciesHasGTF{$species} = 1;
	ParseGTF($species,$gtf);
    }
    
}
close($SpeciesGuide);


# As a bit of a sanity-check, if we don't have any genomes then we bail
if (scalar(keys %SpeciesToGenomes) == 0) {
    die "\n  ERROR:  Failed to locate usable genome mappings in species guide file '$ARGV[1]'\n\n";
}


# If we're writing out to a file, make a directory to store our spliced msas.
# If such a directory already exists, warn about overwriting, but trust that
# the most recent version of the software will give the best output.
my $spliced_dirname = $input_dirname.'Marked-Splice-Sites-MSAs/';
if ($save_msas) {
    if (-d $spliced_dirname) {
	print "\n";
	print "  Warning:  Existing directory of MSAs with marked splice sites located ($spliced_dirname)\n";
	print "            MSAs may be overwritten.\n\n";
    } else {
	CreateDirectory($spliced_dirname);
    }
    $spliced_dirname = ConfirmDirectory($spliced_dirname);
}


# We'll do a preliminary run through 'FinalMSAs' to compute the total number of
# genes, so that we can figure out how to allocate work for each of our processes.
my $FinalMSAs = OpenDirectory($final_results_dirname);
my @GeneList;
while (my $fname = readdir($FinalMSAs)) {
    next if ($fname !~ /(\S+)\.afa$/);
    push(@GeneList,$1);
}
closedir($FinalMSAs);

# Are we asking for too many cpus?
$num_cpus = Min($num_cpus,scalar(@GeneList));


# Look away, children, the processes are spawning!
my $threadID = SpawnProcesses($num_cpus);


# I've been borned! What's my job?
my $thread_portion = int(scalar(@GeneList)/$num_cpus);
my $start_gene_id = $threadID * $thread_portion;
my $end_gene_id = (1+$threadID) * $thread_portion;
$end_gene_id = scalar(@GeneList) if ($threadID == $num_cpus-1);

# Name temporary filenames that we'll want to use (and, while we're at it,
# fill in all of the wild 'n' wacky tblastn arguments we'll be using).
my $nucl_seq_fname = $outgenesdir.'nucl.tmp'.$threadID.'.fa';
my $prot_seq_fname = $outgenesdir.'prot.tmp'.$threadID.'.fa';
my $tbn_out_fname  = $outgenesdir.'tbn.tmp'.$threadID.'.out';
$tblastn = $tblastn.' -subject '.$nucl_seq_fname.' -query '.$prot_seq_fname;
$tblastn = $tblastn.' -out '.$tbn_out_fname.' 1>/dev/null 2>&1';


# TIME FOR THE MAIN EVENT!
#
# We'll do this the easiest (laziest) way, which will just be running
# through each of our final MSAs, loading in *only* those sequences
# which have mappings, and then reverse-engineering exon markers.
#
# Fun!
#
my $total_ghost_exons = 0;
my $total_ghosts_busted = 0;
my @GhostlyGenes;
for (my $gene_id=$start_gene_id; $gene_id<$end_gene_id; $gene_id++) {

    my $gene  = $GeneList[$gene_id];
    my $fname = $final_results_dirname.$gene.'.afa';

    # Pull in an MSA that's been (1) reduced to mapped sequences, and
    # (2) has splice site markers.
    my ($msa_ref,$mapmsa_ref,$seqnames_ref,$num_seqs,$msa_len,$speciestochrs_ref)
	= GetMappedSeqMSA($fname,$gene);

    # If we don't have anything to do with this family, skip to the next
    next if (!$msa_ref);
    
    my @MSA = @{$msa_ref};
    my @MapMSA = @{$mapmsa_ref};
    my @SeqNames = @{$seqnames_ref};
    my %SpeciesToChrs = %{$speciestochrs_ref};

    # At this point, it makes sense to create the result directory for this gene
    my $gene_outdir = CreateDirectory($outgenesdir.$gene);

    # Get your butt into this file, mister!
    RecordSplicedMSA(\@MSA,\@SeqNames,$num_seqs,$msa_len,$gene_outdir.$gene.'-seqs.afa');

    # Now we'll reduce our MSA even further, down to just one sequence per
    # species!
    ($msa_ref,$mapmsa_ref,$seqnames_ref,$num_seqs)
	= ReduceMSAToSpecies(\@MSA,\@MapMSA,\@SeqNames,$num_seqs,$msa_len);
    @MSA = @{$msa_ref};
    @MapMSA = @{$mapmsa_ref};
    my @SpeciesNames = @{$seqnames_ref};
    my $num_species = $num_seqs;

    # Write this one out, too!
    RecordSplicedMSA(\@MSA,\@SpeciesNames,$num_seqs,$msa_len,$gene_outdir.$gene.'-species.afa');
    
    # Now that we have our super-reduced splice-site-ified MSA, let's get real nasty
    # with it (by way of locating exons suggestive of "ghosts")!
    my ($num_ghost_exons, $num_ghosts_busted) =
	FindGhostExons($gene,\@MSA,\@MapMSA,\@SpeciesNames,$num_species,$msa_len,\%SpeciesToChrs);

    # Add onto our overarching tallies
    $total_ghost_exons += $num_ghost_exons;
    $total_ghosts_busted += $num_ghosts_busted;

    # If we didn't bust any ghosts, then there aren't any MSAs to build or
    # GTFs to examine... yay?
    next if (!$num_ghosts_busted);

    # Eeek! Ghosts!
    push(@GhostlyGenes,$gene);

    # Time to build some gorgeous translated MSAs!
    RecordGhostMSAs($gene);
    
}

# How'd I do?  I don't even know!
if ($threadID) {
    my $final_outf = OpenOutputFile($outgenesdir.$threadID.'.final-tally.out');
    print $final_outf "$total_ghosts_busted / $total_ghost_exons\n";
    for (my $i=0; $i<scalar(@GhostlyGenes); $i++) {
	print $final_outf '|' if ($i);
	print $final_outf "$GhostlyGenes[$i]";
    }
    print $final_outf "\n";
    close($final_outf);
}

# The age of threads is coming to a close!
if ($threadID) {
    exit(0);
}
while (wait() != -1) {}


# Woo-hoo!  Janitorial work is my favorite!
for ($threadID=0; $threadID<$num_cpus; $threadID++) {

    # Clear out all these files we don't need
    $nucl_seq_fname = $outgenesdir.'nucl.tmp'.$threadID.'.fa';
    $prot_seq_fname = $outgenesdir.'prot.tmp'.$threadID.'.fa';
    $tbn_out_fname  = $outgenesdir.'tbn.tmp'.$threadID.'.out';

    if (-e $nucl_seq_fname) { system("rm $nucl_seq_fname"); }
    if (-e $prot_seq_fname) { system("rm $prot_seq_fname"); }
    if (-e $tbn_out_fname)  { system("rm $tbn_out_fname");  }

    # How'd ya do, helper?
    if ($threadID) {
	
	my $final_infname = $outgenesdir.$threadID.'.final-tally.out';
	my $final_inf = OpenInputFile($final_infname);
    
	my $line = <$final_inf>;
	$line =~ /(\d+) \/ (\d+)/;
	$total_ghosts_busted += $1;
	$total_ghost_exons += $2;
	
	$line = <$final_inf>;
	$line =~ s/\n|\r//g;
	foreach my $gene (split(/\|/,$line)) {
	    push(@GhostlyGenes,$gene);
	}

	# Now erase every last trace of that darn helper from the Earth!
	close($final_inf);
	system("rm $final_infname");

    }
	
}


# This is what we call "the final judgment": did we find *anything*?
if ($total_ghost_exons == 0) {

    # Wowza yowza bo-bowza!  Really?  Nothing?  Okay...
    system("rm -rf $outdirname");
    print "\n  No potentially unannotated exons detected\n\n";

    exit(0);

}

# We'll write out summary statistics for the full collection of genes
# that had hits
GetMapSummaryStats(\@GhostlyGenes);

# WOOOOOOO, WE FOUND AT LEAST ONE THING TO POSSIBLY NOT CROSS-MAP!
my $bust_rate = int(1000.0*$total_ghosts_busted/$total_ghost_exons)/10.0;
print "\n";
print "  $total_ghost_exons absent exon homolog";
print "s" if ($total_ghost_exons != 1);
print " observed in proteoform set,\n";
if ($total_ghosts_busted == 1) {
    print "  1 ($bust_rate\%) of which has likely coding sequence in its genome\n";
} else {
    print "  $total_ghosts_busted ($bust_rate\%) of which have likely coding sequence in their genomes\n";
}
print "\n";
print "  Results in '$outdirname' (Summary info in file 'Search-Summary.out')\n";
print "\n";


1;






###################
#                 #
#   SUBROUTINES   #
#                 #
###################






###############################################################
#
#  Function: PrintUsage
#
sub PrintUsage
{
    print "\n";
    print "  USAGE :  ./Diviner.pl {OPT.s} [Mirage-Results] [Species-Guide]\n";
    print "\n";
    print "  OPT.s :  -cpus=[int]\n";
    print "           -outdirname=[string]\n";
    die "\n";
}





###############################################################
#
#  Function: ParseArgs
#
sub ParseArgs
{

    my %Options = (
	cpus => 1,
        outdirname => 'Diviner-Results',
        savemsas => 0,
	alicutoff => 0.4,
        );

    &GetOptions( 
        \%Options,
        "help",
	"cpus=i",
        "outdirname=s",
	"savemsas",
	"alicutoff=s"
        )
        || die "\n  ERROR:  Failed to parse command line arguments\n\n";

    if ($Options{help}) {
	die "\n  Help is on the way!\n\n"; # TODO
    }

    return \%Options;
    
}





###############################################################
#
#  Function: ParseGTF
#
sub ParseGTF
{
    my $species = shift;
    my $gtf_fname = shift;

    my $gtf = OpenInputFile($gtf_fname);
    while (my $line = <$gtf>) {

	my $entry_type = '';
        if ($line =~ /^\S+\s+\S+\s+(\S+)/) {
            $entry_type = lc($1);
        }
        
        if ($line =~ /^\#/ || ($entry_type ne 'exon' && $entry_type ne 'cds')) {
            next;
        }

        # Note that startPos is always less than endPos, even when we're
        # indexing into the reverse complement.
        $line =~ /^(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+([\+\-])/;
        my $chr     = $1;
        my $start   = $2;
        my $end     = $3;
        my $revcomp = $4;

	if ($revcomp eq '-') {
	    $chr = $chr.'[revcomp]';
	    my $temp = $start;
	    $start = $end;
	    $end = $temp;
	}

	my $key1 = $species.'|'.$chr.'|'.$start;
	my $key2 = $species.'|'.$chr.'|'.$end;

	# NOTE: This hash is declared top-level of script
	if ($GTFStartsToEnds{$key1}) {
	    $GTFStartsToEnds{$key1} = $GTFStartsToEnds{$key1}.'|'.$end;
	} else {
	    $GTFStartsToEnds{$key1} = $end;
	}

	if ($GTFEndsToStarts{$key2}) {
	    $GTFEndsToStarts{$key2} = $GTFEndsToStarts{$key2}.'|'.$start;
	} else {
	    $GTFEndsToStarts{$key2} = $start;
	}
	
    }
    close($gtf);

}





###############################################################
#
#  Function: GetMappedSeqMSA
#
sub GetMappedSeqMSA
{
    my $msa_infname = shift;
    my $gene = shift;

    # Start off by just reading in the MSA as-is
    my ($msa_ref,$seqnames_ref,$base_num_seqs,$base_msa_len) = ParseAFA($msa_infname);
    my @BaseMSA = @{$msa_ref};
    my @BaseSeqNames = @{$seqnames_ref};

    # Now we'll check each of the present species to see who in this
    # family is mapped.
    my %SpeciesToMapfiles;
    my %SpeciesToChr;
    my %MappedSeqs;
    foreach my $seqname (@BaseSeqNames) {

	# NOTE: We're currently assuming that our names are formatted in
	#       the Mirage-y way.
	# TODO: Change this to accept UniProt formatting
	$seqname =~ /^([^\|]+)\|/;
	my $species = lc($1);

	# If we've already pulled in this species' mapped sequences
	# (or confirmed that this whole species didn't map),
	# then skip right along, ya dinghus!
	next if ($SpeciesToMapfiles{$species});

	# Do we have a mapping file for this gene?
	my $mapfname = $all_species_dirname.$species.'/mappings/'.$gene.'.out';
	if (!(-e $mapfname)) {
	    $SpeciesToMapfiles{$species} = "-1";
	    next;
	}
	$SpeciesToMapfiles{$species} = $mapfname;

	# Oh boy! Time to read a file (right now, all we want is a roster
        # of mapped sequences)
	my $mapf = OpenInputFile($mapfname);
	my $canon_chr;
	while (my $line = <$mapf>) {

	    $line =~ s/\n|\r//g;

	    if ($line =~ /^Canonical Chromosome\: (\S+)/) {

		$canon_chr = $1;
		$SpeciesToChr{$species} = $canon_chr;
		
	    } elsif ($line =~ /^Sequence ID\: (\S+)/) {

		my $mapped_seqname = $1;

		# We'll need to make sure that we've hit on the right chromosome
		$line = <$mapf>; # map method
		$line = <$mapf>; # chromosome
		if ($line =~ /^Chromosome \: (\S+)/) {
		    my $mapped_chr = $1;
		    if ($mapped_chr eq $canon_chr) {
			$MappedSeqs{$mapped_seqname} = 1;
		    }
		}
		
	    }
	    
	}
	close($mapf);
	
    }

    
    # Before we go any further, if there's only one species mapping then we'll
    # bail.  Because 'SpeciesToMapfiles' initially has '-1' for any unmapped
    # species, we'll need to switch those to '0's
    my $num_mapped_species = 0;
    foreach my $species (keys %SpeciesToMapfiles) {
	if ($SpeciesToMapfiles{$species} eq "-1") {
	    $SpeciesToMapfiles{$species} = 0;
	} else {
	    $num_mapped_species++;
	}
    }
    return (0,0,0,0,0,0) if ($num_mapped_species <= 1);
    

    # Awesome!  Now we have a roster of sequences that mapped back to their
    # genomes, so let's take the initial step of reducing our MSA down to just
    # that collection.
    my @MSA;
    my @SeqNames;
    my $num_seqs = 0;
    for (my $i=0; $i<$base_num_seqs; $i++) {
	if ($MappedSeqs{$BaseSeqNames[$i]}) {

	    $SeqNames[$num_seqs] = $BaseSeqNames[$i];

	    for (my $j=0; $j<$base_msa_len; $j++) {
		$MSA[$num_seqs][$j] = $BaseMSA[$i][$j];
	    }
	    $num_seqs++;

	}
    }

    # We'll get rid of any all-gap columns to get our actual msa_len.
    # Note that this will clear any original '/' columns (if we're
    # looking at a dirty MSA), but that's fine since we're re-inserting them.
    my $msa_len = 0;
    for (my $j=0; $j<$base_msa_len; $j++) {
	my $residues = 0;
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($MSA[$i][$j] =~ /[A-Z]/) {
		$residues = 1;
		last;
	    }
	}
	if ($residues) {
	    for (my $i=0; $i<$num_seqs; $i++) {
		$MSA[$i][$msa_len] = $MSA[$i][$j];
	    }
	    $msa_len++;
	}
    }

    # Phew, that was some good MSA constructing!
    #
    # Next up will be associating each residue with its mapping coordinates
    # and determining where there ought to be splice site boundaries marked.
    #
    my @MapMSA;
    foreach my $species (keys %SpeciesToMapfiles) {

	next if (!$SpeciesToMapfiles{$species});

	my $mapf = OpenInputFile($SpeciesToMapfiles{$species});
	while (my $line = <$mapf>) {

	    # Fast-forward to the next canonically-mapped sequence
	    next if ($line !~ /Sequence ID\: (\S+)/);
	    my $seqname = $1;
	    next if (!$MappedSeqs{$seqname});

	    # We got a mappy boi!  Quick -- what row are you?
	    my $seq_row = 0;
	    while ($SeqNames[$seq_row] ne $seqname) {
		$seq_row++;
	    }

	    # How many exons does this sequence have?
	    $line = <$mapf>; # Map method
	    $line = <$mapf>; # Chromosome
	    $line = <$mapf>; # Num Exons!
	    $line =~ /Num Exons  \: (\d+)/;
	    my $num_exons = $1;

	    # We'll read in the full set of chromosome coordinates for
	    # this sequence before attributing them to MSA positions
	    # (just to simplify writing)
	    my $coord_list_str;
	    for (my $i=0; $i<$num_exons; $i++) {

		$line = <$mapf>; # Exon overview
		$line = <$mapf>; # Coords!

		$line =~ s/\n|\r//g;
		if ($coord_list_str) {
		    $coord_list_str = $coord_list_str.','.$line;
		} else {
		    $coord_list_str = $line;
		}
		
	    }

	    # Now we'll split up our coordinate list string into an
	    # array and walk along the MSA attributing each character
	    # to a coordinate
	    my @CoordList = split(/\,/,$coord_list_str);
	    my $next_coord = 0;
	    my $msa_pos = 0;
	    while ($msa_pos < $msa_len) {
		if ($MSA[$seq_row][$msa_pos] =~ /[A-Z]/) {
		    $MapMSA[$seq_row][$msa_pos] = $CoordList[$next_coord];
		    $next_coord++;
		} else {
		    $MapMSA[$seq_row][$msa_pos] = 0;
		}
		$msa_pos++;
	    }
	    

	}
	close($mapf);
	
    }

    
    # Well that sure gives us a 'MapMSA,' but what if we want to have a
    # splicey MSA?  Have no fear, sweet child, I'm gonna ruin you with
    # splice sites.

    
    # NOTE that the way we'll do this is by first identifying where there
    #   are positions in the MapMSA where the distance between two adjacent
    #   residues is >4, and recording all of those positions.
    my @LastCoord;
    for (my $i=0; $i<$num_seqs; $i++) {
	$LastCoord[$i] = 0;
    }
    
    # This will mark the first amino acid in each exon
    my %SpliceCols; 
    for (my $j=0; $j<$msa_len; $j++) {
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($MapMSA[$i][$j]) {
		if (abs($LastCoord[$i]-$MapMSA[$i][$j])>4) {
		    if ($SpliceCols{$j}) { $SpliceCols{$j}++; }
		    else                 { $SpliceCols{$j}=1; }
		}
		$LastCoord[$i] = $MapMSA[$i][$j];
	    }
	}
    }

    # Now we'll resolve any splice-sites that aren't universally agreed upon
    # by the sequences that have residues at the given positions.
    # We'll be (specifically) looking to discount columns where the
    # sequences who want the splice-site are immediately following a gap
    # (these look more like alt 3' splice sites), and then requiring that
    # the majority of non-gap-adjacent residues also call a given splice site.
    my @SpliceColList = sort {$a <=> $b} keys %SpliceCols;
    my @VetoedCols;
    my @VetoedColsVotes;
    for (my $splice_col_id=1; $splice_col_id<scalar(@SpliceColList); $splice_col_id++) { # '0' is mandatory

	my $splice_col = $SpliceColList[$splice_col_id];

	# How many sequences have a residue at *both* this column and the
	# last, and want it to be a splice site?
	my $paired_residue_for = 0;
	my $paired_residue_against = 0;
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($MSA[$i][$splice_col] =~ /[A-Z]/ && $MSA[$i][$splice_col-1] =~ /[A-Z]/) {

                # DEBUGGING
		if (!($MapMSA[$i][$splice_col] && $MapMSA[$i][$splice_col-1])) {
		    print "Mapping issue around col $splice_col in '$SeqNames[$i]'\n";
		    next;
		}
		
		if (abs($MapMSA[$i][$splice_col]-$MapMSA[$i][$splice_col-1])>4) {
		    $paired_residue_for++;
		} else {
		    $paired_residue_against++;
		}

	    }
	}

	# Has this splice site been vetoed?
	if ($paired_residue_against && $paired_residue_against >= $paired_residue_for) {
	    push(@VetoedCols,$splice_col);
	    push(@VetoedColsVotes,$SpliceCols{$splice_col});
	    $SpliceCols{$splice_col} = 0;
	    next;
	}
	
    }

    # If we had a run of vetoed columns that are close (1 or 2 residues apart) and
    # add up to majority support, then we'll re-introduce one of them into the list.
    my $run_start;
    my $run_end=1;
    while ($run_end < scalar(@VetoedCols)-1) {

	# Are we in a run of vetoed columns?
	$run_start = $run_end;
	$run_end++;
	my $run_votes = $VetoedColsVotes[$run_start];
	while ($run_end < scalar(@VetoedCols) && $VetoedCols[$run_end] - $VetoedCols[$run_end-1] < 3) {
	    $run_votes += $VetoedColsVotes[$run_end];
	    $run_end++;
	}

	# If there isn't a run of nearby vetoes, then we'll move right along
	next if ($run_end == $run_start+1);

	# Count the number of sequences that have non-gap characters
	# at one of the vetoed columns
	my $non_gappers = 0;
	for (my $i=0; $i<$num_seqs; $i++) {
	    for (my $col_id=$run_start; $col_id<$run_end; $col_id++) {
		my $msa_pos = $VetoedCols[$col_id];
		if ($MapMSA[$i][$msa_pos] || $MapMSA[$i][$msa_pos-1]) {
		    $non_gappers++;
		    last;
		}
	    }
	}

	# If enough columns have been vetoed to form a majority, then the most
	# used of the group will be designated a splice site.
	next if ($non_gappers >= $run_votes*2);

	my $top_veto_col = $run_start;
	my $top_veto_vote = $VetoedColsVotes[$run_start];
	for (my $i=$run_start+1; $i<$run_end; $i++) {
	    if ($top_veto_vote < $VetoedColsVotes[$i]) {
		$top_veto_col = $i;
		$top_veto_vote = $VetoedColsVotes[$i];
	    }
	}

	# THE PEOPLE HAVE SPOKEN!
	$SpliceCols{$top_veto_col} = $run_votes;
	
    }

    # Get the final list of splice site columns
    @SpliceColList = sort {$a <=> $b} keys %SpliceCols;

    # Now we construct the final splice-site-ified MSAs!
    my @SplicedMSA;
    my @SplicedMapMSA;
    my $spliced_msa_len = 0;
    my $splice_col_id = 0;
    for (my $j=0; $j<$msa_len; $j++) {

	# Have we hit our next splice column?
	if ($splice_col_id < scalar(@SpliceColList) && $SpliceColList[$splice_col_id] == $j) {
	    for (my $i=0; $i<$num_seqs; $i++) {
		$SplicedMSA[$i][$spliced_msa_len] = '/';
		$SplicedMapMSA[$i][$spliced_msa_len] = 0;
	    }
	    $spliced_msa_len++;
	    $splice_col_id++;
	}

	for (my $i=0; $i<$num_seqs; $i++) {
	    $SplicedMSA[$i][$spliced_msa_len] = $MSA[$i][$j];
	    $SplicedMapMSA[$i][$spliced_msa_len] = $MapMSA[$i][$j];
	}
	$spliced_msa_len++;

    }

    # Finally, round things out with a splice site column
    for (my $i=0; $i<$num_seqs; $i++) {
	$SplicedMSA[$i][$spliced_msa_len] = '/';
	$SplicedMapMSA[$i][$spliced_msa_len] = 0;
    }
    $spliced_msa_len++;

    # Pass everything we got back on up the chain of command!
    return (\@SplicedMSA,\@SplicedMapMSA,\@SeqNames,$num_seqs,$spliced_msa_len,\%SpeciesToChr);
    
}





###############################################################
#
#  Function: ParseAFA
#
sub ParseAFA
{
    my $fname = shift;

    my $inf = OpenInputFile($fname);
    my $num_seqs = -1;
    my $msa_len;
    my @MSA;
    my @SeqNames;
    while (my $line = <$inf>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	if ($line =~ /\>(\S+)/) {

	    my $seqname = $1;
	    push(@SeqNames,$seqname);
	    $num_seqs++;
	    $msa_len = 0;

	} else {

	    foreach my $char (split(//,$line)) {
		$MSA[$num_seqs][$msa_len] = uc($char);
		$msa_len++;
	    }

	}
	
    }
    close($inf);

    $num_seqs++;

    return(\@MSA,\@SeqNames,$num_seqs,$msa_len);

}





###############################################################
#
#  Function: RecordSplicedMSA
#
sub RecordSplicedMSA
{

    my $msa_ref = shift;
    my $seqnames_ref = shift;
    my $num_seqs = shift;
    my $msa_len = shift;
    my $outfname = shift;

    my @MSA = @{$msa_ref};
    my @SeqNames = @{$seqnames_ref};
    
    my $outf = OpenOutputFile($outfname);
    for (my $i=0; $i<$num_seqs; $i++) {

	print $outf ">$SeqNames[$i]\n";
	for (my $j=0; $j<$msa_len; $j++) {
	    print $outf "$MSA[$i][$j]";
	    print $outf "\n" if (($j+1) % 60 == 0);
	}
	print $outf "\n" if ($msa_len % 60);
	#print $outf "\n";
	
    }
    close($outf);

}





###############################################################
#
#  Function: ReduceMSAToSpecies
#
sub ReduceMSAToSpecies
{
    my $msa_ref = shift;
    my $mapmsa_ref = shift;
    my $seqnames_ref = shift;
    my $num_seqs = shift;
    my $msa_len = shift;

    my @MSA = @{$msa_ref};
    my @MapMSA = @{$mapmsa_ref};
    my @SeqNames = @{$seqnames_ref};

    # First off, who are our species?  Which sequences belong to which species?
    my %SpeciesToSeqIDs;
    for (my $i=0; $i<$num_seqs; $i++) {
	$SeqNames[$i] =~ /^([^\|]+)\|/;
	my $species = $1;
	if ($SpeciesToSeqIDs{$species}) {
	    $SpeciesToSeqIDs{$species} = $SpeciesToSeqIDs{$species}.','.($i+1);
	} else {
	    $SpeciesToSeqIDs{$species} = ($i+1);
	}
    }

    # Alrighty! Let's get reducing!
    my @SpeciesNames;
    my @SpeciesMSA;
    my @SpeciesMapMSA;
    my $num_species = 0;
    foreach my $species (keys %SpeciesToSeqIDs) {

	push(@SpeciesNames,$species);

	# Get the sequence IDs for every member of this species.
	# Note that we'll need to decrement by one because we previously
	# had to increment by one (to not throw out sequence '0')
	my @SeqIDs = split(/\,/,$SpeciesToSeqIDs{$species});
	for (my $i=0; $i<scalar(@SeqIDs); $i++) {
	    $SeqIDs[$i] = $SeqIDs[$i]-1;
	}

	# Great!  Now it's time to get reducing!  To do this, we'll just
	# walk along the MSA, and wherever there's a residue-containing
	# column we'll take a poll of what the residue should be, and majority
	# wins
	for (my $j=0; $j<$msa_len; $j++) {

	    # If this is a splice-site column, we can take it super easy
	    if ($MSA[0][$j] eq '/') {
		$SpeciesMSA[$num_species][$j] = '/';
		$SpeciesMapMSA[$num_species][$j] = 0;
		next;
	    }
	    
	    my %Residues;
	    foreach my $seq_id (@SeqIDs) {
		if ($MSA[$seq_id][$j] =~ /[A-Z]/) {
		    if ($Residues{$MSA[$seq_id][$j]}) {
			$Residues{$MSA[$seq_id][$j]}++;
		    } else {
			$Residues{$MSA[$seq_id][$j]}=1;
		    }
		}
	    }

	    # If this is a gap column, we can take it sorta easy
	    if (scalar(keys %Residues) == 0) {
		$SpeciesMSA[$num_species][$j] = '-';
		$SpeciesMapMSA[$num_species][$j] = 0;
		next;
	    }

	    # Alrighty then, who's the lucky residue?
	    my @ResidueList = keys %Residues;
	    my $top_residue = $ResidueList[0];
	    my $top_residue_count = $Residues{$top_residue};
	    for (my $i=1; $i<scalar(@ResidueList); $i++) {
		if ($Residues{$ResidueList[$i]} > $top_residue_count) {
		    $top_residue = $ResidueList[$i];
		    $top_residue_count = $Residues{$top_residue};
		}
	    }

	    # Great stuff!  Now, we just need to get a mapping coordinate from
	    # someone who used the top residue.
	    my $map_coord;
	    foreach my $seq_id (@SeqIDs) {
		if ($MSA[$seq_id][$j] eq $top_residue) {
		    $map_coord = $MapMSA[$seq_id][$j];
		    last;
		}
	    }

	    # Record that residue and mapping coordinate!
	    $SpeciesMSA[$num_species][$j] = $top_residue;
	    $SpeciesMapMSA[$num_species][$j] = $map_coord;
	    
	}

	# Another species in the bag!
	$num_species++;

    }

    # That was too easy!
    return(\@SpeciesMSA,\@SpeciesMapMSA,\@SpeciesNames,$num_species);
    
}





###############################################################
#
#  Function: FindGhostExons
#
sub FindGhostExons
{
    my $gene = shift;
    my $msa_ref = shift;
    my $mapmsa_ref = shift;
    my $speciesnames_ref = shift;
    my $num_species = shift;
    my $msa_len = shift;
    my $speciestochrs_ref = shift;

    my @MSA = @{$msa_ref};
    my @MapMSA = @{$mapmsa_ref};
    my @SpeciesNames = @{$speciesnames_ref};
    my %SpeciesToChrs = %{$speciestochrs_ref};

    # First off, let's figure out where the starts (inclusive) and ends (exclusive)
    # of our exons are, as well as how many exons we have
    my @ExonStarts;
    my @ExonEnds;
    my $num_exons = 0;
    for (my $j=0; $j<$msa_len; $j++) {
	if ($MSA[0][$j] eq '/') {
	    $num_exons++;
	    push(@ExonEnds,$j) if ($j);
	    push(@ExonStarts,$j+1) if ($j+1 < $msa_len);
	}
    }
    $num_exons--; # We'll have overcounted by one, but that's okie-dokie

    # Before we get into the nastiness, we'll set some variables to represent the
    # minimum number of amino acids for 'using an exon' and the maximum number
    # of amino acids for 'not using an exon' -- both Inclusive
    #
    # This is intended to pre-empt any occurrences of micro-exons or other
    # weird nonstandard splicing events that might cause an exon to technically
    # have amino acids in it, but clearly only as an artifact of computation.
    #
    # Similarly, we have a minimum number of aminos for a sequence to be worth
    # tblastn-ing.
    my $min_use_aminos = 6;
    my $max_nonuse_aminos = 3;
    my $min_search_aminos = 10;

    # Alright, now for the pain!  We'll do a pairwise comparison of species,
    # seeing which species suggest ghost exons in which other ones.
    # Wherever we find such pairs, we'll record (1) the sequence we're looking for,
    # (2) the species we're looking for it in, (3) the range we're looking for it in,
    # and (4) the species we found it in.
    my @SearchSeqs;
    my @TargetSpecies;
    my @TargetSpeciesRange;
    my @SourceSpecies;
    my @MSAExonRanges;
    my @UsedRegions;
    my $num_ghost_exons = 0;
    for (my $s1=0; $s1<$num_species-1; $s1++) {

	my $species1 = $SpeciesNames[$s1];

	for (my $s2=$s1+1; $s2<$num_species; $s2++) {

	    my $species2 = $SpeciesNames[$s2];

	    # We'll kick things off by computing the number of matches / mismatches
	    # between the two sequences, as well as the total number of amino acids
	    # each sequence has, for each exon.  We'll use these bits of information
	    # to make a determination as to the nastiness of the pairing.
	    my @ExonMatches;
	    my @ExonMismatches;
	    my @ExonPctsID;
	    my @S1ExonAminos;
	    my @S2ExonAminos;
	    my @IsNastyExonPair;
	    my @IsUsedAtAll;
	    my @MicroExonic;
	    for (my $exon_id=0; $exon_id<$num_exons; $exon_id++) {

		$ExonMatches[$exon_id] = 0;
		$ExonMismatches[$exon_id] = 0;
		$S1ExonAminos[$exon_id] = 0;
		$S2ExonAminos[$exon_id] = 0;

		# Let's start with the straight facts
		for (my $j=$ExonStarts[$exon_id]; $j<$ExonEnds[$exon_id]; $j++) {

		    if ($MSA[$s1][$j] =~ /[A-Z]/) {
			$S1ExonAminos[$exon_id]++;
			if ($MSA[$s2][$j] =~ /[A-Z]/) {
			    $S2ExonAminos[$exon_id]++;
			    if ($MSA[$s1][$j] eq $MSA[$s2][$j]) {
				$ExonMatches[$exon_id]++;
			    } else {
				$ExonMismatches[$exon_id]++;
			    }
			}
		    } elsif ($MSA[$s2][$j] =~ /[A-Z]/) {
			$S2ExonAminos[$exon_id]++;
		    }
		    
		}


		# But now it's time for a hot take!
		$ExonPctsID[$exon_id] = 0;
		$IsNastyExonPair[$exon_id] = 0;
		$MicroExonic[$exon_id] = 0;
		$IsUsedAtAll[$exon_id] = 1;

		# Catch micro exons
		if ($S1ExonAminos[$exon_id] <= $max_nonuse_aminos
		    || $S2ExonAminos[$exon_id] <= $max_nonuse_aminos) {
		    $MicroExonic[$exon_id] = 1;
		}

		if ($ExonMismatches[$exon_id] || $ExonMatches[$exon_id]) {
		    $ExonPctsID[$exon_id] = 100.0 * $ExonMatches[$exon_id] / ($ExonMatches[$exon_id]+$ExonMismatches[$exon_id]+0.0);
		} else {
		    $IsNastyExonPair[$exon_id] = 0;
		    if ($S1ExonAminos[$exon_id] == 0 && $S2ExonAminos[$exon_id] == 0) {
			$IsUsedAtAll[$exon_id] = 0;
		    }
		    next;
		}

		# If we have a micro exon we'll go ahead an bail (but we do
		# hang onto the pct id for some reason...)
		next if ($MicroExonic[$exon_id]);

		# Way to nastiness 1: One sequence doesn't use this exon, but
		#   the other does.
		#
		if ($S1ExonAminos[$exon_id] >= $min_use_aminos
		    && $S2ExonAminos[$exon_id] <= $max_nonuse_aminos) {
		    $IsNastyExonPair[$exon_id] = 1;
		    next;
		}
		if ($S2ExonAminos[$exon_id] >= $min_use_aminos
		    && $S1ExonAminos[$exon_id] <= $max_nonuse_aminos) {
		    $IsNastyExonPair[$exon_id] = 1;
		    next;
		}

		# Way to nastiness 2: Both sequences use this exon, but the
		#   total number of aligned columns (matches or mismatches)
		#   is less than half of the length of the shorter of the two
		if (2*($ExonMatches[$exon_id]+$ExonMismatches[$exon_id])
		    < Min($S1ExonAminos[$exon_id],$S2ExonAminos[$exon_id])) {
		    $IsNastyExonPair[$exon_id] = 2;
		    next;
		}

		# Way to nastiness 3: The sequences have a reasonable alignment
		#   region, but the alignment itself is ultra-low-quality.
		my $pct_id = (0.0 + $ExonMatches[$exon_id]) / (0.0 + $ExonMatches[$exon_id] + $ExonMismatches[$exon_id]);
		if ($pct_id <= $bad_ali_cutoff) {
		    $IsNastyExonPair[$exon_id] = 3;
		    next;
		}

	    }

	    # The last thing we'll do to identify nasty exon pairs is check
	    # for drops in alignment quality.  I'm going to bump this over to
	    # a seperate function.
	    my $quality_drops_ref = FindAliQualityDrops(\@ExonPctsID,$num_exons);
	    if ($quality_drops_ref) {
		foreach my $exon_id (@{$quality_drops_ref}) {
		    if ($IsNastyExonPair[$exon_id] == 0) {
			$IsNastyExonPair[$exon_id] = 4;
		    }
		}
	    }

	    # Just kidding about the last time I said 'this is the last thing
	    # we'll do...' Now it's time to let our NastyExons (tm) seep into
	    # any adjacent micro exons.
	    # We'll do this in a forward and a backward pass
	    my $in_nasty_zone = $IsNastyExonPair[0];
	    for (my $exon_id=0; $exon_id<$num_exons; $exon_id++) {
		next if (!$IsUsedAtAll[$exon_id]);
		if ($MicroExonic[$exon_id] && $in_nasty_zone) {
		    $IsNastyExonPair[$exon_id] = 5;
		} else {
		    $in_nasty_zone = $IsNastyExonPair[$exon_id];
		}
	    }

	    for (my $exon_id=$num_exons-1; $exon_id>=0; $exon_id--) {
		next if (!$IsUsedAtAll[$exon_id]);
		if ($MicroExonic[$exon_id] && $in_nasty_zone) {
		    $IsNastyExonPair[$exon_id] = 5;
		} else {
		    $in_nasty_zone = $IsNastyExonPair[$exon_id];
		}
	    }

	    # Radical!  Now we can condense each of our runs of contiguous
	    # nasty exon pairs, extract the relevant sequence, and BOOM!
	    #
	    # NOTE: We don't really need the Matches/Mismatches/S_ExonAminos,
	    #   but it might be good to have at some later point for either
	    #   debugging or more sprawling output, so let's not worry about
	    #   the (probably) unnecessary data hoarding.
	    my @NastyRunStarts;
	    my @NastyRunEnds;
	    my $num_nasty_runs = 0;
	    for (my $exon_id=0; $exon_id<$num_exons; $exon_id++) {
		if ($IsNastyExonPair[$exon_id]) {
		    if ($num_nasty_runs == 0 || $NastyRunEnds[$num_nasty_runs-1] != $exon_id-1) {
			$NastyRunStarts[$num_nasty_runs] = $exon_id;
			$num_nasty_runs++;
		    }
		    $NastyRunEnds[$num_nasty_runs-1] = $exon_id;
		}
	    }

	    # Moment of truth -- are we... NASTY?!
	    next if ($num_nasty_runs == 0);

	    # NASTY! NASTY! NASTY!

	    # Alrighty then, time to take each of our runs and pull out the data
	    # that we'll need to perform our tblastn searches.
	    for (my $run_id=0; $run_id<$num_nasty_runs; $run_id++) {

		my $start_exon = $NastyRunStarts[$run_id];
		my $end_exon = $NastyRunEnds[$run_id];

		# The sequence from seq (species) 1 & the sequence from seq 2
		my $search_seq_1 = '';
		my $search_seq_2 = '';
		for (my $j=$ExonStarts[$start_exon]; $j<$ExonEnds[$end_exon]; $j++) {
		    if ($MSA[$s1][$j] =~ /[A-Z]/) {
			$search_seq_1 = $search_seq_1.$MSA[$s1][$j];
		    }
		    if ($MSA[$s2][$j] =~ /[A-Z]/) {
			$search_seq_2 = $search_seq_2.$MSA[$s2][$j];
		    }
		}

		# We'll be searching for this sequence from species 1 against
		# the genome of species 2.
		if (length($search_seq_1) >= $min_search_aminos) {

		    # [1] Find the nearest genomic coordinates outside of the nasty
		    #     region for sequence 2, starting with the left

		    my $left_bound;
		    if ($start_exon) {
			my $j=$ExonEnds[$start_exon-1]-1;
			while ($j) {
			    if ($MapMSA[$s2][$j]) {
				$left_bound = $MapMSA[$s2][$j];
				last;
			    }
			    $j--;
			}
		    }
		    if (!$left_bound) {
			my $j=0;
			while (!$MapMSA[$s2][$j]) {
			    $j++;
			}
			$left_bound = '[start-of-coding-region:'.$MapMSA[$s2][$j].']';
		    }

		    # [2] Find the coordinates to the right using a method that mirrors
		    #     how we got coordinates to the left

		    my $right_bound;
		    if ($end_exon+1 < $num_exons) {
			my $j=$ExonStarts[$end_exon+1];
			while ($j<$msa_len) {
			    if ($MapMSA[$s2][$j]) {
				$right_bound = $MapMSA[$s2][$j];
				last;
			    }
			    $j++;
			}
		    }
		    if (!$right_bound) {
			# Find the first coordinate for this species, noting it as such
			my $j=$msa_len-1;
			while (!$MapMSA[$s2][$j]) {
			    $j--;
			}
			$right_bound = '[end-of-coding-region:'.$MapMSA[$s2][$j].']';
		    }

		    
		    # [3] We'll need to disqualify any nucleotides currently identified
		    #     as contributing to our mapping, so we'll go exon-by-exon and
		    #     record any already-used regions.

		    my $used_regions_str;
		    for (my $exon_id=$start_exon; $exon_id<=$end_exon; $exon_id++) {

			my $region_start;
			for (my $j=$ExonStarts[$exon_id]; $j<$ExonEnds[$exon_id]; $j++) {
			    if ($MapMSA[$s2][$j]) {
				$region_start = $MapMSA[$s2][$j];
				last;
			    }
			}

			my $region_end;
			for (my $j=$ExonEnds[$exon_id]; $j>$ExonStarts[$exon_id]; $j--) {
			    if ($MapMSA[$s2][$j]) {
				$region_end = $MapMSA[$s2][$j];
				last;
			    }
			}

			if ($region_start && $region_end) {
			    if ($used_regions_str) {
				$used_regions_str = $used_regions_str.'&'.$region_start.'|'.$region_end;
			    } else {
				$used_regions_str = $region_start.'|'.$region_end;
			    }
			}
			
		    }

		    # Just to make sure we're sticking something into our list
		    if (!$used_regions_str) {
			$used_regions_str = '0|0';
		    }


		    # Heck yeah! Let's scream (and shout)!
		    push(@SearchSeqs,lc($search_seq_1));
		    push(@TargetSpecies,$species2);
		    push(@TargetSpeciesRange,$left_bound.'|'.$right_bound);
		    push(@SourceSpecies,$species1);
		    push(@MSAExonRanges,($start_exon+1).'..'.($end_exon+1));
		    push(@UsedRegions,$used_regions_str);
		    $num_ghost_exons++;
		    
		}


		# We'll be searching for this sequence from species 2 against
		# the genome of species 1.
		if (length($search_seq_2) >= $min_search_aminos) {

		    # [1] Find the nearest genomic coordinates outside of the nasty
		    #     region for sequence 1, starting with the left

		    my $left_bound;
		    my $left_bound_msa_pos;
		    if ($start_exon) {
			my $j=$ExonEnds[$start_exon-1]-1;
			while ($j) {
			    if ($MapMSA[$s1][$j]) {
				$left_bound = $MapMSA[$s1][$j];
				$left_bound_msa_pos = $j;
				last;
			    }
			    $j--;
			}
		    }
		    if (!$left_bound) {
			my $j=0;
			while (!$MapMSA[$s1][$j]) {
			    $j++;
			}
			$left_bound = '[start-of-coding-region:'.$MapMSA[$s1][$j].']';
			$left_bound_msa_pos = $j;
		    }

		    # [2] Find the coordinates to the right using a method that mirrors
		    #     how we got coordinates to the left

		    my $right_bound;
		    my $right_bound_msa_pos;
		    if ($end_exon+1 < $num_exons) {
			my $j=$ExonStarts[$end_exon+1];
			while ($j<$msa_len) {
			    if ($MapMSA[$s1][$j]) {
				$right_bound = $MapMSA[$s1][$j];
				$right_bound_msa_pos = $j;
				last;
			    }
			    $j++;
			}
		    }
		    if (!$right_bound) {
			# Find the first coordinate for this species, noting it as such
			my $j=$msa_len-1;
			while (!$MapMSA[$s1][$j]) {
			    $j--;
			}
			$right_bound = '[end-of-coding-region:'.$MapMSA[$s1][$j].']';
			$right_bound_msa_pos = $j;
		    }


		    # [3] We'll need to disqualify any nucleotides currently identified
		    #     as contributing to our mapping, so we'll go exon-by-exon and
		    #     record any already-used regions.

		    my $used_regions_str;
		    for (my $exon_id=$start_exon; $exon_id<=$end_exon; $exon_id++) {

			my $region_start;
			for (my $j=$ExonStarts[$exon_id]; $j<$ExonEnds[$exon_id]; $j++) {
			    if ($MapMSA[$s1][$j]) {
				$region_start = $MapMSA[$s1][$j];
				last;
			    }
			}

			my $region_end;
			for (my $j=$ExonEnds[$exon_id]; $j>$ExonStarts[$exon_id]; $j--) {
			    if ($MapMSA[$s1][$j]) {
				$region_end = $MapMSA[$s1][$j];
				last;
			    }
			}

			if ($region_start && $region_end) {
			    if ($used_regions_str) {
				$used_regions_str = $used_regions_str.'&'.$region_start.'|'.$region_end;
			    } else {
				$used_regions_str = $region_start.'|'.$region_end;
			    }
			}
			
		    }

		    # Just to make sure we're sticking something into our list
		    if (!$used_regions_str) {
			$used_regions_str = '0|0';
		    }


		    # Heck yeah! Let's shout (and scream)!
		    push(@SearchSeqs,lc($search_seq_2));
		    push(@TargetSpecies,$species1);
		    push(@TargetSpeciesRange,$left_bound.'|'.$right_bound);
		    push(@SourceSpecies,$species2);
		    push(@MSAExonRanges,($start_exon+1).'..'.($end_exon+1));
		    push(@UsedRegions,$used_regions_str);
		    $num_ghost_exons++;
		    
		}

		
	    }

	}
	
    }

    # I tried so hard, and got so far, but in the end it still mattered just a lil' bit
    return (0,0) if ($num_ghost_exons == 0);
    
    # NOICE!  Time to get ready for some good 'n' nasty tblastn-ery!

    # We'll tally up the number of successes
    my $ghosts_busted = 0;

    # Prepare to write some stuff to a file!
    my $outf = OpenOutputFile($outgenesdir.$gene.'/search.out');

    # We'll just go hit-by-hit, because that's what's sensible.  'q' for query
    for (my $q=0; $q<$num_ghost_exons; $q++) {

	# Write out the protein sequence to our protein sequence file.
	# We don't use my fancy file functions here because we'll be overwriting.
	open(my $protf,'>',$prot_seq_fname);
	print $protf ">seq\n$SearchSeqs[$q]\n\n";
	close($protf);

	# What's our search chromosome? What are our nucleotide coords?
	my $target_species = $TargetSpecies[$q];
	my $chr = $SpeciesToChrs{$target_species};
	my $revcomp = 0;
	if ($chr =~ /\[revcomp\]/) {
	    $chr =~ s/\[revcomp\]//;
	    $revcomp = 1;
	}

	my @SearchRanges = split(/\|/,$TargetSpeciesRange[$q]);

	# If we're at either terminii of our sequence, pull in an extra 20k (or as
	# much as we can)
	if ($SearchRanges[0] =~ /\:(\d+)/) {
	    my $seq_start = $1;
	    if ($revcomp) {
		$SearchRanges[0] = Min($seq_start+20000,$ChrLensBySpecies{$target_species.'|'.$chr});		
	    } else {
		$SearchRanges[0] = Max($seq_start-20000,1);
	    }
	}
	
	if ($SearchRanges[1] =~ /\:(\d+)/) {
	    my $seq_end = $1;
	    if ($revcomp) {
		$SearchRanges[1] = Max($seq_end-20000,1);
	    } else {
		$SearchRanges[1] = Min($seq_end+20000,$ChrLensBySpecies{$target_species.'|'.$chr});
	    }
	}

	# Well, we sure know what sequence to pull in now, don't we!
	my $sfetch_cmd = $sfetch.' -range '.$SearchRanges[0].'..'.$SearchRanges[1];
	$sfetch_cmd = $sfetch_cmd.' '.$SpeciesToGenomes{$target_species}.' '.$chr;
	$sfetch_cmd = $sfetch_cmd.' > '.$nucl_seq_fname;
	RunSystemCommand($sfetch_cmd);

	# Because of the wonders of filename standardization, we can do this!
	RunSystemCommand($tblastn);

	# Grab the list of coordinates that have already been used for mapping this
	# sequence's analog (not the right word, maybe, but you know what I mean).
	my @PrevUsedRegions = split(/\&/,$UsedRegions[$q]);

	# What did we get?
	my @HitAminoStarts;
	my @HitAminoEnds;
	my @HitNuclStarts;
	my @HitNuclEnds;
	my @HitEVals;
	my $num_tbn_hits = 0;
	my $tbnf = OpenInputFile($tbn_out_fname);
	while (my $line = <$tbnf>) {
	    if ($line) {
		
		my @HitData = split(/\s+/,$line);
		my $amino_start = $HitData[6];
		my $amino_end   = $HitData[7];
		my $nucl_start  = $HitData[8];
		my $nucl_end    = $HitData[9];
		my $e_val       = $HitData[10];

		# We'll take notice of any hits with an e-val with a '-'
		next if ($e_val !~ /e\-/);

		# What are the actual nucleotide coords?
		if ($revcomp) {
		    $nucl_start = $SearchRanges[0] - $nucl_start;
		    $nucl_end   = $SearchRanges[0] - $nucl_end;
		} else {
		    $nucl_start += $SearchRanges[0];
		    $nucl_end   += $SearchRanges[0];
		}

		# We need to confirm that this isn't overlapping with any of
		# the nucleotide regions that have been previously used to
		# map this sequence.
		# It also can't fully contain an exon, ding-dong!
		my $was_prev_used = 0;
		foreach my $prev_used_region (@PrevUsedRegions) {

		    $prev_used_region =~ /^(\d+)\|(\d+)$/;
		    my $prev_start = $1;
		    my $prev_end = $2;

		    if ($revcomp) {
			if ($nucl_start <= $prev_start && $nucl_start >= $prev_end) {
			    $was_prev_used = 1;
			    last;
			}
			if ($nucl_end <= $prev_start && $nucl_end >= $prev_end) {
			    $was_prev_used = 1;
			    last;
			}
			if ($nucl_start >= $prev_start && $nucl_end <= $prev_end) {
			    $was_prev_used = 1;
			    last;
			}
		    } else {
			if ($nucl_start >= $prev_start && $nucl_start <= $prev_end) {
			    $was_prev_used = 1;
			    last;
			}
			if ($nucl_end >= $prev_start && $nucl_end <= $prev_end) {
			    $was_prev_used = 1;
			    last;
			}
			if ($nucl_start <= $prev_start && $nucl_end >= $prev_end) {
			    $was_prev_used = 1;
			    last;
			}
		    }

		}

		next if ($was_prev_used);

		push(@HitAminoStarts,$amino_start);
		push(@HitAminoEnds,$amino_end);
		push(@HitNuclStarts,$nucl_start);
		push(@HitNuclEnds,$nucl_end);
		push(@HitEVals,$e_val);
		$num_tbn_hits++;

	    }
	}
	close($tbnf);

	# For outputting, let's get the textual representation of direction
	# into the chromosome name (after recording the chromosome length of
	# our target species, for exon overlap checking).
	my $chr_len = $ChrLensBySpecies{$target_species.'|'.$chr};
	$chr = $chr.'[revcomp]' if ($revcomp);

	# Speaking of outputting, let's just have these tidbits on-hand, too
	my $source_species = $SourceSpecies[$q];

	$MSAExonRanges[$q] =~ /(\d+)\.\.(\d+)/;
	my $start_exon = $1;
	my $end_exon = $2;
	my $exon_str = 'Exon';
	if ($start_exon != $end_exon) {
	    $exon_str = $exon_str.'s '.$start_exon.'..'.$end_exon;
	} else {
	    $exon_str = $exon_str.' '.$end_exon;
	}

	my $target_info = "MSA Ali Region  : $exon_str\n";
	$target_info = $target_info."    Target Genome   : $target_species ($chr:$SearchRanges[0]..$SearchRanges[1])\n";
	$target_info = $target_info."    Source Species  : $source_species\n";
	
	# Is it an especially elusive ghost we're chasing?
	if ($num_tbn_hits == 0) {
	    print $outf "[ ] Search failure (no tblastn hits)\n";
	    print $outf "    $target_info";
	    print $outf "    Search Sequence : $SearchSeqs[$q]\n\n";
	    next;
	}
	
	# Oh, this is a most profitable Ghost Adventure indeed!
	$ghosts_busted++;

	# Let's illustrate how much of the sequence has been covered.
	# NOTE: We're assuming that our hits are consistent with one another,
	#   but we may want to double-check that in the future...
	my @MappedSeq = split(//,lc($SearchSeqs[$q]));
	for (my $hit=0; $hit<$num_tbn_hits; $hit++) {
	    for (my $pos=$HitAminoStarts[$hit]-1; $pos<$HitAminoEnds[$hit]; $pos++) {
		$MappedSeq[$pos] = uc($MappedSeq[$pos]);
	    }
	}
	my $mapped_seq = '';
	foreach my $char (@MappedSeq) {
	    $mapped_seq = $mapped_seq.$char;
	}

	# Sing it to high heaven!
	print $outf "[+] Search success!\n";
	print $outf "    $target_info";
	print $outf "    Search Sequence : $mapped_seq\n";
	print $outf "    Num tblastn Hits: $num_tbn_hits\n";

	# I'm going to take this 'underlining' out for now, and let the
	# upper / lower case distinction speak for itself.
	#
	#print $outf "    ";
	#foreach my $char (@MappedSeq) {
	#if ($char eq uc($char)) { print $outf '-'; }
	#else                    { print $outf ' '; }
	#}
	#print $outf "\n";

	for (my $hit=0; $hit<$num_tbn_hits; $hit++) {

	    print $outf "    + Aminos $HitAminoStarts[$hit]..$HitAminoEnds[$hit] ";
	    print $outf "mapped to $target_species $chr:$HitNuclStarts[$hit]..$HitNuclEnds[$hit] ";
	    print $outf "($HitEVals[$hit])\n";

	    # Does this hit overlap with a known (GTF-recorded) coding region?
	    my $start_nucl = $HitNuclStarts[$hit];
	    my $end_nucl = $HitNuclEnds[$hit];

	    my $keybase = $target_species.'|'.$chr.'|';
	    my $gtf_overlap = 0;

	    my $increment = 1;
	    $increment = -1 if ($revcomp);
	    for (my $i=$start_nucl; $i<=$end_nucl; $i+=$increment) {

		my $key = $keybase.$i;

		if ($GTFStartsToEnds{$key}) {

		    my @EndList;
		    if ($revcomp) {
			@EndList = sort { $a <=> $b } split(/\|/,$GTFStartsToEnds{$key});
		    } else {
			@EndList = sort { $b <=> $a } split(/\|/,$GTFStartsToEnds{$key});
		    }

		    $gtf_overlap = $i.'..'.$EndList[0];
		    last;
		    
		}

		if ($GTFEndsToStarts{$key}) {

		    my @StartList;
		    if ($revcomp) {
			@StartList = sort { $b <=> $a } split(/\|/,$GTFEndsToStarts{$key});
		    } else {
			@StartList = sort { $a <=> $b } split(/\|/,$GTFEndsToStarts{$key});
		    }
		    
		    $gtf_overlap = $StartList[0].'..'.$i;
		    last;
		    
		}
		
	    }

	    # Aw rats, overlap!
	    if ($gtf_overlap) {
		print $outf "      - Overlaps with annotated exon $chr:$gtf_overlap\n";
		next;
	    }

	    # One last thing that we'll do is scan backwards and forwards to the
	    # next annotated exon to make sure the range we've found isn't square
	    # in the middle of an exon.
	    if ($revcomp) {

		my $low_scan = $end_nucl-1;
		while ($low_scan) {
		    if ($GTFEndsToStarts{$low_scan}) {
			my @StartList = sort { $b <=> $a } split(/\|/,$GTFEndsToStarts{$low_scan});
			my $highest_start = $StartList[0];
			if ($highest_start > $end_nucl) {
			    $gtf_overlap = $highest_start.'..'.$low_scan;
			}
			last;
		    }
		    $low_scan--;
		}

	    } else {

		my $low_scan = $start_nucl-1;
		while ($low_scan) {
		    if ($GTFStartsToEnds{$low_scan}) {
			my @EndList = sort { $b <=> $a } split(/\|/,$GTFStartsToEnds{$low_scan});
			my $highest_end = $EndList[0];
			if ($highest_end > $start_nucl) {
			    $gtf_overlap = $low_scan.'..'.$highest_end;
			}
			last;
		    }
		    $low_scan--;
		}

	    }

	    # Aw Rats, OverLap!
	    if ($gtf_overlap) {
		print $outf "      - Overlaps with annotated exon $chr:$gtf_overlap\n";
		next;
	    }

	    # FINAL SCAN! Finding a high point and seeing if it's lowest corresponding
	    # low point is problematic (i.e., implies that this exon is known)
	    if ($revcomp) {

		my $high_scan = $start_nucl-1;
		while ($high_scan < $chr_len) {
		    if ($GTFStartsToEnds{$high_scan}) {
			my @EndList = sort { $a <=> $b } split(/\|/,$GTFStartsToEnds{$high_scan});
			my $lowest_end = $EndList[0];
			if ($lowest_end < $start_nucl) {
			    $gtf_overlap = $high_scan.'..'.$lowest_end;
			}
			last;
		    }
		    $high_scan++;;
		}

	    } else {

		my $high_scan = $end_nucl+1;
		while ($high_scan < $chr_len) {
		    if ($GTFEndsToStarts{$high_scan}) {
			my @StartList = sort { $a <=> $b } split(/\|/,$GTFEndsToStarts{$high_scan});
			my $lowest_start = $StartList[0];
			if ($lowest_start < $end_nucl) {
			    $gtf_overlap = $lowest_start.'..'.$high_scan;
			}
			last;
		    }
		    $high_scan++;
		}

	    }
	    
	    # AW RATS, OVERLAP!
	    if ($gtf_overlap) {
		print $outf "      - Overlaps with annotated exon $chr:$gtf_overlap\n";
		next;
	    }

	    # AW GOOD RATS, NO OVERLAP!!!
	    print $outf "      + No observed overlap with annotated exons\n";
	    
	}
	
	print $outf "\n";
	
    }

    close($outf);

    # How'd we do there?
    return ($num_ghost_exons, $ghosts_busted);

}






###############################################################
#
#  Function: FindAliQualityDrops
#
sub FindAliQualityDrops
{
    my $pcts_id_ref = shift;
    my $num_exons = shift;

    my @PctsID = @{$pcts_id_ref};

    # At what drop from a local maximum exon pair's quality
    # do we get concerned?
    my $q_drop_threshold = 25; # rec. 35+ for the alt. approach.

    # To avoid a bookkeeping headache, we'll just do two passes
    # (left-to-right, right-to-left, DUH).
    my @HasQDrop;
    for (my $i=0; $i<$num_exons; $i++) {
	$HasQDrop[$i] = 0;
    }

    my $exon_id = 0;
    while ($PctsID[$exon_id] == 0) {
	$exon_id++;
    }
    my $last_exon_pct_id = $PctsID[$exon_id++];

    while ($exon_id < $num_exons) {
	if ($PctsID[$exon_id]) {
	    if ($last_exon_pct_id - $PctsID[$exon_id] >= $q_drop_threshold) {
		$HasQDrop[$exon_id] = 1;
	    }
	    $last_exon_pct_id = $PctsID[$exon_id];
	}
	$exon_id++;
    }

    # Reverse course!
    # Note that exon_id starts out-of-bounds, so we decrement before anything else
    while ($exon_id) {
	$exon_id--;
	if ($PctsID[$exon_id]) {
	    if ($last_exon_pct_id - $PctsID[$exon_id] >= $q_drop_threshold) {
		$HasQDrop[$exon_id] = 1;
	    }
	    $last_exon_pct_id = $PctsID[$exon_id];
	}	
    }

    # We can imagine a case where we have something like 100/50/45/50/100,
    # but if we stick with what we have now, that middle exon isn't going
    # to be seen as low-quality.
    # To remedy this, we'll run along from each low-quality exon labeling any
    # exons that don't improve by >5% as also low-quality.
    # Again, to simplify this we'll just do a left-to-right and a right-to-left
    # scan.
    while ($exon_id < $num_exons) {
	if ($HasQDrop[$exon_id]) {
	    my $pct_to_beat = $PctsID[$exon_id] + 5.0;
	    $exon_id++;
	    while ($exon_id < $num_exons) {
		if ($PctsID[$exon_id]) {
		    if ($PctsID[$exon_id] < $pct_to_beat) {
			$HasQDrop[$exon_id] = 1;
		    } else {
			last;
		    }
		}
		$exon_id++;
	    }
	} else {
	    $exon_id++;
	}
    }

    while ($exon_id > 0) {
	$exon_id--;
	if ($HasQDrop[$exon_id]) {
	    my $pct_to_beat = $PctsID[$exon_id] + 5.0;
	    $exon_id--;
	    while ($exon_id >= 0) {
		if ($PctsID[$exon_id]) {
		    if ($PctsID[$exon_id] < $pct_to_beat) {
			$HasQDrop[$exon_id] = 1;
		    } else {
			last;
		    }
		}
		$exon_id--;
	    }
	}
    }

    # Phew! So close to being done!
    # Now we just need to get a list of the QDroppers and be on
    # our merry way.
    my @LowQExonPairs;
    for (my $i=0; $i<$num_exons; $i++) {
	if ($HasQDrop[$i]) {
	    push(@LowQExonPairs,$i);
	}
    }

    if (scalar(@LowQExonPairs) == 0) {
	return 0;
    }
    return \@LowQExonPairs;





    ###########################################
    #                                         #
    #  ALT APPROACH: DROP FROM LOCAL MAXIMUM  #
    #                                         #
    ###########################################
    #
    # NOTE: An issue with this method was that we would see some
    #       extremely high-quality exon alignments (e.g., 100%),
    #       from which the alignment quality would slowly walk down
    #       to something like 68%, thus flagging that lower-scoring
    #       exon even if it were flanked by 70s-ish exons.
    #


    # We'll start off by computing the avg. pct ID (ignoring 0s)
    # and initializing an array to store 'quality ratings'
    my $avg_pct_id = 0.0;
    my $num_nonzero = 0;
    my @QualityRatings;
    foreach my $pct_id (@PctsID) {
	push(@QualityRatings,0);
	if ($pct_id) {
	    $avg_pct_id += $pct_id;
	    $num_nonzero++;
	}
    }
    
    if ($num_nonzero) {
	$avg_pct_id /= $num_nonzero;
    } else {
	return 0;
    }

    # Next up, we'll give a '1' to every exon pair that has an
    # alignment quality equal to or above the average and a '-1'
    # to every pair below the average
    my @TopShelfExons;
    for (my $i=0; $i<$num_exons; $i++) {
	if ($PctsID[$i] >= $avg_pct_id) {
	    $QualityRatings[$i] = 1;
	} elsif ($PctsID[$i]) {
	    $QualityRatings[$i] = -1;
	}
    }

    # Find the local maxima.  Note that this isn't "every peak,"
    # but rather a maximum among a set of relatively strong exons.
    my @LocalMaxima;
    for (my $i=0; $i<$num_exons; $i++) {

        if ($QualityRatings[$i] == 1) {

	    my $max_pct = $PctsID[$i];
	    my $max_pos = $i;

	    $i++;
	    while ($i<$num_exons && $QualityRatings[$i] != -1) {
		if ($PctsID[$i] > $max_pct) {
		    $max_pct = $PctsID[$i];
		    $max_pos = $i;
		}
		$i++;
	    }

	    $QualityRatings[$max_pos] = 2;
	    push(@LocalMaxima,$max_pos);
	    
	}
	
    }


    # Now we can run in both directions from each local maximum
    # checking for drops below our q-drop threshold
    foreach my $max_pos (@LocalMaxima) {

	my $pct_id_cutoff = $PctsID[$max_pos] - $q_drop_threshold;

	# Scan left
	for (my $i=$max_pos-1; $i>=0; $i--) {
	    last if ($QualityRatings[$i] == 2);
	    if ($PctsID[$i] && $PctsID[$i] <= $pct_id_cutoff) {
		$QualityRatings[$i] = -2;
	    }
	}
	
	# Scan right
	for (my $i=$max_pos+1; $i<$num_exons; $i++) {
	    last if ($QualityRatings[$i] == 2);
	    if ($PctsID[$i] && $PctsID[$i] <= $pct_id_cutoff) {
		$QualityRatings[$i] = -2;
	    }
	}
	
    }

    # Finally, we can condense our below-threshold exons into a list
    my @QDroppers;
    for (my $i=0; $i<$num_exons; $i++) {
	push(@QDroppers,$i) if ($QualityRatings[$i] == -2);
    }

    return \@QDroppers;
    
}






###############################################################
#
#  Function:  RecordGhostMSAs
#
sub RecordGhostMSAs
{
    my $gene = shift;

    # First things first, we're going to need to grab ahold of this gene's dir
    my $genedir = ConfirmDirectory($outgenesdir.$gene);

    # Open up the file with all of the hits for this gene and read in all
    # successful maps
    my $inf = OpenInputFile($genedir.'search.out');

    my %TargetSpeciesToHits;
    while (my $line = <$inf>) {

	next if ($line !~ /Search success/);

	$line = <$inf>; # MSA position info.
	$line =~ /(Exons? \S+)/;
	my $source_exons = $1;
	
	$line = <$inf>; # Full target search region
	$line =~ /\: (\S+) \(([^\:]+)\:/;
	my $target_species = $1;
	my $target_chr = $2;
	my $revcomp = 0;
	$revcomp = 1 if ($target_chr =~ /\[revcomp\]/);

	$line = <$inf>; # The source (amino) species
	$line =~ /\: (\S+)/;
	my $source_species = $1;

	$line = <$inf>; # The sequence searched against the target genome
	$line =~ /\: (\S+)/;
	my $source_seq = uc($1);

	$line = <$inf>; # How many separate tblastn hits did we get?
	$line =~ /\: (\d+)/;
	my $num_tbn_hits = $1;

	my $wide_start = 0;
	my $wide_end = 0;
	my $novel_exon = 1;
	while ($num_tbn_hits) {

	    $line = <$inf>;
	    $line =~ /\:(\d+)\.\.(\d+)/;
	    my $hit_start = $1;
	    my $hit_end = $2;

	    if ($revcomp) {
		if (!$wide_start || $hit_start > $wide_start) {
		    $wide_start = $hit_start;
		}
		if (!$wide_end || $hit_end < $wide_end) {
		    $wide_end = $hit_end;
		}
	    } else {
		if (!$wide_start || $hit_start < $wide_start) {
		    $wide_start = $hit_start;
		}
		if (!$wide_end || $hit_end > $wide_end) {
		    $wide_end = $hit_end;
		}
	    }

	    $line = <$inf>;
	    $novel_exon = 0 if ($line =~ /\- Overlaps/);

	    $num_tbn_hits--;
	    
	}

	# Time to record this bad boi!
	my $hash_val = $target_chr.':'.$wide_start.'..'.$wide_end;
	$hash_val = $hash_val.'|'.$source_species.':'.$source_seq.':'.$source_exons;
	$hash_val = $hash_val.'|'.$novel_exon;

	if ($TargetSpeciesToHits{$target_species}) {
	    $TargetSpeciesToHits{$target_species} = $TargetSpeciesToHits{$target_species}.'&'.$hash_val;
	} else {
	    $TargetSpeciesToHits{$target_species} = $hash_val;
	}

    }

    close($inf);

    # Make an output directory for our alignment visualizations
    my $gene_ali_dir = CreateDirectory($genedir.'alignments');
    
    # Now we can run through our species actually building up some dang MSAs!
    foreach my $target_species (keys %TargetSpeciesToHits) {

	my @AllSpeciesHits = split(/\&/,$TargetSpeciesToHits{$target_species});

	# We'll need to make sure that our hits are split into exon groups
	my $num_exons = 0;
	my @ExonTargetSpecies;
	my @ExonChrs; # Just to be super sure!
	my @ExonNuclStarts;
	my @ExonNuclEnds;
	my @ExonHits;
	foreach my $hit (@AllSpeciesHits) {

	    $hit =~ /^([^\:]+)\:(\d+)\.\.(\d+)/;
	    my $hit_chr = $1;
	    my $hit_start = $2;
	    my $hit_end = $3;

	    $hit =~ /^[^\|]+\|([^\:]+)\:([^\|]+)\|/;
	    my $source_species = $1;
	    my $source_seq_str = $2;

	    my $revcomp = 0;
	    $revcomp = 1 if ($hit_chr =~ /\[revcomp\]/);

	    if (!$num_exons) {
		push(@ExonTargetSpecies,$target_species);
		push(@ExonChrs,$hit_chr);
		push(@ExonNuclStarts,$hit_start);
		push(@ExonNuclEnds,$hit_end);
		$ExonHits[0] = $hit;
		$num_exons++;
		next;
	    }

	    # All overlapping hits (w.r.t. the target species' genome) are grouped
	    # together as an "exon group"
	    my $exon_group = -1;
	    for (my $i=0; $i<$num_exons; $i++) {

		next if ($hit_chr ne $ExonChrs[$i]);

		if ($revcomp) {

		    if (($hit_start >= $ExonNuclStarts[$i] && $hit_end <= $ExonNuclStarts[$i])
			|| ($hit_end <= $ExonNuclEnds[$i] && $hit_start >= $ExonNuclEnds[$i])) {
			$exon_group = $i;
			$ExonNuclStarts[$i] = Max($hit_start,$ExonNuclStarts[$i]);
			$ExonNuclEnds[$i] = Min($hit_end,$ExonNuclEnds[$i]);
			last;
		    }

		} else {
		    
		    if (($hit_start <= $ExonNuclStarts[$i] && $hit_end >= $ExonNuclStarts[$i])
			|| ($hit_end >= $ExonNuclEnds[$i] && $hit_start <= $ExonNuclEnds[$i])) {
			$exon_group = $i;
			$ExonNuclStarts[$i] = Min($hit_start,$ExonNuclStarts[$i]);
			$ExonNuclEnds[$i] = Max($hit_end,$ExonNuclEnds[$i]);
			last;
		    }

		}
		
	    }

	    if ($exon_group == -1) {
		push(@ExonTargetSpecies,$target_species);
		push(@ExonChrs,$hit_chr);
		push(@ExonNuclStarts,$hit_start);
		push(@ExonNuclEnds,$hit_end);
		push(@ExonHits,$hit);
		$num_exons++;
	    } else {
		$ExonHits[$exon_group] = $ExonHits[$exon_group].'&'.$hit;
	    }
	    
	}

	
	# Them's some exons!  Now we can go through each exon-group and generate
	# an MSA representing the translated target sequence and each of the source
	# species' amino acid sequences
	#
	# Reminder: Each "exon" is a group of tblastn hits from "source" species to an
	#   an overlapping region of the "target" genome
	#
	# We'll write out a file with our MSA visualizations for each species
	my $outfname = $gene_ali_dir.$target_species.'.MSAs.out';
	my $outf = OpenOutputFile($outfname);
	for (my $i=0; $i<$num_exons; $i++) {

	    # Start off by getting access to the specific source species matches
	    my @SourceSpecies;
	    my @SourceSeqs;
	    my @SourceExons;
	    my $novel_exon = 1;
	    foreach my $hit (split(/\&/,$ExonHits[$i])) {

		$hit =~ /^[^\|]+\|([^\:]+)\:([^\:]+)\:([^\|]+)\|/;
		push(@SourceSpecies,$1);
		push(@SourceSeqs,$2);
		push(@SourceExons,lc($3));

		$hit =~ /\|(\d)$/;
		$novel_exon *= $1;

	    }
	    my $num_source_species = scalar(@SourceSpecies);

	    my $chr = $ExonChrs[$i];
	    my $revcomp = 0;
	    if ($chr =~ /\[revcomp\]/) {
		$chr =~ s/\[revcomp\]//;
		$revcomp = 1;
	    }
	    
	    my $search_start = $ExonNuclStarts[$i];
	    my $search_end = $ExonNuclEnds[$i];

	    # We'll pull in just  a lil' bit of extra genomic sequence, for the
	    # good of America.
	    my $extra_window = 10;
	    if ($revcomp) {
		$search_start += $extra_window;
		$search_end -= $extra_window;
	    } else {
		$search_start -= $extra_window;
		$search_end += $extra_window;
	    }

	    my $sfetch_cmd = $sfetch.' -range '.$search_start.'..'.$search_end;
	    $sfetch_cmd = $sfetch_cmd.' '.$SpeciesToGenomes{$target_species}.' '.$chr;
	    my $nucl_inf = OpenSystemCommand($sfetch_cmd);
	    my $header_line = <$nucl_inf>;
	    my $nucl_seq = '';
	    while (my $line = <$nucl_inf>) {
		$line =~ s/\n|\r//g;
		$nucl_seq = $nucl_seq.uc($line);
	    }
	    close($nucl_inf);

	    my @NuclSeq = split(//,$nucl_seq);

	    # Check which reading frame looks like it's the one we're supposed to be
	    # working with...
	    my $best_frame_num;
	    my $best_frame_score = 0;
	    my $best_frame_trans;
	    my @BestFrameStarts;
	    my @BestFrameEnds;
	    my @FrameScores;
	    my @FrameTranslations;
	    for (my $frame=0; $frame<3; $frame++) {

		# Pull in this reading frame
		my $frame_str = '';
		my $trans_str = '';
		for (my $i=$frame; $i+2<scalar(@NuclSeq); $i+=3) {

		    my $codon = $NuclSeq[$i].$NuclSeq[$i+1].$NuclSeq[$i+2];
		    $frame_str = $frame_str.$codon;

		    my $trans_aa = TranslateCodon($codon);
		    $trans_str = $trans_str.$trans_aa;

		}

		push(@FrameTranslations,$trans_str);

		# OUTDATED: Just get score from match / mismatch scoring
		## Perform a quick alignment to each source amino sequence, summing the
		## scores.
		#my $sum_score = 0;
		#foreach my $source_seq_str (@SourceSeqs) {
		#$sum_score += MatchMismatchScore($source_seq_str,$trans_str);
		#}

		# NEW APPROACH: Get the best local score along with start / end
		#   coordinates, relative to the search amino sequence.
		my $sum_score = 0;
		for (my $source_id=0; $source_id<scalar(@SourceSeqs); $source_id++) {

		    my ($lmm_score,$lmm_t_start,$lmm_t_end,$lmm_s_start,$lmm_s_end) =
			LocalMatchMismatchAli($trans_str,$SourceSeqs[$source_id]);

		    $sum_score += $lmm_score;

		    $BestFrameStarts[$frame][$source_id] = $lmm_s_start;
		    $BestFrameEnds[$frame][$source_id] = $lmm_s_end;
		    $FrameScores[$frame][$source_id] = $lmm_score;

		}

		if ($sum_score > $best_frame_score) {
		    $best_frame_score = $sum_score;
		    $best_frame_num = $frame;
		    $best_frame_trans = $trans_str;
		}
		
	    }

	    # I'm not sure why this is possible, but it seems that we have a
	    # handful of edge cases where we're slurping bits of non-ghost amino
	    # sequence into our searches, which leaves us searching for something
	    # we've already found in a place adjacent to (but not including)
	    # where we found it.
	    #
	    # That's a long-winded way of saying that sometimes our searches yield
	    # 0 good outputs, so we need to be able to catch cases where our best
	    # score is bad
	    next if ($best_frame_score <= 0);

	    # In case any of the sources didn't prefer the 'best' frame, we'll
	    # want to (1.) kick them off the team, and (2.) report that there
	    # might be something funky going on
	    my @MatchedSourceIDs;
	    my @UnmatchedSourceIDs;
	    my @UnmatchedFramePrefs;
	    for (my $source_id=0; $source_id<$num_source_species; $source_id++) {

		if ($FrameScores[$best_frame_num][$source_id] > 0) {

		    push(@MatchedSourceIDs,$source_id);

		} else {

		    # Which frame is better for this source amino sequence?
		    my $preferred_frame = 0;
		    if ($FrameScores[1] > $FrameScores[0]) {
			if ($FrameScores[1] > $FrameScores[2]) {
			    $preferred_frame = 1;
			} else {
			    $preferred_frame = 2;
			}
		    } elsif ($FrameScores[2] > $FrameScores[0]) {
			$preferred_frame = 2;
		    }

		    push(@UnmatchedSourceIDs,$source_id);
		    push(@UnmatchedFramePrefs,$preferred_frame);

		}
	    }

	    my $num_matched = scalar(@MatchedSourceIDs);
	    my $num_unmatched = scalar(@UnmatchedSourceIDs);

	    # Oh, dear, it looks like we need to annouce a disagreement on the
	    # proper frame...
	    if ($num_unmatched) {
		RecordFrameConflict($gene_ali_dir.'frame-disagreements.out',
				    $target_species,$chr,$revcomp,$search_start,
				    $search_end,$nucl_seq,$best_frame_num,
				    \@FrameTranslations,\@MatchedSourceIDs,
				    \@UnmatchedSourceIDs,\@UnmatchedFramePrefs,
				    \@SourceSpecies,\@SourceSeqs);
	    }

	    # Now that we have our best frame (and associated data) figured out,
	    # time to actually get alignin'!

	    # Starting off by priming with the first source sequence
	    my $source_id = $MatchedSourceIDs[0];
	    my @SourceSeqChars = split(//,$SourceSeqs[$source_id]);
	    my @AminoMSA;
	    for (my $char_id=$BestFrameStarts[$best_frame_num][$source_id];
		 $char_id<=$BestFrameEnds[$best_frame_num][$source_id];
		 $char_id++) {
		push(@AminoMSA,$SourceSeqChars[$char_id]);
	    }
	    
	    # And now for the rest of the crew...
	    for (my $meta_id=1; $meta_id<$num_matched; $meta_id++) {
		$source_id = $MatchedSourceIDs[$meta_id];
		@SourceSeqChars = split(//,$SourceSeqs[$source_id]);
		my @SourceAliChars;
		for (my $char_id=$BestFrameStarts[$best_frame_num][$source_id];
		     $char_id<=$BestFrameEnds[$best_frame_num][$source_id];
		     $char_id++) {
		    push(@SourceAliChars,$SourceSeqChars[$char_id]);
		}
		my $amino_msa_ref = MultiAminoSeqAli(\@AminoMSA,\@SourceAliChars);
		@AminoMSA = @{$amino_msa_ref};
	    }
	    
	    # We align the target sequence last so that it's (perhaps) more of an
	    # approximation of aligning to an "exon family profile"
	    my @TargetTrans = split(//,$best_frame_trans);

	    my $amino_msa_ref = MultiAminoSeqAli(\@TargetTrans,\@AminoMSA);
	    @AminoMSA = @{$amino_msa_ref};
	    my $amino_msa_len = scalar(@AminoMSA);

	    
	    # What are the actual nucleotide bounds of our putative coding region?
	    my $true_nucl_start = $search_start;
	    my $true_nucl_end;
	    if ($revcomp) {
		$true_nucl_start -= $best_frame_num;
		$true_nucl_end = $true_nucl_start+1 - (3 * length($best_frame_trans));
	    } else {
		$true_nucl_start += $best_frame_num;
		$true_nucl_end = $true_nucl_start-1 + (3 * length($best_frame_trans));
	    }

	    # If we have translated sequence aligned to nothing, we'll scrape it off
	    
	    # 1. Checking the left side
	    #
	    my $start_col=0;
	    while ($start_col<$amino_msa_len) {

		my @Col = split(//,$AminoMSA[$start_col]);

		my $trim_it = 1;
		for (my $i=1; $i<=$num_matched; $i++) {
		    if ($Col[$i] ne '-') {
			$trim_it = 0;
			last;
		    }
		}
		last if (!$trim_it);

		$start_col++;
		if ($revcomp) { $true_nucl_start -= 3; }
		else          { $true_nucl_start += 3; }

	    }

	    # If we didn't scrape anything off, we'll need to see if extending our
	    # nucleotide pull outwards makes sense
	    if ($start_col == 0) {
		for (my $col_id=0; $col_id<$amino_msa_len; $col_id++) {
		    if ($AminoMSA[$col_id] =~ /^\-/) {
			$AminoMSA[$col_id] =~ s/^\-/ /;
			if ($revcomp) { $true_nucl_start += 3; }
			else          { $true_nucl_start -= 3; }
		    } else {
			last;
		    }
		}		
	    }

	    # 2. Checking the right side
	    #
	    my $end_col=$amino_msa_len-1;
	    while ($end_col>=0) {

		my @Col = split(//,$AminoMSA[$end_col]);

		my $trim_it = 1;
		for (my $i=1; $i<=$num_matched; $i++) {
		    if ($Col[$i] ne '-') {
			$trim_it = 0;
			last;
		    }
		}
		last if (!$trim_it);

		$end_col--;
		if ($revcomp) { $true_nucl_end += 3; }
		else          { $true_nucl_end -= 3; }

	    }

	    # If we didn't scrape anything off, we'll need to see if extending our
	    # nucleotide pull outwards makes sense
	    if ($end_col == $amino_msa_len-1) {
		for (my $col_id=$amino_msa_len-1; $col_id>=0; $col_id--) {
		    if ($AminoMSA[$col_id] =~ /^\-/) {
			$AminoMSA[$col_id] =~ s/^\-/ /;
			if ($revcomp) { $true_nucl_end -= 3; }
			else          { $true_nucl_end += 3; }
		    } else {
			last;
		    }
		}
	    }

	    # Before we extend out, record the true start of the translated sequence
	    my $translation_start = $true_nucl_start;
	    my $translation_end   = $true_nucl_end;

	    # The last thing we're going to do is extend out 60 nucls on each side
	    # of the alignment...
	    if ($revcomp) {
		$true_nucl_start += 60;
		$true_nucl_end   -= 60;
	    } else {
		$true_nucl_start -= 60;
		$true_nucl_end   += 60;
	    }
	    
	    # Great!  Now that we have our final nucleotide region, let's grab
	    # those nucleotides.
	    $sfetch_cmd = $sfetch.' -range '.$true_nucl_start.'..'.$true_nucl_end;
	    $sfetch_cmd = $sfetch_cmd.' '.$SpeciesToGenomes{$target_species}.' '.$chr;
	    $nucl_inf = OpenSystemCommand($sfetch_cmd);
	    $header_line = <$nucl_inf>;
	    $nucl_seq = '';
	    while (my $line = <$nucl_inf>) {
		$line =~ s/\n|\r//g;
		next if (!$line);
		$nucl_seq = $nucl_seq.uc($line);
	    }
	    close($nucl_inf);
	    @NuclSeq = split(//,$nucl_seq);


	    # FINALLY TIME TO SKETCH OUR FINAL MSA
	    my @MSA;
	    my $msa_len=0;
	    
	    # 1. The lead-in nucleotides
	    my $nucl_seq_pos = 0;
	    while ($nucl_seq_pos < 60) {
		$MSA[0][$msa_len] = ' ';
		$MSA[1][$msa_len] = lc($NuclSeq[$nucl_seq_pos]);
		for (my $i=0; $i<$num_matched; $i++) {
		    $MSA[$i+2][$msa_len] = ' ';
		}
		$nucl_seq_pos++;
		$msa_len++;
	    }

	    # 2. The amino MSA
	    #    (Where we absolutely want to record %ID-able info)
	    my @SourceMatches;
	    my @SourceMismatches;
	    for (my $i=0; $i<$num_matched; $i++) {
		$SourceMatches[$i] = 0;
		$SourceMismatches[$i] = 0;
	    }
		 
	    for (my $col_id=$start_col; $col_id<=$end_col; $col_id++) {

		my @Col = split(//,$AminoMSA[$col_id]);

		# 2.a. The translated target amino sequence
		$MSA[0][$msa_len]   = ' ';
		$MSA[0][$msa_len+1] = $Col[0];
		$MSA[0][$msa_len+2] = ' ';
		
		# 2.b. The target nucleotides
		if ($Col[0] eq '-') {
		    $MSA[1][$msa_len]   = '-';
		    $MSA[1][$msa_len+1] = '-';
		    $MSA[1][$msa_len+2] = '-';
		} elsif ($Col[0] eq ' ') {
		    $MSA[1][$msa_len]   = lc($NuclSeq[$nucl_seq_pos++]);
		    $MSA[1][$msa_len+1] = lc($NuclSeq[$nucl_seq_pos++]);
		    $MSA[1][$msa_len+2] = lc($NuclSeq[$nucl_seq_pos++]);
		} else {
		    $MSA[1][$msa_len]   = $NuclSeq[$nucl_seq_pos++];
		    $MSA[1][$msa_len+1] = $NuclSeq[$nucl_seq_pos++];
		    $MSA[1][$msa_len+2] = $NuclSeq[$nucl_seq_pos++];
		}

		# 3.b. The source amino sequence(s)
		for (my $i=0; $i<$num_matched; $i++) {

		    $MSA[$i+2][$msa_len]   = ' ';

		    # Approach 1: Uppercase for matches, lowercase for mismatches
		    #$MSA[$i+2][$msa_len+1] = $Col[$i+1]; # Skip 0 (the target)
		    #$MSA[$i+2][$msa_len+1] = lc($Col[$i+1]) if ($Col[0] ne $Col[$i+1]);

		    # Approach 2: Periods for matches, lowercase for mismatches
		    if ($Col[$i+1] =~ /[A-Z]/ && $Col[$i+1] eq $Col[0]) {
			$MSA[$i+2][$msa_len+1] = '.';
			$SourceMatches[$i]++;
		    } else {
			$MSA[$i+2][$msa_len+1] = lc($Col[$i+1]);
			$SourceMismatches[$i]++;
		    }

		    $MSA[$i+2][$msa_len+2] = ' ';

		}

		# PROGRESS!
		$msa_len += 3;
		
	    }

	    # 3. The lead-out nucleotides
	    while ($nucl_seq_pos < scalar(@NuclSeq)) {
		$MSA[0][$msa_len] = ' ';
		$MSA[1][$msa_len] = lc($NuclSeq[$nucl_seq_pos]);
		for (my $i=0; $i<$num_matched; $i++) {
		    $MSA[$i+2][$msa_len] = ' ';
		}
		$nucl_seq_pos++;
		$msa_len++;		
	    }


	    # THAT'S IT!
	    # Now the only remaining work is the final formatting of the string!
	    my $longest_name_len = length($target_species);
	    foreach my $source_id (@MatchedSourceIDs) {
		my $species = $SourceSpecies[$source_id];
		if (length($species) > $longest_name_len) {
		    $longest_name_len = length($species);
		}
	    }
	    $longest_name_len += 4; # Two spaces on either side

	    my @FormattedNames;
	    $FormattedNames[0] = '  '.$target_species.'  ';
	    while (length($FormattedNames[0]) < $longest_name_len) {
		$FormattedNames[0] = ' '.$FormattedNames[0];
	    }

	    $FormattedNames[1] = ' ';
	    while (length($FormattedNames[1]) < $longest_name_len) {
		$FormattedNames[1] = ' '.$FormattedNames[1];
	    }

	    for (my $i=0; $i<$num_matched; $i++) {
		$FormattedNames[$i+2] = '  '.$SourceSpecies[$MatchedSourceIDs[$i]].'  ';
		while (length($FormattedNames[$i+2]) < $longest_name_len) {
		    $FormattedNames[$i+2] = ' '.$FormattedNames[$i+2];
		}
	    }

	    # Buffer in the alignment string and let 'er rip!
	    my $ali_str = "\n\n";
	    my $chars_per_line = 60;
	    my $msa_pos = 0;
	    while ($msa_pos < $msa_len) {

		my $next_stop = Min($msa_len,$msa_pos+$chars_per_line);

		for (my $i=0; $i<$num_matched+2; $i++) {

		    $ali_str = $ali_str.$FormattedNames[$i];

		    my $pos = $msa_pos;
		    while ($pos < $next_stop) {
			$ali_str = $ali_str.$MSA[$i][$pos++];
		    }
		    $ali_str = $ali_str."\n";

		}

		$ali_str = $ali_str."\n\n";

		$msa_pos += 60;
		
	    }
	    $ali_str = $ali_str."\n";

	    # Before we spit out our alignment string, we'll also make a string with
	    # hit metadata.

	    my @SourcePctsID;
	    my @SourceRatios;
	    for (my $i=0; $i<$num_matched; $i++) {

		my $ratio = $SourceMatches[$i]+$SourceMismatches[$i];
		my $pct_id = int(1000.0 * $SourceMatches[$i] / $ratio);
		$ratio = '('.$SourceMatches[$i].'/'.$ratio.')';

		# Formatting: Will we need to add a '.0'
		if ($pct_id % 10 == 0) {
		    $pct_id = $pct_id / 10.0;
		    $pct_id = $pct_id.'.0%';
		} else {
		    $pct_id = $pct_id / 10.0;
		    $pct_id = $pct_id.'%';
		}
		$pct_id = $pct_id.' alignment identity';

		push(@SourcePctsID,$pct_id);
		push(@SourceRatios,$ratio);
		
	    }

	    # Metadata item 1: Target sequence info.
	    my $meta_str = "\n";
	    $meta_str = $meta_str."  Target : $target_species $chr";
	    $meta_str = $meta_str.'[revcomp]' if ($revcomp);
	    $meta_str = $meta_str.":$translation_start..$translation_end\n";

	    if ($novel_exon) {
		$meta_str = $meta_str."         : Novel exon (no GTF overlaps)\n";
	    } else {
		$meta_str = $meta_str."         : Overlaps with GTF entry\n";
	    }

	    # Metadata item 2: Source sequence info.
	    $meta_str = $meta_str."  Source";
	    if ($num_matched > 1) { $meta_str = $meta_str.'s'; }
	    else                  { $meta_str = $meta_str.' '; }
	    $source_id = $MatchedSourceIDs[0];
	    $meta_str  = $meta_str.": $SourceSpecies[$source_id] $SourceExons[$source_id]\n";
	    $meta_str  = $meta_str."         : $SourcePctsID[0] $SourceRatios[0]\n";
	    for (my $i=1; $i<$num_matched; $i++) {
		$source_id = $MatchedSourceIDs[$i];
		$meta_str  = $meta_str."         : $SourceSpecies[$source_id] $SourceExons[$source_id]\n";
		$meta_str  = $meta_str."         : $SourcePctsID[$i] $SourceRatios[$i]\n";
	    }

	    # Print the alignment!!!
	    print $outf "\n\n-----------------------------------------------\n" if ($i);
	    print $outf "$meta_str";
	    print $outf "$ali_str";
	    
	}

	# We're officially done with this target species!
	close($outf);

	# Do a check to see if we actually reported any hits...
	if (!(-s $outfname)) { RunSystemCommand("rm \"$outfname\""); }

    }

}






###############################################################
#
#  Function:  MatchMismatchScore
#
sub MatchMismatchScore
{
    my $str1 = shift;
    my $str2 = shift;

    my @Seq1 = split(//,uc($str1));
    my $len1 = scalar(@Seq1);

    my @Seq2 = split(//,uc($str2));
    my $len2 = scalar(@Seq2);

    my @Matrix;
    for (my $i=0; $i<=$len1; $i++) { $Matrix[$i][0] = 0-$i; }
    for (my $j=0; $j<=$len2; $j++) { $Matrix[0][$j] = 0-$j; }

    for (my $i=1; $i<=$len1; $i++) {
	for (my $j=1; $j<=$len2; $j++) {
	    $Matrix[$i][$j]
		= Max(Max($Matrix[$i-1][$j],$Matrix[$i][$j-1])-1,$Matrix[$i-1][$j-1]);
	    if ($Seq1[$i-1] eq $Seq2[$j-1]) {
		$Matrix[$i][$j]++;
	    } else {
		$Matrix[$i][$j]--;
	    }
	}
    }

    return $Matrix[$len1][$len2];
    
}





###############################################################
#
#  Function:  LocalMatchMismatchAli
#
sub LocalMatchMismatchAli
{
    my $str1 = shift;
    my $str2 = shift;

    my @Seq1 = split(//,uc($str1));
    my $len1 = scalar(@Seq1);

    my @Seq2 = split(//,uc($str2));
    my $len2 = scalar(@Seq2);

    my @Matrix;
    for (my $i=0; $i<=$len1; $i++) { $Matrix[$i][0] = 0; }
    for (my $j=0; $j<=$len2; $j++) { $Matrix[0][$j] = 0; }

    my $mismatch = -1;
    my $match = 1;
    my $gap = -1;
    
    # First off, we'll find the highest scoring coordinate under a
    # local alignment scheme
    my $max_i;
    my $max_j;
    my $max_score = 0;
    for (my $i=1; $i<=$len1; $i++) {
	for (my $j=1; $j<=$len2; $j++) {

	    my $cell_score = $mismatch;
	    $cell_score = $match if ($Seq1[$i-1] eq $Seq2[$j-1]);
	    $cell_score += $Matrix[$i-1][$j-1];
	    
	    $Matrix[$i][$j] = Max(Max($Matrix[$i-1][$j],$Matrix[$i][$j-1])+$gap,
				  $cell_score);
	    $Matrix[$i][$j] = 0 if ($Matrix[$i][$j] < 0);

	    if ($max_score < $Matrix[$i][$j]) {
		$max_score = $Matrix[$i][$j];
		$max_i = $i;
		$max_j = $j;
	    }
	    
	}
    }

    # If the max score indicates that we weren't able to get a reasonably exon-y
    # alignment going, jump off
    if ($max_score < $match * 5) {
	return(-1,0,0,0,0);
    }

    # Next, we'll re-compute the top-left and bottom-right quadrants
    # of the matrix so that we have a global path that leads to our
    # highest-scoring region

    # Top-left
    for (my $i=1; $i<=$max_i; $i++) {
	for (my $j=1; $j<=$max_j; $j++) {

	    my $cell_score = $mismatch;
	    $cell_score = $match if ($Seq1[$i-1] eq $Seq2[$j-1]);
	    $cell_score += $Matrix[$i-1][$j-1];
	
	    $Matrix[$i][$j] = Max(Max($Matrix[$i-1][$j],$Matrix[$i][$j-1])+$gap,
				  $cell_score);
	    
	}
    }
    my $score_save = $Matrix[$max_i][$max_j]; # This cell gets overwritten

    # Bottom-right
    for (my $i=$max_i; $i<=$len1; $i++) { $Matrix[$i][$max_j] = ($i-$max_i) * $gap; }
    for (my $j=$max_j; $j<=$len2; $j++) { $Matrix[$max_i][$j] = ($j-$max_j) * $gap; }

    for (my $i=$max_i+1; $i<=$len1; $i++) {
	for (my $j=$max_j+1; $j<=$len2; $j++) {
	    
	    my $cell_score = $mismatch;
	    $cell_score = $match if ($Seq1[$i-1] eq $Seq2[$j-1]);
	    $cell_score += $Matrix[$i-1][$j-1];

	    $Matrix[$i][$j] = Max(Max($Matrix[$i-1][$j],$Matrix[$i][$j-1])+$gap,
				  $cell_score);
	    
	}
    }

    # Traceback!
    my @ITrace;
    my @JTrace;
    push(@ITrace,$len1);
    push(@JTrace,$len2);

    my $i=$len1;
    my $j=$len2;
    while ($i>$max_i && $j>$max_j) {

	my $cell_score = $mismatch;
	$cell_score = $match if ($Seq1[$i-1] eq $Seq2[$j-1]);
	$cell_score += $Matrix[$i-1][$j-1];
	
	if ($Matrix[$i][$j] == $cell_score ) {
	    $i--;
	    $j--;
	} elsif ($Matrix[$i][$j] == $Matrix[$i-1][$j]+$gap) {
	    $i--;
	} else {
	    $j--;
	}

	push(@ITrace,$i);
	push(@JTrace,$j);

    }

    while ($i>$max_i) {
	$i--;
	push(@ITrace,$i);
	push(@JTrace,$j);
    }

    while ($j>$max_j) {
	$j--;
	push(@ITrace,$i);
	push(@JTrace,$j);
    }

    # Halfway(-ish) there!  Time to trace our way back through the top-left quad
    $Matrix[$i][$j] = $score_save;

    while ($i && $j) {

	my $cell_score = $mismatch;
	$cell_score = $match if ($Seq1[$i-1] eq $Seq2[$j-1]);
	$cell_score += $Matrix[$i-1][$j-1];
	
	if ($Matrix[$i][$j] == $cell_score ) {
	    $i--;
	    $j--;
	} elsif ($Matrix[$i][$j] == $Matrix[$i-1][$j]+$gap) {
	    $i--;
	} else {
	    $j--;
	}

	push(@ITrace,$i);
	push(@JTrace,$j);

    }

    while ($i) {
	$i--;
	push(@ITrace,$i);
	push(@JTrace,$j);
    }

    while ($j) {
	$j--;
	push(@ITrace,$i);
	push(@JTrace,$j);
    }

    # Awesome!  Before we make any more SWEET progress, note that our traceback
    # is backwards, so we need to flip it around
    my $trace_len = scalar(@ITrace);
    for (my $pos=0; $pos<$trace_len/2; $pos++) {

	my $flip_pos = ($trace_len-1) - $pos;

	my $tmp = $ITrace[$pos];
	$ITrace[$pos] = $ITrace[$flip_pos];
	$ITrace[$flip_pos] = $tmp;

	$tmp = $JTrace[$pos];
	$JTrace[$pos] = $JTrace[$flip_pos];
	$JTrace[$flip_pos] = $tmp;

    }
    $trace_len--;

    # Determine the "score contribution" for each cell
    my @TraceScore;
    $TraceScore[0] = 0;
    my $key_pos; # Where do we pass through [max_i][max_j]?
    for (my $pos=1; $pos<$trace_len; $pos++) {

	# Knock this check out first
	if ($ITrace[$pos] == $max_i && $JTrace[$pos] == $max_j) {
	    $key_pos = $pos;
	}

	# What was this cell's contribution to the score of the maximum path?
	if ($ITrace[$pos] == $ITrace[$pos-1]+1 && $JTrace[$pos] == $JTrace[$pos-1]+1) {
	    if ($Seq1[$ITrace[$pos]-1] eq $Seq2[$JTrace[$pos]-1]) {
		$TraceScore[$pos] = $match;
	    } else {
		$TraceScore[$pos] = $mismatch;
	    }
	} else {
	    $TraceScore[$pos] = $gap;
	}
    }

    # We'll set our condition for killing the alignment as being a window of 8
    # aminos where 6 of the positions are gaps or mismatches.
    my $window_size = 8;
    my $min_matches = 2;
    my $kill_trigger # Scores below this terminate our walk
	= ($window_size-$min_matches) * Max($gap,$mismatch) + $min_matches*$match;

    # Scan left
    my $left_end_pos = $key_pos - int($window_size/2);
    if ($left_end_pos > 0 && $left_end_pos + $window_size < $trace_len) {

	my $window_score = 0;
	for (my $pos=0; $pos<$window_size; $pos++) {
	    $window_score += $TraceScore[$left_end_pos+$pos];
	}

	while ($window_score >= $kill_trigger && $left_end_pos) {
	    $left_end_pos--;
	    $window_score -= $TraceScore[$left_end_pos+$window_size];
	    $window_score += $TraceScore[$left_end_pos];
	}

    } else {
	$left_end_pos = 0;
    }

    # Scan right
    my $right_end_pos = $key_pos + int($window_size/2);
    if ($right_end_pos < $trace_len-1 && $right_end_pos - $window_size >= 0) {

	my $window_score = 0;
	for (my $pos=0; $pos<$window_size; $pos++) {
	    $window_score += $TraceScore[$right_end_pos-$pos];
	}

	while ($window_score >= $kill_trigger && $right_end_pos < $trace_len-1) {
	    $right_end_pos++;
	    $window_score -= $TraceScore[$right_end_pos-$window_size];
	    $window_score += $TraceScore[$right_end_pos];
	}
	
    } else {
	$right_end_pos = $trace_len-1;
    }

    # Eat inwards until we hit a match position
    while ($TraceScore[$left_end_pos] != $match) {
	$left_end_pos++;
    }
    while ($TraceScore[$right_end_pos] != $match) {
	$right_end_pos--;
    }

    # FINALLY!  Note that we're returning the original max local score,
    # which may not be representative of where we've trimmed the alignment.
    # NOTE that we need to reduce by 1 because our matrix corresponds to
    #   1-indexed sequences.
    return ($max_score,$ITrace[$left_end_pos]-1,$ITrace[$right_end_pos]-1,
	    $JTrace[$left_end_pos]-1,$JTrace[$right_end_pos]-1);

}






###############################################################
#
#  Function:  GetB62Score
#
sub GetB62Score
{
    my $seqstr1 = shift;
    my $seqstr2 = shift;

    # If either of our sequences involves a stop codon, we're unhappy
    if ($seqstr1 =~ /\*/ || $seqstr2 =~ /\*/) {
	return -100.0;
    }

    my @Chars1 = split(//,$seqstr1);
    my @Chars2 = split(//,$seqstr2);
    my $len1 = scalar(@Chars1);
    my $len2 = scalar(@Chars2);

    my $score = 0.0;
    my $score_mult = 1.0 / ($len1 * $len2);
    for (my $i=0; $i<$len1; $i++) {

	my $char1 = $Chars1[$i];
	next if (!$AminoIndex{$char1});
	$char1 = $AminoIndex{$char1} * 21;
	
	for (my $j=0; $j<$len2; $j++) {

	    my $char2 = $Chars2[$j];
	    next if (!$AminoIndex{$char2});
	    $char2 = $AminoIndex{$char2};

	    $score += $Blosum62[$char1+$char2] * $score_mult;

	}
    }

    return $score;
    
}






###############################################################
#
#  Function:  MultiAminoSeqAli
#
sub MultiAminoSeqAli
{
    my $seqs_1_ref = shift;
    my $seqs_2_ref = shift;

    my @Seqs1 = @{$seqs_1_ref};
    my @Seqs2 = @{$seqs_2_ref};

    my $len1 = scalar(@Seqs1);
    my $len2 = scalar(@Seqs2);

    # Let's be good and proper and use affine gapping
    my $gap_open = -1.5;
    my $gap_end  = -1.5;
    my $gap_ext  = -1.0;
    my @Match;
    my @HorizGap;
    my @VertGap;
    $Match[0][0]    = 0.0;
    $HorizGap[0][0] = $gap_open;
    $VertGap[0][0]  = $gap_open;
    for (my $i=1; $i<=$len1; $i++) {
	$HorizGap[$i][0] = $gap_open;
	$Match[$i][0]    = -100.0;
	$VertGap[$i][0]  = -100.0;
    }
    for (my $j=1; $j<=$len2; $j++) {
	$VertGap[0][$j]  = $gap_open;
	$Match[0][$j]    = -100.0;
	$HorizGap[0][$j] = -100.0;
    }

    for (my $i=1; $i<=$len1; $i++) {
	for (my $j=1; $j<=$len2; $j++) {

	    my $b62 = GetB62Score($Seqs1[$i-1],$Seqs2[$j-1]);

	    $Match[$i][$j] = Max(Max($HorizGap[$i-1][$j-1]+$gap_end+$b62,
				     $VertGap[$i-1][$j-1]+$gap_end+$b62),
				 $Match[$i-1][$j-1]+$b62);
	    
	    $HorizGap[$i][$j] = Max($HorizGap[$i-1][$j]+$gap_ext,
				    $Match[$i-1][$j]+$gap_open);

	    $VertGap[$i][$j] = Max($VertGap[$i][$j-1]+$gap_ext,
				   $Match[$i][$j-1]+$gap_open);
	    
	}
    }

    my $gapstr1 = '-';
    while (length($gapstr1) < length($Seqs1[0])) { $gapstr1 = $gapstr1.'-'; }

    my $gapstr2 = '-';
    while (length($gapstr2) < length($Seqs2[0])) { $gapstr2 = $gapstr2.'-'; }

    # During the traceback we'll need to know which state we're in
    my $s;
    if ($Match[$len1][$len2] > $HorizGap[$len1][$len2]) {
	if ($Match[$len1][$len2] > $VertGap[$len1][$len2]) {
	    $s='m';
	} else {
	    $s='v';
	}
    } elsif ($HorizGap[$len1][$len2] > $VertGap[$len1][$len2]) {
	$s='h';
    } else {
	$s='v';
    }

    # Time to back-trace!
    my @Ali;
    my $i=$len1;
    my $j=$len2;
    while ($i && $j) {

	if ($s eq 'm') {
	    
	    push(@Ali,$Seqs1[$i-1].$Seqs2[$j-1]);

	    my $b62 = GetB62Score($Seqs1[$i-1],$Seqs2[$j-1]);

	    if ($Match[$i][$j] == $Match[$i-1][$j-1]+$b62) {
		$s='m';
	    } elsif ($Match[$i][$j] == $HorizGap[$i-1][$j-1]+$gap_end+$b62) {
		$s='h';
	    } else {
		$s='v';
	    }

	    $i--;
	    $j--;

	} elsif ($s eq 'h') {

	    push(@Ali,$Seqs1[$i-1].$gapstr2);

	    if ($HorizGap[$i][$j] == $Match[$i-1][$j]+$gap_open) {
		$s='m';
	    } # else $s='h'
	    
	    $i--;

	} else { # $s=='v'

	    push(@Ali,$gapstr1.$Seqs2[$j-1]);

	    if ($VertGap[$i][$j] == $Match[$i][$j-1]+$gap_open) {
		$s='m';
	    } # else $s='v'

	    $j--;

	}

    }

    while ($i) {
	push(@Ali,$Seqs1[$i-1].$gapstr2);
	$i--;
    }

    while ($j) {
	push(@Ali,$gapstr1.$Seqs2[$j-1]);
	$j--;
    }

    # Uh-oh!  That alignment is BACKWARDS!!!
    for (my $i=0; $i<scalar(@Ali)/2; $i++) {
	my $tmp = $Ali[$i];
	$Ali[$i] = $Ali[scalar(@Ali)-1-$i];
	$Ali[scalar(@Ali)-1-$i] = $tmp;
    }

    return \@Ali;
    
}







###############################################################
#
#  Function:  RecordFrameConflict
#
sub RecordFrameConflict
{
    my $fname = shift;
    
    my $target_species = shift;
    my $chr = shift;
    my $revcomp = shift;
    my $search_start = shift;
    my $search_end = shift;
    my $nucl_seq = shift;
    my $best_frame_num = shift;
    my $frame_trans_ref = shift;
    
    my @FrameTranslations = @{$frame_trans_ref};

    my $matched_ids_ref = shift;
    my $unmatched_ids_ref = shift;
    my $unmatched_frames_ref = shift;
    my $source_species_ref = shift;
    my $source_seqs_ref = shift;

    my @MatchedSourceIDs = @{$matched_ids_ref};
    my @UnmatchedSourceIDs = @{$unmatched_ids_ref};
    my @UnmatchedFramePrefs = @{$unmatched_frames_ref};
    my @SourceSpecies = @{$source_species_ref};
    my @SourceSeqs = @{$source_seqs_ref};

    my $file_exists = 0;
    if (-e $fname) {
	$file_exists = 1;
    }

    open(my $outf,'>>',$fname) || die "\n  ERROR: Failed to open output file '$fname'\n\n";

    if ($file_exists) {
	print $outf "===========================================================\n\n\n";
    }
    
    if ($revcomp) {
	$chr = $chr.'[revcomp]';
    }

    print $outf "Target Species : $target_species\n";
    print $outf "Search Range   : $chr:$search_start..$search_end\n";
    print $outf "Nucleotides    : ";
    my @Nucls = split(//,$nucl_seq);
    for (my $i=0; $i<scalar(@Nucls); $i++) {
	print $outf "$Nucls[$i]";
	if ($i+1 < scalar(@Nucls) && ($i+1) % 60 == 0) {
	    print $outf "\n                 ";
	}
    }
    print $outf "\n\n";

    print $outf "MSA Frame ($best_frame_num)  : ";
    my @Seq = split(//,$FrameTranslations[$best_frame_num]);
    for (my $i=0; $i<scalar(@Seq); $i++) {
	print $outf "$Seq[$i]";
	if ($i+1 < scalar(@Seq) && ($i+1) % 60 == 0) {
	    print $outf "\n                 ";
	}
    }
    print $outf "\n\n";

    foreach my $source_id (@MatchedSourceIDs) {
	print $outf "                 + $SourceSpecies[$source_id]\n";
    }
    print $outf "\n";
    
    my %PrefFramesToSpecies;
    my %SpeciesToIDs;
    for (my $i=0; $i<scalar(@UnmatchedFramePrefs); $i++) {
	my $pref_frame = $UnmatchedFramePrefs[$i];
	my $source_species = $SourceSpecies[$UnmatchedSourceIDs[$i]];
	$SpeciesToIDs{$source_species} = $i+1;
	if ($PrefFramesToSpecies{$pref_frame}) {
	    $PrefFramesToSpecies{$pref_frame} = $PrefFramesToSpecies{$pref_frame}.','.$source_species;
	} else {
	    $PrefFramesToSpecies{$pref_frame} = $source_species;
	}
    }

    for (my $frame=0; $frame<3; $frame++) {

	next if ($frame == $best_frame_num);

	print $outf "    Frame  $frame   : ";
	my @Seq = split(//,$FrameTranslations[$frame]);
	for (my $i=0; $i<scalar(@Seq); $i++) {
	    print $outf "$Seq[$i]";
	    if ($i+1 < scalar(@Seq) && ($i+1) % 60 == 0) {
		print $outf "\n                 ";
	    }
	}
	print $outf "\n\n";

	if (!$PrefFramesToSpecies{$frame}) {
	    print $outf "                 - Not preferred by any source species\n\n";
	    next;
	}

	foreach my $source_species (split(/\,/,$PrefFramesToSpecies{$frame})) {
	    my $source_id = $SpeciesToIDs{$source_species}-1;
	    print $outf "                 + $source_species\n";
	    print $outf "                   $SourceSeqs[$source_id]\n\n";
	}
	
    }
    
    print $outf "\n\n";
    close($outf);
    
}






###############################################################
#
#  Function:  GetMapSummaryStats
#
sub GetMapSummaryStats
{

    my $ghostlygenes_ref = shift;
    my @GhostlyGenes = @{$ghostlygenes_ref};

    my $outf = OpenOutputFile($outdirname.'Search-Summary.out');

    my %TargetSpeciesToSuggested;
    my %TargetSpeciesToAnnotated;
    my @FullHitList;
    
    foreach my $gene (sort @GhostlyGenes) {
	
	my $infname = $outgenesdir.$gene.'/search.out';
	my $inf = OpenInputFile($infname);
	
	print $outf "\n  $gene\n";
	
	# We'll want to know how many ranges had mappings for this gene
	# (e.g., distinct blocks in the MSA).
	my $num_mapped_ranges = 0;
	my %MapsByExonRange;
	while (my $line = <$inf>) {
	    
	    # The sure-fire sign of a good time!
	    if ($line =~ /\[\+\]/) {
		
		# Where in the MSA is this mapping from?
		$line = <$inf>;
		$line =~ /Exons? (\S+)/;
		my $exon_range = $1;
		
		# Which genome was being searched against?
		$line = <$inf>;
		$line =~ /\: (\S+) \((\S+)\)\s*$/;
		my $target_species = $1;
		my $target_region = $2;
		
		# Which species provided the protein sequence?
		$line = <$inf>;
		$line =~ /\: (\S+)/;
		my $source_species = $1;
		
		# What was the protein sequence we were searching with?
		# Note that mapped residues are uppercase, so we can
		#   compute the percentage of the search sequence that
		#   was mapped.
		$line = <$inf>;
		$line =~ /\: (\S+)/;
		my $num_mapped_chars = 0;
		my $num_unmapped_chars = 0;
		foreach my $char (split(//,$line)) {
		    if ($char =~ /[A-Z]/) { $num_mapped_chars++;   }
		    else                  { $num_unmapped_chars++; }
		}
		my $search_seq_len = $num_mapped_chars + $num_unmapped_chars;
		my $pct_seq_mapped = int(1000.0 * $num_mapped_chars / $search_seq_len) / 10.0;
		
		# tblastn hit count? (this plays into the final part of our hash val)
		$line = <$inf>;
		$line =~ /\: (\d+)/;
		my $num_hits = $1;
		
		# We'll store this data in a hash, which we can then sort by
		# MSA position (further broken down by target species).
		my $hash_val = $target_species.'|'.$target_region.'|'.$source_species;
		$hash_val = $hash_val.'|'.$search_seq_len.'|'.$pct_seq_mapped;
		
		# Add the list of hit regions to 
		for (my $i=0; $i<$num_hits; $i++) {
		    $line = <$inf>;
		    $line =~ /mapped to \S+ \S+\:(\d+\.\.\d+)/;
		    my $map_range = $1;
		    $line = <$inf>;
		    my $known_exon_overlap = 1;
		    $known_exon_overlap = 0 if ($line =~ /\+ No observed overlap/);
		    $hash_val = $hash_val.'|'.$map_range.'/'.$known_exon_overlap;
		}
		
		if ($MapsByExonRange{$exon_range}) {
		    $MapsByExonRange{$exon_range} = $MapsByExonRange{$exon_range}.'&'.$hash_val;
		} else {
		    $MapsByExonRange{$exon_range} = $hash_val;
		    $num_mapped_ranges++;
		}
		
	    }
	    
	}
	
	close($inf);
	
	# Print out the number of mapped ranges 
	print $outf "  $num_mapped_ranges exon range";
	print $outf "s" if ($num_mapped_ranges > 1);
	print $outf " in MSA produced ghost exon mappings\n";
	
	# We now have all the info. we could need to grab from the file,
	# so it's time to organize it!
	foreach my $exon_range (sort keys %MapsByExonRange) {
	    
	    my %MapsByTargetSpecies;
	    my $num_targets_in_range = 0;
	    foreach my $map_in_range (split(/\&/,$MapsByExonRange{$exon_range})) {

		$map_in_range =~ /^([^\|]+)\|/;
		my $species = $1;
		
		if ($MapsByTargetSpecies{$species}) {
		    $MapsByTargetSpecies{$species} = $MapsByTargetSpecies{$species}.'&'.$map_in_range;
		} else {
		    $MapsByTargetSpecies{$species} = $map_in_range;
		    $num_targets_in_range++;
		}
		
	    }
	    
	    # Announce this range
	    print $outf "\n    > MSA exon";
	    print $outf "s" if ($exon_range =~ /\.\./);
	    print $outf " $exon_range: Found mappings to $num_targets_in_range genome";
	    print $outf "s" if ($num_targets_in_range > 1);
	    print $outf "\n";
	    
	    # We can now consider searches within this exon range by order
	    # of target species.
	    foreach my $target_species (sort keys %MapsByTargetSpecies) {
		
		my @MapsToSpecies = split(/\&/,$MapsByTargetSpecies{$target_species});
		
		# We'll go ahead and grab the region of the target species' genome that
		# was mapped to (should be consistent...)
		$MapsToSpecies[0] =~ /^[^\|]+\|([^\|]+)\|/;
		my $target_species_region = $1;
		
		my $num_maps_to_species = scalar(@MapsToSpecies);
		
		# As is customary in this region, we must declare the number of mappings
		# made to this species' genome.
		print $outf "\n        $num_maps_to_species mapping";
		print $outf "s" if ($num_maps_to_species > 1);
		print $outf " to target species $target_species ($target_species_region)\n";
		
		# We'll want to record all of the genomic coordinate ranges that we hit
		# to so that we can infer the number of putative exons suggested by the
		# sum of mappings.
		my @RangeList;
		my @KnownExon;
		
		foreach my $map_to_species (@MapsToSpecies) {
		    
		    # Recover the details of our mapping
		    my @MappingDetails = split(/\|/,$map_to_species);
		    my $source_species = $MappingDetails[2];
		    my $search_seq_len = $MappingDetails[3];
		    my $pct_seq_mapped = $MappingDetails[4];
		    
		    # Pull in all of the genomic ranges we mapped to
		    for (my $i=5; $i<scalar(@MappingDetails); $i++) {
			$MappingDetails[$i] =~ /^([^\/]+)\/(\d)$/;
			push(@RangeList,$1);
			push(@KnownExon,$2);
		    }
		    
		    # Yell about 'em!
		    print $outf "          + $source_species:";
		    print $outf " $pct_seq_mapped\% of $search_seq_len-amino search sequence mapped\n";
		    
		}
		
		# Before we move along, we'll use the genomic coordinates to determine
		# how many distinct putative exons are suggested by our data
		my ($num_suggested_exons,$num_annotated_exons)
		    = CollapseAndCountOverlaps(\@RangeList,\@KnownExon);
		print $outf "          = $num_suggested_exons unique exon";
		print $outf "s" if ($num_suggested_exons > 1);
		print $outf " suggested ($num_annotated_exons ";
		if ($num_annotated_exons == 1) { print $outf "has";  }
		else                           { print $outf "have"; }
		print $outf " GTF annotations)\n";

		if ($TargetSpeciesToSuggested{$target_species}) {
		    $TargetSpeciesToSuggested{$target_species} += $num_suggested_exons;
		} else {
		    $TargetSpeciesToSuggested{$target_species}  = $num_suggested_exons;
		}

		if ($TargetSpeciesToAnnotated{$target_species}) {
		    $TargetSpeciesToAnnotated{$target_species} += $num_annotated_exons;
		} else {
		    $TargetSpeciesToAnnotated{$target_species}  = $num_annotated_exons;
		}

		# Finally, record this gene's hit ratio
		my $hit_info = $target_species.'|'.$gene.'|'.$num_suggested_exons.'|'.$num_annotated_exons;
		push(@FullHitList,$hit_info);
		
	    }
	    
	}
	
    }

    close($outf);

    # Now we'll run through each of our target species and give a little bit of info.
    foreach my $species (keys %TargetSpeciesToSuggested) {

	$outf = OpenOutputFile($outdirname.$species.'-summary.out');

	my $num_suggested_exons = $TargetSpeciesToSuggested{$species};
	my $num_annotated_exons = 0;
	$num_annotated_exons = $TargetSpeciesToAnnotated{$species};

	print $outf "Species: $species\n";
	print $outf "Total tblastn-suggested Exons: $num_suggested_exons\n";
	print $outf "Number of GTF-annotated Exons: $num_annotated_exons\n";
	print $outf "\n";

	# To properly format our output, we're going to need to know the
	# longest gene name...
	my $longest_gene_name = 4; # 'Gene'
	foreach my $hit (@FullHitList) {

	    my @HitData = split(/\|/,$hit);
	    next if ($HitData[0] ne $species);

	    if (length($HitData[1]) > $longest_gene_name) {
		$longest_gene_name = length($HitData[1]);
	    }

	}

	my $fmted_header_1 = 'Gene';
	my $fmted_header_2 = '----';
	while (length($fmted_header_1) < $longest_gene_name) {
	    $fmted_header_1 = $fmted_header_1.' ';
	    $fmted_header_2 = $fmted_header_2.'-';
	}
	print $outf "$fmted_header_1  # Suggested  # Annotated\n";
	print $outf "$fmted_header_2  -----------  -----------\n";

	foreach my $hit (@FullHitList) {

	    my @HitData = split(/\|/,$hit);
	    next if ($HitData[0] ne $species);

	    my $field_1 = $HitData[1];
	    while (length($field_1) < $longest_gene_name) {
		$field_1 = $field_1.' ';
	    }

	    my $field_2 = $HitData[2];
	    while (length($field_2) < length('# Suggested')) {
		$field_2 = ' '.$field_2;
	    }

	    my $field_3 = $HitData[3];
	    while (length($field_3) < length('# Annotated')) {
		$field_3 = ' '.$field_3;
	    }

	    print $outf "$field_1  $field_2  $field_3\n";
	    
	}

	print $outf "\n";

	close($outf);
	
    }
    
}






#################################################################
#
#  Function:  CollapseAndCountOverlaps
#
sub CollapseAndCountOverlaps
{
    my $range_list_ref = shift;
    my @RangeList = @{$range_list_ref};

    my $known_exon_ref = shift;
    my @KnownExon = @{$known_exon_ref};
    
    my @RangeStarts;
    my @RangeEnds;
    foreach my $range (@RangeList) {
	$range =~ /^(\d+)\.\.(\d+)$/;
	push(@RangeStarts,$1);
	push(@RangeEnds,$2);
    }
    
    my $num_ranges = scalar(@RangeList);
    
    # Because we're just identifying overlaps, we can swap things around
    # so that start < end (i.e., de-revcomp)
    if ($RangeStarts[0] > $RangeEnds[0]) {
	for (my $i=0; $i<$num_ranges; $i++) {
	    my $temp = $RangeStarts[$i];
	    $RangeStarts[$i] = $RangeEnds[$i];
	    $RangeEnds[$i] = $temp;
	}
    }
    
    my %StartsToEnds;
    my %StartsToKnownExons;
    for (my $i=0; $i<$num_ranges; $i++) {
	my $start = $RangeStarts[$i];
	my $end = $RangeEnds[$i];
	if (!$StartsToEnds{$start}
	    || ($StartsToEnds{$start} && $StartsToEnds{$start} < $end)) {
	    $StartsToEnds{$start} = $end;
	}
	if ($KnownExon[$i]) {
	    $StartsToKnownExons{$start} = 1;
	}
    }
    
    my @CollapsedStarts = sort { $a <=> $b } keys %StartsToEnds;
    my @CollapsedEnds;
    foreach my $start (@CollapsedStarts) {
	push(@CollapsedEnds,$StartsToEnds{$start});
    }

    my $num_known = 0;
    my $is_known = 0;
    $is_known = 1 if ($StartsToKnownExons{$CollapsedStarts[0]});
    my $num_exons = 1;
    my $current_end = $CollapsedEnds[0];
    for (my $i=1; $i<scalar(@CollapsedEnds); $i++) {
	if ($current_end < $CollapsedStarts[$i]) {
	    $num_exons++;
	    $num_known += $is_known;
	    $is_known = 0;
	}
	$is_known = 1 if ($StartsToKnownExons{$CollapsedStarts[$i]});
	$current_end = Max($current_end,$CollapsedEnds[$i]);
    }

    $num_known += $is_known;

    return ($num_exons,$num_known);
    
}






# EOF

#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;
use Cwd;
use Getopt::Long;
use Time::HiRes;

# I AM UNAVOIDABLE
sub GetThisDir { my $lib = $0; $lib =~ s/\/Oasis.pl$//; return $lib; }
use lib GetThisDir();
use BureaucracyMirage;
use DisplayProgress;

# Subroutines
sub ParseArgs;
sub GetMappedSeqMSA;
sub ParseAFA;
sub RecordSplicedMSA;
sub ReduceMSAToSpecies;
sub FindGhostExons;
sub FindAliQualityDrops;





##############
#            #
#   SCRIPT   #
#            #
##############



if (@ARGV < 2) { die "\n  USAGE:  ./Oasis.pl [Mirage-Results] [Species-Guide]\n\n"; }



# Figure out what the location of the Mirage src directory is
my $location = $0;
$location =~ s/Oasis\.pl$//;

# We're going to need these friends
my $sindex = $location.'../inc/hsi/sindex';
my $sfetch = $location.'../inc/hsi/sfetch';
my $sstat  = $location.'../inc/hsi/sstat';

# Another friend we'll need is blat... but which one?!
my $blat = $location.'../inc/blat/';
my $UnameCmd = OpenSystemCommand('uname -a |');
my $uname = <$UnameCmd>;
close($UnameCmd);
if    (uc($uname) =~ /^LINUX /)  { $blat = $blat.'blat.linux.x86_64';  }
elsif (uc($uname) =~ /^DARWIN /) { $blat = $blat.'blat.macOSX.x86_64'; }
else                             { $blat = $blat.'blat.macOSX.i386';   }


# TODO: Make these options available as commandline arguments
my $options_ref = ParseArgs();
my %Options = %{$options_ref};
my $num_cpus = $Options{cpus};
my $outdirname = CreateDirectory($Options{outdirname});
my $save_msas = $Options{savemsas}; # Do we want to write our spliced MSAs to files?
my $bad_ali_cutoff = $Options{alicutoff};

# Confirm that the input directory looks like the real deal
my $input_dirname = ConfirmDirectory($ARGV[0]);
my $final_results_dirname = ConfirmDirectory($input_dirname.'Final-MSAs');
my $all_species_dirname = ConfirmDirectory($input_dirname.'Species-MSAs');


# An astute observer will notice that these aren't the same settings as Quilter
# uses, which is because this isn't frickin' Quilter, geez.
# It's rad that our standardized filenames let us play like this!
$blat = $blat.' -tileSize=5 -minIdentity=90 -maxIntron=0';
$blat = $blat.' -t=dnax -q=prot -out=blast8 1>/dev/null 2>&1';


# Start off the real work by parsing the species guide, which will give
# us genome locations and chromosome lengths.
my $tildedir_check = OpenSystemCommand('echo ~');
my $tildedir = <$tildedir_check>;
$tildedir =~ s/\n|\r//g;
$tildedir = ConfirmDirectory($tildedir);
close($tildedir_check);
my $SpeciesGuide = OpenInputFile($ARGV[1]);
my %SpeciesToGenomes;
my %ChrLensBySpecies;
while (my $line = <$SpeciesGuide>) {

    $line =~ s/\n|\r//g;
    next if ($line !~ /(\S+)\s+(\S+)\s+\S+/);
    my $species = lc($1);
    my $genome = $2;
    $genome =~ s/\~\//$tildedir/;

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
# fill in all of the wild 'n' wacky blat arguments we'll be using).
my $nucl_seq_fname = $outdirname.'nucl.tmp'.$threadID.'.fa';
my $prot_seq_fname = $outdirname.'prot.tmp'.$threadID.'.fa';
my $blat_out_fname = $outdirname.'blat.tmp'.$threadID.'.out';
$blat = $blat.' '.$nucl_seq_fname.' '.$prot_seq_fname.' '.$blat_out_fname;


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

    # Do we want to write this out to a spliced msa file?
    if ($save_msas) {
	RecordSplicedMSA($spliced_dirname.$gene.'.afa',\@MSA,\@SeqNames,$num_seqs,$msa_len);
    }

    # Now we'll reduce our MSA even further, down to just one sequence per
    # species!
    ($msa_ref,$mapmsa_ref,$seqnames_ref,$num_seqs)
	= ReduceMSAToSpecies(\@MSA,\@MapMSA,\@SeqNames,$num_seqs,$msa_len);
    @MSA = @{$msa_ref};
    @MapMSA = @{$mapmsa_ref};
    my @SpeciesNames = @{$seqnames_ref};
    my $num_species = $num_seqs;
    
    # Now that we have our super-reduced splice-site-ified MSA, let's get real nasty
    # with it (by way of locating exons suggestive of "ghosts")!
    my ($num_ghost_exons, $num_ghosts_busted) =
	FindGhostExons($gene,\@MSA,\@MapMSA,\@SpeciesNames,$num_species,$msa_len,\%SpeciesToChrs);

    # Add onto our overarching tallies
    $total_ghost_exons += $num_ghost_exons;
    $total_ghosts_busted += $num_ghosts_busted;

    # If we busted any ghosts, take note of this gene
    push(@GhostlyGenes,$gene) if ($num_ghosts_busted);
    
}

# How'd I do?  I don't even know!
if ($threadID) {
    my $final_outf = OpenOutputFile($outdirname.$threadID.'.final-tally.out');
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
    $nucl_seq_fname = $outdirname.'nucl.tmp'.$threadID.'.fa';
    $prot_seq_fname = $outdirname.'prot.tmp'.$threadID.'.fa';
    $blat_out_fname = $outdirname.'blat.tmp'.$threadID.'.out';

    if (-e $nucl_seq_fname) { system("rm $nucl_seq_fname"); }
    if (-e $prot_seq_fname) { system("rm $prot_seq_fname"); }
    if (-e $blat_out_fname) { system("rm $blat_out_fname"); }

    # How'd ya do, helper?
    if ($threadID) {
	
	my $final_infname = $outdirname.$threadID.'.final-tally.out';
	my $final_inf = OpenInputFile($final_infname);
    
	my $line = <$final_inf>;
	$line =~ /(\d+) \/ (\d+)/;
	$total_ghosts_busted += $1;
	$total_ghost_exons += 2;
	
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
    print "\n  No potentially-unannotated exons detected\n\n";

} else {

    # We'll list off all of our "busted" genes in a special file
    my $final_outfname = $outdirname.'Genes-With-Ghost-Exons.out';
    my $final_outf = OpenOutputFile($final_outfname);
    my $num_busted_genes = scalar(@GhostlyGenes);
    print $final_outf "$num_busted_genes genes with at least one mapped ghost exon:\n";
    foreach my $gene (sort(@GhostlyGenes)) {
	print $final_outf "  $gene\n";
    }
    print $final_outf "\n";
    close($final_outf);

    # WOOOOOOO, WE FOUND AT LEAST ONE THING TO POSSIBLY NOT CROSS-MAP!
    my $bust_rate = int(1000.0*$total_ghosts_busted/$total_ghost_exons)/10;
    print "\n";
    print "  $total_ghost_exons possible unannotated exons detected,\n";
    print "  $total_ghosts_busted of which have some sequence-level support ($bust_rate\%)\n";
    print "\n";
    print "  Results in '$outdirname'\n";
    print "  List of genes with successfully mapped ghost exons in '$final_outfname'\n";
    print "\n";

}


1;






###################
#                 #
#   SUBROUTINES   #
#                 #
###################






###############################################################
#
#  Function: ParseArgs
#
sub ParseArgs
{

    my %Options = (
	cpus => 1,
        outdirname => 'Oasis-Results',
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
    # Note that this will clear any original '*' columns (if we're
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
    # TODO: Maybe have an option to spit out a 'dirtified' MSA to a file, once
    #       we've drawn in our splice site positions...
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
	    $run_end++;
	    $run_votes += $VetoedColsVotes[$run_end];
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
		$SplicedMSA[$i][$spliced_msa_len] = '*';
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
	$SplicedMSA[$i][$spliced_msa_len] = '*';
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

    my $outfname = shift;
    my $msa_ref = shift;
    my $seqnames_ref = shift;
    my $num_seqs = shift;
    my $msa_len = shift;

    my @MSA = @{$msa_ref};
    my @SeqNames = @{$seqnames_ref};
    
    # Duplicates are *not* invited to my birthday party
    RunSystemCommand("rm $outfname") if (-e $outfname);
    
    my $outf = OpenOutputFile($outfname);
    for (my $i=0; $i<$num_seqs; $i++) {

	print $outf ">$SeqNames[$i]\n";
	for (my $j=0; $j<$msa_len; $j++) {
	    print $outf "$MSA[$i][$j]";
	    print $outf "\n" if (($j+1) % 60 == 0);
	}
	print $outf "\n" if ($msa_len % 60);
	print $outf "\n";
	
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
	    if ($MSA[0][$j] eq '*') {
		$SpeciesMSA[$num_species][$j] = '*';
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
	if ($MSA[0][$j] eq '*') {
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
    # Blat-ing.
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
	    # that we'll need to perform our Blat searches.
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
    
    # NOICE!  Time to get ready for some good 'n' nasty blattery!
    # (blattery will get you everywhere)

    # We'll tally up the number of successes
    my $ghosts_busted = 0;

    # We'll just go hit-by-hit, because that's what's sensible.  'q' for query
    my $outf = OpenOutputFile($outdirname.$gene.'.out');
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
	RunSystemCommand($blat);

	# Grab the list of coordinates that have already been used for mapping this
	# sequence's analog (not the right word, maybe, but you know what I mean).
	my @PrevUsedRegions = split(/\&/,$UsedRegions[$q]);

	# What did we get?
	my @HitAminoStarts;
	my @HitAminoEnds;
	my @HitNuclStarts;
	my @HitNuclEnds;
	my @HitEVals;
	my $num_blat_hits = 0;
	my $blatf = OpenInputFile($blat_out_fname);
	while (my $line = <$blatf>) {
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
		    } else {
			if ($nucl_start >= $prev_start && $nucl_start <= $prev_end) {
			    $was_prev_used = 1;
			    last;
			}
			if ($nucl_end >= $prev_start && $nucl_end <= $prev_end) {
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
		$num_blat_hits++;

	    }
	}
	close($blatf);

	# For outputting, let's get the textual representation of direction
	# into the chromosome name
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

	my $target_info = "MSA Ali Region : $exon_str\n";
	$target_info = $target_info."    Target Genome  : $target_species ($chr:$SearchRanges[0]..$SearchRanges[1])\n";
	$target_info = $target_info."    Source Species : $source_species\n";
	
	# Is it an especially elusive ghost we're chasing?
	if ($num_blat_hits == 0) {
	    print $outf "[ ] Search failure (no BLAT hits)\n";
	    print $outf "    $target_info";
	    print $outf "    Search Sequence: $SearchSeqs[$q]\n\n";
	    next;
	}
	
	# Oh, this is a most profitable Ghost Adventure indeed!
	$ghosts_busted++;

	# Let's illustrate how much of the sequence has been covered.
	# NOTE: We're assuming that our hits are consistent with one another,
	#   but we may want to double-check that in the future...
	my @MappedSeq = split(//,lc($SearchSeqs[$q]));
	for (my $hit=0; $hit<$num_blat_hits; $hit++) {
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
	print $outf "    Search Sequence: $mapped_seq\n";
	print $outf "    Num BLAT Hits  : $num_blat_hits\n";

	# I'm going to take this 'underlining' out for now, and let the
	# upper / lower case distinction speak for itself.
	#
	#print $outf "    ";
	#foreach my $char (@MappedSeq) {
	#if ($char eq uc($char)) { print $outf '-'; }
	#else                    { print $outf ' '; }
	#}
	#print $outf "\n";

	for (my $hit=0; $hit<$num_blat_hits; $hit++) {
	    print $outf "    + Aminos $HitAminoStarts[$hit]..$HitAminoEnds[$hit] ";
	    print $outf "mapped to $target_species $chr:$HitNuclStarts[$hit]..$HitNuclEnds[$hit] ";
	    print $outf "($HitEVals[$hit])\n";
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





# EOF

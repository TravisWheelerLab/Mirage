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
sub GetMappedSeqMSA;
sub ParseAFA;
sub RecordSplicedMSA;



##############
#            #
#   SCRIPT   #
#            #
##############



if (@ARGV != 2) { die "\n  USAGE:  ./Oasis.pl [Mirage-Results] [Species-Guide]\n\n"; }

# Figure out what the location of the Mirage src directory is
my $location = $0;
$location =~ s/Oasis\.pl$//;

# We're going to need these friends
my $sindex = $location.'../inc/hsi/sindex';
my $sfetch = $location.'../inc/hsi/sfetch';
my $sstat  = $location.'../inc/hsi/sstat';


# Confirm that the input directory looks like the real deal
my $input_dirname = ConfirmDirectory($ARGV[0]);
my $final_results_dirname = ConfirmDirectory($input_dirname.'Final-MSAs');
my $all_species_dirname = ConfirmDirectory($input_dirname.'Species-MSAs');


# TODO: Parse other arguments
my $write_spliced_msas = 1; # Do we want to write our spliced MSAs to files?


# Start off by parsing the species guide, which will give
# us genome locations and chromosome lengths.
my $SpeciesGuide = OpenInputFile($ARGV[1]);
my %SpeciesToGenomes;
my %ChrLensBySpecies;
while (my $line = <$SpeciesGuide>) {

    $line =~ s/\n|\r//g;
    next if ($line !~ /(\S+)\s+(\S+)\s+\S+/);
    my $species = lc($1);
    my $genome  = $2;

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
if ($write_spliced_msas) {
    if (-d $spliced_dirname) {
	print "\n";
	print "  Warning:  Directory of MSAs with marked splice sites ($spliced_dirname) located.\n";
	print "            MSAs may be over-written.\n\n";
    } else {
	CreateDirectory($spliced_dirname);
    }
}


# Cool!  Now that we've imported the species guide, we can get down
# to the real business!
#
# We'll do this the easiest (laziest) way, which will just be running
# through each of our final MSAs, loading in *only* those sequences
# which have mappings, and then reverse-engineering exon markers.
#
# Fun!
#
my $FinalMSAs = OpenDirectory($final_results_dirname);
while (my $fname = readdir($FinalMSAs)) {

    next if ($fname !~ /(\S+)\.afa$/);
    my $gene = $1;
    $fname = $final_results_dirname.$fname;

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
    if ($write_spliced_msas) {
	RecordSplicedMSA($spliced_dirname.$gene.'.afa',\@MSA,\@SeqNames,$num_seqs,$msa_len);
    }
    
}
closedir($FinalMSAs);

1;




###################
#                 #
#   SUBROUTINES   #
#                 #
###################





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
	$seqname =~ /\>([^\|]+)\|/;
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
	my $mapf = OpenInputfile($mapfname);
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
    foreach my $species (keys %SpeciesToMapfiles) {
	if ($SpeciesToMapfiles{$species} eq "-1") {
	    $SpeciesToMapfiles{$species} = 0;
	}
    }
    return (0,0,0,0,0,0) if (scalar(keys %SpeciesToMapfiles) <= 1);
    

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
		$MSA[$num_seqs][$j] = $MSA[$i][$j];
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





# EOF

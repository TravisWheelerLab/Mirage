#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;
use Cwd;
use Getopt::Long;
use Time::HiRes;

# Now tell me, who are you?
sub GetThisDir { my $lib = $0; $lib =~ s/\/Bazaar.pl$//; return $lib; }
use lib GetThisDir();
use BureaucracyMirage;
use DisplayProgress;


# SUBROUTINES
sub ParseSeqName;
sub OppositeStrands;



#############
#           #
#  SCRIPTS  #
#           #
#############



if (@ARGV != 2) { die "\n  USAGE:  ./Bazaar.pl [Mirage-Results] [Species-Guide]\n\n"; }


# What's life without your best friend (sfetch)?
my $sfetch = $0;
$sfetch =~ s/Bazaar\.pl/\.\.\/inc\/hsi\/sfetch/;
if (!(-e $sfetch)) {
    die "\n  ERROR:  Failed to location sfetch (expected to be at '$sfetch')\n\n";
}

# If there isn't a folder holding our misses, then there's nothing
# bizarre about this output.
my $miragedirname = ConfirmDirectory($ARGV[0]);
my $missdirname = $miragedirname.'Mapping-Misses/';
if (!(-d $missdirname)) {
    die "\n  No chimeras detected among sequences in '$miragedirname'\n\n";
}

# Cool!  Misses!  Which ones are chimeras?
my $chimeragrep = OpenSystemCommand("grep 'Chimeric Mapping' $missdirname\*");
my @ChimericSeqNames;
while (my $line = <$chimeragrep>) {
    if ($line =~ /\.misses\:([^\#]+)/) {
	my $seqname = $1;
	$seqname =~ s/\s+$//;
	push(@ChimericSeqNames,$seqname);
    }
}
close($chimeragrep);

# Last chance to get away with not doing work today...
if (scalar(@ChimericSeqNames) == 0) {
    die "\n  No chimeras detected among sequences in '$miragedirname'\n\n";
}


# CHIMERAS AHOY!  Next up, we'll need to parse our species guide
# to learn where all our happy little genomes live
my %SpeciesToGenome;
my $SpeciesGuide = OpenInputFile($ARGV[1]);
while (my $line = <$SpeciesGuide>) {
    if ($line =~ /^(\S+)\s+(\S+)/) {
	my $species = lc($1);
	my $genome  = $2;
	if (!(-e $genome)) {
	    print "\n  Warning: Can't find $species genome (looking at '$genome')\n\n";
	} else {
	    $SpeciesToGenome{$species} = $genome;
	}
    }
}
close($SpeciesGuide);


# How foolish it would be to not create an output directory
my $outdirname = CreateDirectory('Bazaar-Results');


# Now we'll run through our chimeras learning their secrets!
# NOTE: Just in case the user updated their species guide and has
#   lost some species we'll need to double-check genome existence...
my $current_species = '';
my $current_genome = 0;
my $current_species_dirname = '';

my $current_num_chimeras;
my $current_num_badorder;
my $current_num_multichr;
my $current_num_revexons;
my $current_num_long_introns;
my @SpeciesChimeras;
my @BadOrderSeqs;
my @MultiChrSeqs;
my @RevExonSeqs;
my @LongIntronSeqs;

for (my $seq_id=0; $seq_id<scalar(@ChimericSeqNames); $seq_id++) {

    my $seqname = $ChimericSeqNames[$seq_id];
    my ($species,$gene,$iso_id) = ParseSeqName($seqname);

    if ($species ne $current_species) {

	# Write species summary out to file
	if ($current_genome && $current_num_chimeras) {
	    my $outf = OpenOutputFile($current_species_dirname.'SUMMARY.out');
	    print $outf "Species            : $species\n";
	    print $outf "Num. Chimeric Seq.s: $current_num_chimeras\n";
	    print $outf "\n";
	    print $outf "All Chimeric Sequences\n";
	    print $outf "----------------------\n";
	    foreach my $chimeric_seq (@SpeciesChimeras) {
		print $outf "  $chimeric_seq\n";
	    }
	    print $outf "\n";

	    if ($current_num_badorder) {
		print $outf "Sequences with Biologically Inconsistent Splicing\n";
		print $outf "-------------------------------------------------\n";
		foreach my $chimeric_seq (@BadOrderSeqs) {
		    print $outf "  $chimeric_seq\n";
		}
		print $outf "\n";
	    }
	    
	    if ($current_num_multichr) {
		print $outf "Sequences Mapped to Multiple Chromosomes\n";
		print $outf "----------------------------------------\n";
		foreach my $chimeric_seq (@MultiChrSeqs) {
		    print $outf "  $chimeric_seq\n";
		}
		print $outf "\n";
	    }
	    
	    if ($current_num_revexons) {
		print $outf "Sequences with Exons that Flip Strand Direction\n";
		print $outf "-----------------------------------------------\n";
		foreach my $chimeric_seq (@RevExonSeqs) {
		    print $outf "  $chimeric_seq\n";
		}
		print $outf "\n";
	    }
	    
	    if ($current_num_long_introns) {
		print $outf "Sequences with Questionably Long Introns\n";
		print $outf "----------------------------------------\n";
		foreach my $chimeric_seq (@LongIntronSeqs) {
		    print $outf "  $chimeric_seq\n";
		}
		print $outf "\n";
	    }

	    print $outf "\n";
	    close($outf);

	}

	$current_species = $species;

	$current_num_chimeras = 0;
	$current_num_badorder = 0;
	$current_num_multichr = 0;
	$current_num_revexons = 0;
	$current_num_long_introns = 0;
	@SpeciesChimeras = ();
	@BadOrderSeqs = ();
	@MultiChrSeqs = ();
	@RevExonSeqs = ();
	@LongIntronSeqs = ();

	if ($SpeciesToGenome{$current_species}) {

	    $current_genome = $SpeciesToGenome{$current_species};
	    $current_species_dirname
		= ConfirmDirectory($miragedirname.'Species-MSAs/'.$current_species);

	} else {

	    print "\n  Warning: Provided 'Species-Guide' file differs from that used by Mirage\n";
	    print "           (species '$species' has mappings, but no listed genome)\n\n";
	    $current_genome = 0;

	}
    }

    next if (!$current_genome);

    
    # Alrighty, it looks like we've got a genome -- let's pull in the
    # protein's amino acid sequence.    
    my $msafname = $current_species_dirname.'alignments/'.$gene.'.afa';
    my $msaf = OpenInputFile($msafname);
    while (my $line = <$msaf>) {
	$line =~ s/\n|\r//g;
	last if ($line =~ /\>$seqname$/);
    }

    # UH-OH!
    if (eof($msaf)) {
	close($msaf);
	die "\n  ERROR:  Failed to find sequence '$seqname' in file '$msafname'\n\n";
    }

    my $seqstr = '';
    while (my $line = <$msaf>) {
	last if ($line =~ /\>/);
	$line =~ s/\n|\r|\s|\*//g;
	$seqstr = $seqstr.uc($line);
    }
    close($msaf);

    if (!$seqstr) {
	die "\n  ERROR:  Sequence '$seqname' is empty in file '$msafname'?!\n\n";
    }

    my @Aminos = split(//,$seqstr);

    
    # Woohoo!  That's what I call a sequence!
    # Next up, pull in the mapping (both the ranges 'n' the codon centers)

    
    my $mapfname = $current_species_dirname.'mappings/'.$gene.'.afa';
    my $mapf = OpenInputFile($mapfname);
    my $canon_chr = <$mapf>;
    $canon_chr =~ s/^Canonical Chromosome\:\s+//;
    while (my $line = <$mapf>) {
	$line =~ s/\n|\r//g;
	if ($line =~ /Sequence ID\: $seqname$/) {
	    $line = <$mapf>; # Eat mapping method
	    $line = <$mapf>; # Eat chromosome(s)
	    last;
	}
    }

    # 'Nother chance to catch a wuh-woh
    if (eof($mapf)) {
	close($mapf);
	die "\n  ERROR:  Failed to find mapping for '$seqname' in '$mapfname'\n\n";
    }

    my $num_exons = <$mapf>;
    $num_exons =~ s/Num Exons  \: (\d+)//;
    $num_exons = $1;

    my @ExonAminoStarts;
    my @ExonAminoEnds;
    my @ExonChrs;
    my @ExonChrStarts;
    my @ExonChrEnds;
    my @ExonCoordListStrs;
    for (my $i=0; $i<$num_exons; $i++) {

	# Summary info
	my $line = <$mapf>;
	$line =~ /Aminos (\d+)\.\.(\d+)\, ([^\:]+)\:(\d+)\.\.(\d+)/;
	push(@ExonAminoStarts,$1);
	push(@ExonAminoEnds,$2);
	push(@ExonChrs,$3);
	push(@ExonChrStarts,$4);
	push(@ExonChrEnds,$5);

	# Coords
	$line = <$mapf>;
	$line =~ s/\n|\r//g;
	push(@ExonCoordListStrs,$line);
	
    }

    close($mapf);


    my $canon_revcomp = 0;
    $canon_revcomp = 1 if ($canon_chr =~ /\[revcomp\]/);

    # Now that we've read in our mappings, we can look to see what's making
    # this sequence look so cuh-razy!
    # We'll group the exons together into internally-coherent groups.
    my $saw_badorder = 0;
    my $saw_multichr = 0;
    my $saw_revexon = 0;
    my $saw_long_intron = 0;
    my @ExonGroupStarts;
    my @ExonGroupEnds;
    my @ExonGroupChrs;
    my @ExonGroupInfo;
    my $num_exon_groups = -1;
    for (my $i=0; $i<$num_exons; $i++) {

	# Are we changing chromosomes / strands? (or looking at the first exon?)
	if ($ExonChrs[$i-1] ne $ExonChrs[$i] || $i==0) {

	    $num_exon_groups++;
	    $ExonGroupStarts[$num_exon_groups] = $i;
	    $ExonGroupEnds[$num_exon_groups] = $i;
	    $ExonGroupChrs[$num_exon_groups] = $ExonChrs[$i];

	    if ($ExonChrs[$i] eq $canon_chr) {
		$ExonGroupInfo[$num_exon_groups] = 'std';
	    } elsif (OppositeStrands($canon_chr,$ExonChrs[$i])) {
		$ExonGroupInfo[$num_exon_groups] = 'revchr';
		$saw_revexon = 1;
	    } else {
		$ExonGroupInfo[$num_exon_groups] = 'noncanon';
		$saw_multichr = 1;
	    }

	} elsif (abs($ExonChrEnds[$i-1]-$ExonChrStarts[$i]) > 2000000) {

	    my $intron_len = abs($ExonChrEnds[$i-1]-$ExonChrStarts[$i]);

	    # Long intron!  Not as long as Mirage needed to call an intron
	    # as long, but might still be worth pointing out...
	    $num_exon_groups++;
	    $ExonGroupStarts[$num_exon_groups] = $i;
	    $ExonGroupEnds[$num_exon_groups] = $i;
	    $ExonGroupChrs[$num_exon_groups] = $ExonChrs[$i];
	    $ExonGroupInfo[$num_exon_groups] = 'long-intron:'.$intron_len;
	    $saw_long_intron = 1;
	    
	} elsif (($ExonChrs[$i] =~ /\[revcomp\]/ && $ExonChrStarts[$i] > $ExonChrStarts[$i-1])
		 || ($ExonChrs[$i] !~ /\[revcomp\]/ && $ExonChrStarts[$i] < $ExonChrStarts[$i-1])) {
	    
	    # We're continuing along the same chromsome, but we're pulling our
	    # nucls from a "biologically inconsistent" region
	    $num_exon_groups++;
	    $ExonGroupStarts[$num_exon_groups] = $i;
	    $ExonGroupEnds[$num_exon_groups] = $i;
	    $ExonGroupChrs[$num_exon_groups] = $ExonChrs[$i];
	    $ExonGroupInfo[$num_exon_groups] = 'badorder';
	    $saw_badorder = 1;
	    
	} else {
	    
	    # Maybe we're weird, but we're weird together!
	    $ExonGroupEnds[$num_exon_groups] = $i;
	    
	}
	
    }

    # We'll be one short in our count, so fix that
    $num_exon_groups++;

    # Record the weirdnesses that we saw for our summary statistics.
    $current_num_chimeras++;
    push(@SpeciesChimeras,$seqname);
    
    if ($saw_badorder) {
	$current_num_badorder++;
	push(@BadOrderSeqs,$seqname);
    }
    if ($saw_multichr) {
	$current_num_multichr++;
	push(@MultiChrSeqs,$seqname);
    }
    if ($saw_revexon) {
	$current_num_revexons++;
	push(@RevExonSeqs,$seqname);
    }
    if ($saw_long_intron) {
	$current_num_long_introns++;
	push(@LongIntronSeqs,$seqname);
    }

    # Now we'll pull in the DNA sequence for each of our exons
    my @ExonAminoAliStrs;
    my @ExonNuclAliStrs;
    my $big_nucl_str;
    my $aa_pos = 0;
    for (my $exon_id=0; $exon_id<$num_exons; $exon_id++) {

	my $search_chr = $ExonChrs[$exon_id];
	my $strand_dir = 1;
	if ($search_chr =~ /\[revcomp\]/) {
	    $strand_dir = -1;
	    $search_chr =~ s/\[revcomp\]//;
	}
	
	my $sfetchcmd = $sfetch.' -range '.$ExonChrStarts[$exon_id].'..'.$ExonChrEnds[$exon_id];
	$sfetchcmd = $sfetchcmd.' '.$current_genome.' '.$search_chr;

	my $sfetchf = OpenSystemCommand($sfetchcmd);
	my $line = <$sfetchf>; # Eat the file name
	my $nuclseq = '';
	while ($line = <$sfetchf>) {
	    $line =~ s/\n|\r//g;
	    $nuclseq = $nuclseq.uc($line) if ($line);
	}
	close($sfetchf);
	
	# NOICE!  Them's some dope-phresche nucleotides.  Time to make a
	# good ol' fashioned alignment outta that madness.
	my @Nucls = split(//,$nuclseq);
	my @CodonCenters = split(/\,/,$ExonCoordListStrs[$exon_id]);
	my @ExonAminoAli;
	my @ExonNuclAli;
	my $alipos = 0;
	my $nuclpos = 0;
	my $chrpos = $ExonChrStarts[$exon_id];
	my $codon_id = 0;
	while ($alipos < scalar(@Nucls)) {

	    $ExonNuclAli[$alipos] = $Nucls[$nuclpos++];
	    $chrpos += $strand_dir;

	    if ($codon_id < scalar(@CodonCenters) && $CodonCenters[$codon_id] == $chrpos) {

		$ExonAminoAli[$alipos++] = $Aminos[$aa_pos++];
		$codon_id++;

		# If we're inserting (relative to genome), then cruise through the
		# insertion region.
		if ($codon_id < scalar(@CodonCenters) && $CodonCenters[$codon_id] == $CodonCenters[$codon_id-1]) {

		    # Pull in the next nucleotide, to round out the last codon
		    $ExonAminoAli[$alipos]  = ' ';
		    $ExonNuclAli[$alipos++] = $Nucls[$nuclpos++];
		    $chrpos += $strand_dir;

		    while ($codon_id < scalar(@CodonCenters) && $CodonCenters[$codon_id] == $CodonCenters[$codon_id-1]) {

			$ExonAminoAli[$alipos]  = ' ';
			$ExonNuclAli[$alipos++] = '-';
			
			$ExonAminoAli[$alipos]  = $Aminos[$aa_pos++];
			$ExonNuclAli[$alipos++] = '-';
			$codon_id++;

			$ExonAminoAli[$alipos]  = ' ';
			$ExonNuclAli[$alipos++] = '-';

		    }

		}
		
	    } else {

		$ExonAminoAli[$alipos] = ' ';
		$alipos++;

	    }
	    
	}

	# PHEW, that was some good aligning!
	#
	# Before we can get around to producing a translated alignment,
	# if there are any frameshift-y things happening we'll need to
	# add some gaps into the nucleotide sequence.
	#
	my $amino_ali_str = '';
	my $nucl_ali_str = '';
	for (my $i=0; $i<$alipos; $i++) {

	    $amino_ali_str = $amino_ali_str.$ExonAminoAli[$i];
	    $nucl_ali_str  = $nucl_ali_str.$ExonNuclAli[$i];

	    if ($ExonAminoAli[$i] =~ /[A-Z]/) {
		
		# Check if either of the next two positions are also translated
		if ($i+1 < $alipos && $ExonAminoAli[$i+1] =~ /[A-Z]/) {
		    $amino_ali_str = $amino_ali_str.'  ';
		    $nucl_ali_str  = $nucl_ali_str.'--';
		} elsif ($i+2 < $alipos && $ExonAminoAli[$i+2] =~ /[A-Z]/) {
		    $amino_ali_str = $amino_ali_str.'  ';
		    $nucl_ali_str  = $nucl_ali_str.$ExonNuclAli[++$i].'-';
		}

	    }
	    
	}

	# Aligned!
	$ExonAminoAliStrs[$exon_id] = $amino_ali_str;
	$ExonNuclAliStrs[$exon_id]  = $nucl_ali_str;
	$big_nucl_str = $big_nucl_str.$nucl_ali_str;
	
    }

    # Finally, do some of that excellent translating that you've become famous for!
    my @AllNucls = split(//,$big_nucl_str);
    my @ExonTransAliStrs;
    my $all_nucls_pos = 0;
    for (my $exon_id=0; $exon_id<$num_exons; $exon_id++) {

	my $trans_ali_str = '';
	foreach my $char (split(//,$ExonAminoAliStrs[$exon_id])) {

	    if ($char =~ /[A-Z]/) {

		my $codon = TranslateCodon($AllNucls[$all_nucls_pos-1].$AllNucls[$all_nucls_pos].$AllNucls[$all_nucls_pos+1]);
		$trans_ali_str = $trans_ali_str.$codon;

	    } else {
		
		$trans_ali_str = $trans_ali_str.' ';

	    }
	    
	    $all_nucls_pos++;

	}

	$ExonTransAliStrs[$exon_id] = $trans_ali_str;
	
    }


    #
    #  --> ALIGNMENTS COMPLETE! 
    #
    #  Time to write this sucker out to a file!
    #

    my $outf = OpenOutputFile($current_species_dirname.$gene.'.'.$iso_id.'.out');

    print $outf "Sequence Name  : $seqname\n";
    print $outf "Num Exon Groups: $num_exon_groups\n";
    for (my $i=0; $i<$num_exon_groups; $i++) {

	my $start_exon = $ExonGroupStarts[$i]+1;
	my $end_exon   = $ExonGroupEnds[$i]+1;
	my $group_info = $ExonGroupInfo[$i];
	my $group_id   = $i+1;

	print $outf "  + Exon Group $i: Exons $start_exon\.\.$end_exon\n";
	if ($group_info eq 'std') {
	    
	    if ($start_exon == $end_exon) {
		print $outf "    Canonical exon for this gene family\n";
	    } else {
		print $outf "    Canonical exons for this gene family\n";
	    }
	    
	} elsif ($group_info eq 'revchr') {
	    
	    print $outf "    Encoded on the reverse complement strand (relative to canon exons)\n";
	    
	} elsif ($group_info eq 'noncanon') {

	    print $outf "    Encoded on a non-canonical chromosome ($ExonGroupChrs[$i])\n";

	} elsif ($group_info =~ 'long-intron:(\d+)') {

	    my $intron_len = $1;
	    print $outf "    Same chromosome as previous exon, but long intervening intron ($intron_len nucleotides)\n";

	} elsif ($group_info eq 'badorder') {

	    print $outf "    Placed in a biologically inconsistent location relative to previous exon\n";

	} else {
	    print "\n  OOPS!   Alex forgot about exons grouped as '$group_info'\n\n";
	}
    }
    print $outf "\n\n";

    # Now we'll provide the alignments for each group (along with a reminder of the
    # reason for our grouping
    for (my $i=0; $i<$num_exon_groups; $i++) {

	my $start_exon = $ExonGroupStarts[$i];
	my $end_exon   = $ExonGroupEnds[$i];
	my $num_exons  = ($end_exon - $start_exon) + 1;
	my $group_info = $ExonGroupInfo[$i];
	my $group_id   = $i+1;

	print $outf "\n";
	print $outf "--------------------------------------------------------------\n";
	print $outf "\n\n";
	print $outf "  Group ID    : $group_id\n";
	print $outf "  Group Class : ";
	if    ($group_info eq 'std'        ) { print $outf "Canonical";          }
	elsif ($group_info eq 'revchar'    ) { print $outf "Strand Flip";        }
	elsif ($group_info eq 'noncanon'   ) { print $outf "Alt. Chromosome";    }
	elsif ($group_info =~ 'long-intron') { print $outf "Long Intron";        }
	elsif ($group_info eq 'badorder'   ) { print $outf "Inconsistent Order"; }
	print $outf "\n";
	print $outf "  Chromosome  : $ExonGroupChrs[$group_id]\n";
	print $outf "  Num Exons   : $num_exons\n";
	for (my $exon_id=$start_exon; $exon_id<=$end_exon; $exon_id++) {
	    my $output_id = $exon_id+1;
	    print $outf "  + Exon $output_id: ";
	    print $outf "Aminos $ExonAminoStarts[$exon_id]..$ExonAminoEnds[$exon_id],";
	    print $outf " Nucls $ExonChrStarts[$exon_id]..$ExonChrEnds[$exon_id]\n";
	}
	print $outf "\n";

	my $group_amino_ali = $ExonAminoAliStrs[$start_exon];
	my $group_nucl_ali  = $ExonNuclAliStrs[$start_exon];
	my $group_trans_ali = $ExonTransAliStrs[$start_exon];
	for (my $exon_id=$start_exon+1; $exon_id<=$end_exon; $exon_id++) {
	    $group_amino_ali = $group_amino_ali.'|'.$ExonAminoAliStrs[$exon_id];
	    $group_nucl_ali  =  $group_nucl_ali.'|'.$ExonNuclAliStrs[$exon_id];
	    $group_trans_ali = $group_trans_ali.'|'.$ExonTransAliStrs[$exon_id];
	}
	my @GroupAminos = split(//,$group_amino_ali);
	my @GroupNucls  = split(//,$group_nucl_ali);
	my @GroupTrans  = split(//,$group_trans_ali);

	my $line_len = 60;
	for (my $line_start = 0; $line_start < scalar(@GroupNucls); $line_start += $line_len) {
	    print $outf "\n";

	    print $outf "  Aminos: ";
	    for (my $line_pos = $line_start; $line_pos < Min($line_start+$line_len,scalar(@GroupNucls)); $line_pos++) {
		print $outf "$GroupAminos[$line_pos]";
	    }
	    print $outf "\n";

	    print $outf "  Nucl.s: ";
	    for (my $line_pos = $line_start; $line_pos < Min($line_start+$line_len,scalar(@GroupNucls)); $line_pos++) {
		print $outf "$GroupNucls[$line_pos]";
	    }
	    print $outf "\n";

	    print $outf "  Trans.: ";
	    for (my $line_pos = $line_start; $line_pos < Min($line_start+$line_len,scalar(@GroupNucls)); $line_pos++) {
		print $outf "$GroupTrans[$line_pos]";
	    }
	    print $outf "\n";

	    print $outf "\n";
	}

	print $outf "\n\n";
	
    }
    
    close($outf);
    
}


1;




#################
#               #
#  SUBROUTINES  #
#               #
#################





##########################################################################
#
#  Function: ParseSeqName
#
sub ParseSeqName
{
    my $seqname = shift;

    my $species;
    my $gene;
    my $iso_id;
    if ($seqname =~ /^([^\|]+)\|([^\|]+)\|([^\|]+)\s*$/) {

	$species = lc($1);
	$iso_id = $3;
	my @GeneList = split(/\//,$2);
	$gene = $GeneList[0];

    } else {

	$seqname =~ /OS\=([^\=]+\=?)/;
	$species = lc($1);

	$seqname =~ /GN\=([^\=]+\=?)/;
	$gene = lc($1);

	# We'll use the accession as our isoform id
	$seqname =~ /^[^\|]+\|([^\|]+)\|/;
	$iso_id = $1;

    }

    if (!($iso_id && $gene && $species)) {
	die "\n  ERROR:  Failed to parse name of sequence '$seqname'\n\n";
    }

    return($species,$gene,$iso_id);
    
}






##########################################################################
#
#  Function: OppositeStrands
#
sub OppositeStrands
{
    my $chr1 = shift;
    my $chr2 = shift;
    if ($chr1.'[revcomp]' eq $chr2 || $chr2.'[revcomp]' eq $chr1) {
	return 1;
    }
    return 0;
}












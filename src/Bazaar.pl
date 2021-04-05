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
sub OppositeStrands;



#############
#           #
#  SCRIPTS  #
#           #
#############



if (@ARGV != 2) { die "\n  USAGE:  ./Bazaar.pl [Mirage-Results] [Species-Guide]\n\n"; }


# What's life without your best friend (sfetch)?
my $sfetch = $0;
my $sfetch =~ s/Bazaar\.pl/\.\.\/inc\/hsi\/sfetch/;
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
    if ($line =~ /\.misses\:(\S+)/) {
	push(@ChimericSeqNames,$1);
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

my $current_num_chimera;
my $current_num_badorder;
my $current_num_multichr;
my $current_num_revexons;
my $current_num_long_introns;
my @BadOrderSeqs;
my @MultiChrSeqs;
my @RevExonSeqs;
my @LongIntronSeqs;

for (my $seq_id=0; $seq_id<scalar(@ChimericSeqNames); $seq_id++) {

    my $seqname = $ChimericSeqNames[$seq_id];
    $seqname =~ /^([^\|]+)\|/;
    my $species = lc($1);

    if ($species ne $current_species) {

	# Write species summary out to file
	if ($current_genome) {
	    my $species_outf = OpenOutputFile($current_species_dirname.'SUMMARY.out');
	    # TO DO
	    close($species_outf);
	}

	$current_species = $species;

	$current_num_chimera = 0;
	$current_num_badorder = 0;
	$current_num_multichr = 0;
	$current_num_revexons = 0;
	$current_num_long_introns = 0;
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

    
    # Get the (primary) gene for this sequence and locate its final MSA
    $seqname =~ /^[^\|]+\|([^\|]+)\|/;
    my $gene_list_str = $1;
    my @GeneList = split(/\//,$gene_list_str);
    my $gene = lc($GeneList[0]);

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

    $seqstr = '';
    while (my $line = <$msaf>) {
	last if ($line =~ /\>/);
	$line =~ s/\n|\r|\s|\*//g;
	$seqstr = $seqstr.uc($line);
    }
    close($msaf);

    if (!$seqstr) {
	die "\n  ERROR:  Sequence '$seqname' is empty in file '$msafname'?!\n\n";
    }

    my @Seq = split(//,$seqstr);

    
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


    # Now that we've read in our mappings, we can look to see what's making
    # this sequence look so cuh-razy!
    # We'll group the exons together into internally-coherent groups.
    my @ExonGroupStarts;
    my @ExonGroupEnds;
    my @ExonGroupInfo;
    my $num_exon_groups = -1;
    for (my $i=0; $i<$num_exons; $i++) {

	# Are we changing chromosomes / strands? (or looking at the first exon?)
	if ($ExonChrs[$i-1] ne $ExonChrs[$i] || $i==0) {

	    $num_exon_groups++;
	    $ExonGroupStarts[$num_exon_groups] = $i;
	    $ExonGroupEnds[$num_exon_groups] = $i;

	    if ($ExonChrs[$i] eq $canon_chr) {
		$ExonGroupInfo[$num_exon_groups] = 'std';
	    } elsif (OppositeStrands($canon_chr,$ExonChrs[$i])) {
		$ExonGroupInfo[$num_exon_groups] = 'revchr';
	    } else {
		$ExonGroupInfo[$num_exon_groups] = 'noncanon';
	    }

	} elsif (abs($ExonChrEnds[$i-1]-$ExonChrStarts[$i]) > 2000000) {

	    # Long intron!  Not as long as Mirage needed to call an intron
	    # as long, but might still be worth pointing out...
	    $num_exon_groups++;
	    $ExonGroupStarts[$num_exon_groups] = $i;
	    $ExonGroupEnds[$num_exon_groups] = $i;
	    $ExonGroupInfo[$num_exon_groups] = 'long-intron';
	    
	} else {

	    # Maybe we're weird, but we're weird together!
	    $ExonGroupEnds[$num_exon_groups] = $i;
	    
	}
	
    }

    # We'll be one short in our count, so fix that
    $num_exon_groups++;


    # Now we'll pull in the DNA sequence for each of our exons
    my @ExonNuclSeq;
    for (my $i=0; $i<$num_seqs; $i++) {
	
    }
    
}


1;




#################
#               #
#  SUBROUTINES  #
#               #
#################




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












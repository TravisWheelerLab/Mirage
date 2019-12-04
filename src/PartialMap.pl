#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub ReadMSA;
sub BuildMapMSA;
sub SquishMSAs;


# NOTES (to keep me from losing my mind)
#
#   * This isn't going to be the 'whole-shebang' for identifying partial maps.  This script will focus on
#       identifying segments that we can confidently call partial maps without using BLAT.
#
#   * If at least one sequence mapped, we can identify where in the MSA there's a quality change between the
#       mapping sequences and each non-mapping sequence.  We can then split the non-mapping sequence at that point
#       and see if we can get one or both halves to map to the region, noting which positions remain unmapped.
#
#   * If no sequences mapped, BUT we have some amount of consensus on a region to test based on near-hits,
#       then we 
#
#   * If no sequences mapped, AND we have no idea where to look.... Hmmmm.....
#       My inclination is to look for a high-identity span and then re-BLAT with a sequence made from just that
#       region, but we'd need to unify all such sequences to avoid getting bogged down in a billion hours of BLAT.
#


if (@ARGV < 4 || (@ARGV - 4) % 3 != 0) { 
    die "\n  USAGE: PartialMap.pl [fam.afa] [fam.hits.out] [ProteinDB.fa] [Genome.fa] { [seq] [chr] [range] }\n\n"; 
}

my $location = $0;
$location =~ s/PartialMap\.pl//;

# Get those essentials on hand
my $spaln      = $location.'../inc/spaln2.3.3/src/spaln';
my $eslsfetch  = $location.'../inc/easel/miniapps/esl-sfetch';
my $eslseqstat = $location.'../inc/easel/miniapps/esl-seqstat';

my $inmsafile = $ARGV[0];
my $inhitfile = $ARGV[1];
my $protdb    = $ARGV[2];
my $genome    = $ARGV[3];

# We're (hopefully) going to need to know the lengths of the chromosomes
# in our genome, so let's get our hands on that file Quilter generated.
open(my $chrlenfile,'<',$genome.'.chr_len_index');
my %ChrLens;
while (my $line = <$chrlenfile>) {
    $line =~ s/\n|\r//g;
    if ($line =~ /^(\S+) (\d+)$/) {
	$ChrLens{$1} = $2;
    }
}
close($chrlenfile);

my @NearHitSeqs;
my @NearHitChrs;
my @NearHitRanges;
my $i = 4;
while ($i<scalar(@ARGV)) {
    push(@NearHitSeqs,$ARGV[$i++]);
    push(@NearHitChrs,$ARGV[$i++]);
    push(@NearHitRanges,$ARGV[$i++]);
}

my ($msaref,$seqnamesref,$num_seqs,$msa_len) = ReadMSA($inmsafile);
my @MSA = @{$msaref};
my @SeqNames = @{$seqnamesref};

# In case we've recorded any ARF info, we'll need the original names
my @OrigNames;
for ($i=0; $i<$num_seqs; $i++) {
    my $OrigNames[$i] = $SeqNames[$i];
    $OrigNames[$i] =~ s/\|ARFs?\:[^\|]+\|/\|/;
}

# Use the hitfile to identify where specific positions in the MSA 
# can be attributable to mapping coordinates
my ($mapmsaref,$usesmapref,$num_maps,$map_chr,$min_genome_index,$max_genome_index) 
    = BuildMapMSA(\@MSA,\@OrigNames,$num_seqs,$msa_len,$inhitfile);
my @MapMSA = @{$mapmsaref};
my @UsesMap = @{$usesmapref};

# Alrighty, lads.  Now we have the info. that we need to determine
# which of the paths we're on.
if ($num_maps) {

    # RAD! Our MSA is already built around there being some way to map
    # to a particular range of a chromosome, so we'll work off of that.
    
    # First step is to squish our MSAs into consensus-like arrays (really, more like profiles)
    my @MappedSeqIDs;
    my @UnmappedSeqIDs;
    for (my $i=0; $i<$num_seqs; $i++) {
	if ($UsesMap[$i]) { push(@MappedSeqIDs,$i);   }
	else              { push(@UnmappedSeqIDs,$i); }
    }

    # The format of the MSA profile is a comma-separated list of character:count pairs,
    # and the format of the mapping profile is the same, but with positions instead of characters.
    # If indels (for non-mappers) created a column with no mapped sequences, it'll have a '0'
    my ($msasquish_ref,$mapmsasquish_ref) = SquishMSAs(\@MSA,\@MapMSA,\@MappedSeqIDs,$num_seqs,$num_maps,$msa_len);
    my @MSASquish = @{$msasquish_ref};
    my @MapMSASquish = @{$mapmsasquish_ref};

    # Because we're working with a specific sequence range, we can extract the specific region of the
    # genome that we're restricted to considering (we'll pull in an extra 100Kb for good measure, though)
    my $chr = $map_chr;
    $chr =~ s/\[revcomp\]//;

    my $low_genome_pos  = $min_genome_index - 100000;
    $low_genome_pos     = 1 if ($low_genome_pos < 1);
    my $high_genome_pos = $max_genome_index + 100000;
    $high_genome_pos    = $ChrLens{$chr} if ($high_genome_pos > $ChrLens{$chr});

    # NOTE: We'll go ahead and flip things if we're in revcomp land
    if ($map_chr =~ /\[revcomp\]/) {
	my $temp_genome_pos = $low_genome_pos;
	$low_genome_pos = $high_genome_pos;
	$high_genome_pos = $temp_genome_pos;
    }

    # We should be all set for extraction!
    my $dnafilename = $inmsafile;
    $dnafilename =~ s/\.afa/\.dna\.fa/;
    my $sfetchCmd = $eslsfetch." -c $low_genome_pos\.\.$high_genome_pos \"$genome\" \"$chr\" > \"$dnafilename\"";
    if (system($sfetchCmd)) { die "\n  Failed to pull in DNA region ($sfetchCmd)\n\n"; }

    # Let's also get a handle on where the exon start and stop positions are.
    my @ExonStarts;
    my @ExonEnds;
    my $num_exons = 0;
    push(@ExonStarts,1);
    for (my $i=2; $i<$msa_len-1; $i++) {
	if ($MSASquish[$i] eq '*') {
	    $ExonEnds[$num_exons] = $i-1;
	    $num_exons++;
	    $ExonStarts[$num_exons] = $i+1;
	    $i+=2;
	}
    }
    $ExonEnds[$num_exons] = $msa_len-1;
    $num_exons++;

    # For each of the non-mapping sequences, compute percent ID they have with each of the exons that
    # they have some amount of sequence in.
    my @ExonPctsID;
    my @SplitPos;
    foreach (my $j=0; $j<scalar(@UnmappedSeqIDs); $j++) {

	my $seqid = $UnmappedSeqIDs[$j];
	my @PctList;
	my @Content; # How many characters did this sequence align to this exon? (matches + mismatches)

	for (my $exon=0; $exon<$num_exons; $exon++) {

	    # We'll call something a match if it matches any character represented in our profile (of mapped seqs)
	    my $num_matches    = 0;
	    my $num_mismatches = 0;

	    for (my $i=$ExonStarts[$exon]; $i<=$ExonEnds[$exon]; $i++) {
		if ($MSA[$seqid][$i] =~ /[A-Za-z]/) {
		    if ($MSASquish[$i]) {
			my $capChar = uc($MSA[$seqid][$i]);
			if ($MSASquish[$i] =~ /$capChar\:\d/) {
			    $num_matches++;
			} else {
			    $num_mismatches++;
			}
		    } else {
			$num_mismatches++; # Inserted characters count as mismatches to the exon
		    }
		}
	    }

	    # Is there a percent identity to think about here?
	    if ($num_matches || $num_mismatches) {
		my $num_chars = $num_mismatches + $num_matches + 0.0;
		my $pctID = ($num_matches + 0.0) / $num_chars;
		$ExonPctsID[$j][$exon] = $pctID;
		push(@PctList,$pctID);
		push(@Content,$num_chars);
	    } else {
		$ExonPctsID[$j][$exon] = '-'; # This exon is unused by this sequence.
	    }

	}

	# Now that we have the percents-identity we can find if there's a particular part of the
	# sequence that's interfering with our ability to map to the genome beautifully.
	#
	# We allow a sequence to be as low-identity as 80% without being "bad"
	#
	my @GoodOrBad;
	foreach my $pct (@PctList) {
	    if ($pct >= 0.8) { push(@GoodOrBad,'g'); }
	    else             { push(@GoodOrBad,'b'); }
	}

	# We'll simplify by figuring out the lengths of the runs of good/bad alignments
	my @StateRuns;
	my $current_state = $GoodOrBad[0];
	my $state_run_len = 1;
	my $state_chars   = $Content[0];
	for (my $i=0; $i<scalar(@GoodOrBad); $i++) {
	    if ($GoodOrBad[$i] eq $current_state) {
		$state_run_len++;
		$state_chars += $Content[$i];
	    } else {
		push(@StateRuns,$current_state.':'.$state_run_len.'/'.$state_chars);
		$current_state = $GoodOrBad[$i];
		$state_run_len = 1;
		$state_chars   = $Content[$i];
	    }
	}
	push(@StateRuns,$current_state.':'.$state_run_len.'/'.$state_chars);


	# Once we have our state runs, we'll need to make some executive decisions based around how much
	# variation there is in the quality of aligned regions.
	if (scalar(@StateRuns) > 1) {
	    
	    # NOTE: We'll handle the single state run (homogenously good/bad) differently in several ways,
	    #  hence these somewhat goofus embedded conditionals.
	    
	    # So, now our task is to see if there's a clear division between a good part and a bad part...
	    my $last_left_exon = -1; # This will be w.r.t. exons that this sequence aligned to
	    
	    if (scalar(@StateRuns) == 2) {
		
		$StateRuns[0] =~ /\:(\d+)\//;
		$last_left_exon = $1;
		
		
	    } elsif (scalar(@StateRuns) > 2) {
		
		
		# What we'll need to do is figure out which way of splitting the sequence into two parts
		# does the best job of segregating good parts to one side and bad parts to the other.
		#
		#        L   R
		#      .---.---.
		#    b |   |   |
		#      |---+---|
		#    g |   |   |
		#      '---'---'
		#
		# Maximize  abs(Lb - Lg) + abs(Rb - Rg)
		#
		
		# Initialize the matrix
		my @CountMatrix;
		$CountMatrix[0][0] = 0;
		$CountMatrix[0][1] = 0;
		$CountMatrix[1][0] = 0;
		$CountMatrix[1][1] = 0;
		
		# Start off with the last left exon being the furthest left
		$StateRuns[0] =~ /(\S)\:\d+\/(\d+)$/;
		$current_state = $1;
		$state_chars = $2;
		if ($current_state eq 'b') { $CountMatrix[0][1] = $state_chars; }
		else                       { $CountMatrix[0][0] = $state_chars; }
		
		# Run along filling out the rest of the matrix (assuming all other exons belong to the right half)
		for (my $i=1; $i<scalar(@StateRuns); $i++) {
		    $StateRuns[$i] =~ /(\S)\:\d+\/(\d+)$/;
		    $current_state = $1;
		    $state_chars = $2;
		    if ($current_state eq 'b') { $CountMatrix[1][1] += $state_chars; }
		    else                       { $CountMatrix[1][0] += $state_chars; }
		}
		
		# Figure out how this placement treats us
		my $last_left_state = 0;
		my $top_val = abs($CountMatrix[0][0]-$CountMatrix[0][1]) + abs($CountMatrix[1][1]-$CountMatrix[1][0]);
		
		# Run along, moving characters to the left side from the right
		for (my $i=1; $i<scalar(@StateRuns)-1; $i++) {
		    
		    $StateRuns[$i] =~ /(\S)\:\d+\/(\d+)$/;
		    $current_state = $1;
		    $state_chars = $2;
		    if ($current_state eq 'b') {
			$CountMatrix[0][1] += $state_chars;
			$CountMatrix[1][1] -= $state_chars;
		    } else {
			$CountMatrix[0][0] += $state_chars;
			$CountMatrix[1][0] -= $state_chars;
		    }
		    
		    # Does this break point work better than our current fav?
		    my $new_val = abs($CountMatrix[0][0]-$CountMatrix[0][1]) + abs($CountMatrix[1][1]-$CountMatrix[1][0]);
		    if ($new_val > $top_val) {
			$last_left_state = $i;
			$top_val = $new_val;
		    }
		    
		}
		
		# We now know which of the entries in 'StateRuns' tells us where to split things,
		# but we'll still need to do a tiny bit of work to figure out the actual number of
		# exons inwards that we have to count.
		$last_left_exon = 0;
		for (my $i=0; $i<=$last_left_state; $i++) {
		    $StateRuns[$i] =~ /\:(\d+)\//;
		    $state_run_len = $1;
		    $last_left_exon += $state_run_len;
		}
		
		
	    }

	    # Find where in the MSA corresponds to the identified exon
	    my $break_pos = -1;
	    for (my $i=0; $i<$last_left_exon; $i++) {
		# Advance to the next exon with content for this sequence
		for (my $j=$break_pos+1; $j<$num_exons; $j++) {
		    if ($ExonPctsID[$j][$exon] ne '-') {
			$break_pos = $j;
			last;
		    }
		}
	    }
	    
	    # At this point, 'break_pos' is the index of the exon that marks the
	    # end of the left-half of the sequence (vis-a-vis quality break).
	    # Translate that into the position of the dividing '*' in the MSA.
	    $break_pos = $ExonEnds[$break_pos] + 1;
	    
	    
	    #
	    # The next task will be to rip this sequence apart and see if we can
	    # map one half (or both halves?!) back to the genome
	    #
	    
	    
	    $SeqNames[$seqid] =~ /^[^\|]+\|([^\|]+)\|/;
	    my $iso = $1;
	    
	    my $Lprotfname = $infmsaname;
	    $Lprotfname =~ s/\.afa$/\./;
	    $Lprotfname = $Lprotfname.$iso.'.left-prot.fa';
	    
	    my $Rprotfname = $infmsaname;
	    $Rprotfname =~ s/\.afa$/\./;
	    $Rprotfname = $Rprotfname.$iso.'.right-prot.fa';
	    
	    # Run along the sequence until we hit the breakpoint, constructing the
	    # sequence that we'll write out to our file.
	    #
	    # More importantly, why the HECK did I make 'j' the sequence index?
	    #
	    my $Lstr = '>Left-'.$gene.'-'.$iso."\n";
	    my $newline = 60;
	    for (my $i=0; $i<$break_pos; $i++) {
		if ($MSA[$seqid][$i] =~ /[A-Za-z]/) {
		    $Lstr = $Lstr.$MSA[$seqid][$i];
		    $newline--;
		    if (!$newline) {
			$Lstr = $Lstr."\n";
			$newline = 60;
		    }
		}
	    }
	    $Lstr = $Lstr."\n" if ($newline != 60);
	    
	    open(my $Lprotf,'>',$Lprotfname);
	    print $Lprotf "$Lstr\n";
	    close($Lprotf);
	    
	    
	    my $Rstr = '>Right-'.$gene.'-'.$iso."\n";
	    $newline = 60;
	    for (my $i=$break_pos; $i<$msa_len; $i++) {
		if ($MSA[$seqid][$i] =~ /[A-Za-z]/) {
		    $Rstr = $Rstr.$MSA[$seqid][$i];
		    $newline--;
		    if (!$newline) {
			$Rstr = $Rstr."\n";
			$newline = 60;
		    }
		}
	    }
	    $Rstr = $Rstr."\n" if ($newline != 60);
	    
	    open(my $Rprotf,'>',$Rprotfname);
	    print $Rprotf "$Rstr\n";
	    close($Rprotf);
	    
	    
	    # OK, well what are you waiting for? Do you seriously thing Spaln is just going
	    # to run itself? You dork, you gotta call those commands, what is this, your first
	    # day of programming?
	    my $Lspalnfname = $Lprotfname;
	    $Lspalnfname =~ s/\-prot\.fa/\-spaln\.out/;
	    
	    my $Rspalnfname = $Rprotfname;
	    $Rspalnfname =~ s/\-prot\.fa/\-spaln\.out/;
	    
	    my $spalnCmd = $spaln.' -Q3 -O1 -S3 -ya3 "'.$dnafilename.'" "[prot]" 1>"[out]" 2>/dev/null';
	    
	    my $LspalnCmd = $spalnCmd;
	    $LspalnCmd =~ s/\[prot\]/$Lprotfname/;
	    $LspalnCmd =~ s/\[out\]/$Lspalnfname/;
	    if (system($LspalnCmd)) { die "\n  Spaln command '$LspalnCmd' failed\n\n"; }
	    
	    my $RspalnCmd = $spalnCmd;
	    $RspalnCmd =~ s/\[prot\]/$Rprotfname/;
	    $RspalnCmd =~ s/\[out\]/$Rspalnfname/;
	    if (system($RspalnCmd)) { die "\n  Spaln command '$RspalnCmd' failed\n\n"; }

	    # Now that we have our two separate SPALNs, let's see if they give us *ANY* decent mapping.
	    # This is largely just a copy of the now infamous 'ParseSPALNOutput' from Quilter, but with
	    # some of the finer points shaved down, 'cuz that's the way I do it, baby.
	    
	    
	} else {

	    #
	    # This is where we find ourselves when there's only one state (good or bad).
	    # In this case we just re-run SPALN and examine whether there's just a little
	    # goof preventing us from calling it a good mapping...
	    #
	    
	}

    }
	
    # Alrighty, lads -- we find ourselves faced with 2 options.
    #
    # 1. Do a full-protein Spaln to see if we're in a situation where there's an *okay* hit for the
    #    full length, but with some minor issues that (now that we have a scaffold of sorts) we're
    #    willing to overlook.
    #
    # 2. Try to find a difference between the high-quality parts of the alignment and low-quality parts,
    #    and split up our protein for separate Spalns of the parts.
    #
    # ---> OBSERVATION: If we're in a case like (1), shouldn't (2) still work?
    #                   It isn't elegant, but presumably the result should be the same, right?

} elsif (@ARGV > 4) {

    # If we don't have a context to work in, mayhaps we can draw consensus on a chromosome
    # (and, based on that consensus, a range on that chromosome to work in).
    my %ChrCounts;
    for (my $i=0; $i<scalar(@NearHitChrs); $i++) {
	my $next_chr = $NearHitChrs[$i];
	if ($ChrCounts{$next_chr}) { $ChrCounts{$next_chr}++; }
	else                       { $ChrCounts{$next_chr}=1; }
    }

    my $top_chr = '';
    my $top_chr_counts = 0;
    foreach my $candidate_chr (keys %ChrCounts) {
	if ($ChrCounts{$candidate_chr} > $top_chr_counts) {
	    $top_chr = $candidate_chr;
	    $top_chr_counts = $ChrCounts{$candidate_chr};
	}
    }

    # Is there some semblance of consensus around the particular chromosome we're likely to find
    # success with? (requiring at least 50% - one over for odd numbers of near-hitters).
    if ($top_chr_counts >= int((scalar(@NearHitChrs)/2.0))+0.5) {

	# HELLA SWAG!
	my $chr = $top_chr;

	# Now, where on EARTH are we looking on that chromosome?!
	my $low_genome_pos  = 10 ** 10;
	my $high_genome_pos = 0;

	for (my $i=0; $i<scalar(@NearHitChrs); $i++) {
	    if ($NearHitChrs[$i] eq $top_chr) {
		$NearHitRanges[$i] =~ /(\d+)\.\.(\d+)/;
		my $range_start = $1;
		my $range_end   = $2;
		# I'm not sure if 'start' is always lower, so let's be chill for like one second
		if ($range_start < $low_genome_pos ) { $low_genome_pos  = $range_start; }
		if ($range_start > $high_genome_pos) { $high_genome_pos = $range_start; }
		if ($range_end   < $low_genome_pos ) { $low_genome_pos  = $range_end;   }
		if ($range_end   > $high_genome_pos) { $high_genome_pos = $range_end;   }
	    }
	}

	# Let's pull in that range, with a spare 1Mb on either side
	$low_genome_pos  -= 1000000;
	$low_genome_pos   = 1 if ($low_genom_pos < 1);
	$high_genome_pos += 1000000;
	$high_genome_pos  = $ChrLens{$chr} if ($high_genome_pos > $ChrLens{$chr});

	my $dnafilename = $inmsafile;
	$dnafilename =~ s/\.afa/\.dna\.fa/;
	my $sfetchCmd = $eslsfetch." -c $low_genome_pos\.\.$high_genome_pos \"$genome\" \"$chr\" > \"$dnafilename\"";
	if (system($sfetchCmd)) { die "\n  Failed to pull in DNA region ($sfetchCmd)\n\n"; }

    }

} else {
    # No "else" since we require some bit of info. better than just redoing Quilter's work
    # ... At least not yet.  We might try to ID a happy chunk and restrict ourselves just to
    #     partial sequences.
}


1;






##########################################################################
#
#  FUNCTION:  ReadMSA
#
sub ReadMSA
{
    my $msaname = shift;
    open(my $msafile,'<',$msaname);

    my @MSA;
    my @SeqNames;
    my $num_seqs = -1;
    my $msa_len  = 0;

    while (my $line = <$msafile>) {

	$line =~ s/\n|\r//g;
	if ($line =~ /\>(\S+)/) {
	    my $seqname = $1;
	    push(@SeqNames,$seqname);
	    $num_seqs++;
	    $msa_len = 0;
	} else {
	    foreach my $char (split(//,$line)) {
		$MSA[$num_seqs][$msa_len] = $char;
		$msa_len++;
	    }
	}

    }

    close($msafile);

    return(\@MSA,\@SeqNames,$num_seqs,$msa_len);

}







##########################################################################
#
#  FUNCTION:  BuildMapMSA
#
sub BuildMapMSA
{
    my $msaref = shift;
    my $seqnamesref = shift;
    my $num_seqs = shift;
    my $msa_len = shift;
    my $hitfilename = shift;

    my @MSA = @{$msaref};
    my @SeqNames = @{$seqnamesref}; # (secretly 'OrigNames' -- no snitchin' now!)

    # Prime the MSA
    my @MapMSA;
    my @UsesMap;
    for (my $i=0; $i<$num_seqs; $i++) {
	$UsesMap[$i] = 0;
	for (my $j=0; $j<$msa_len; $j++) {
	    $MapMSA[$i][$j] = 0;
	}
    }

    my $num_maps = 0;
    my $map_chr;

    my $min_genome_index = 10 ** 10;
    my $max_genome_index = 0;

    open(my $hitfile,'<',$hitfilename);
    while (my $line = <$hitfile>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	if ($line =~ /Isoform ID \: (\S+)/) {

	    my $seqname = $1;

	    my $i = 0;
	    $i++ while ($SeqNames[$i] ne $seqname);
	   
	    $line = <$hitfile>; # Nobody ever cares about the method line :'(

	    $line = <$hitfile>;
	    $line =~ /Chromosome \: (\S+)/;
	    my $chr = $1;

	    $map_chr = $chr;
	    
	    $line = <$hitfile>;
	    $line =~ /Match Pos\.s\: (\S+)/;
	    my $hit_coord_str = $1;
	    $hit_coord_str =~ s/\,\*//g;
	    my @HitCoords = split(/\,/,$hit_coord_str);

	    my $range_low  = $HitCoords[0];
	    my $range_high = $HitCoords[scalar(@HitCoords)-1];

	    if ($range_low  < $min_genome_index) { $min_genome_index = $range_low;  }
	    if ($range_high > $max_genome_index) { $max_genome_index = $range_high; }

	    $UsesMap[$i] = $chr.':'.$range_low.'..'.$range_high;
	    $num_maps++;
 
	    my $z = 0;
	    for (my $j=0; $j<$msa_len; $j++) {
		if ($MSA[$i][$j] =~ /[A-Za-z]/) {
		    $MapMSA[$i][$j] = $HitCoords[$z];
		    $z++;
		}
	    }
	    
	} elsif ( $line =~ /\- UNUSED MAPPINGS \-/) {
	    last;
	}

    }
    close($hitfile);

    return (\@MapMSA,\@UsesMap,$num_maps,$map_chr,$min_genome_index,$max_genome_index);

}







##########################################################################
#
#  FUNCTION:  SquishMSAs
#
sub SquishMSAs
{

    my $msa_ref = shift;
    my $mapmsa_ref = shift;
    my $mappedseqids_ref = shift;

    my $num_seqs = shift;
    my $num_maps = shift;
    my $msa_len  = shift;
    
    my @MSA = @{$msa_ref};
    my @MapMSA = @{$mapmsa_ref};
    my @MappedSeqIDs = @{$mappedseqids_ref};

    my @MSASquish;
    my @MapMSASquish;

    for (my $j=0; $j<$msa_len; $j++) {

	my %CharCounts;
	my %MapPosCounts;

	my $is_splice_site = 0;

	# Special case: splice site
	if ($MSA[$MappedSeqIDs[0]][$j] eq '*') {
	    $MSASquish[$j]    = '*';
	    $MapMSASquish[$j] = '*'; # I don't think this will ever matter...
	    next;
	}

	
	for (my $i=0; $i<$num_maps; $i++) {
	    
	    my $seq_id = $MappedSeqIDs[$i];

	    my $char = uc($MSA[$seq_id][$j]);
	    if ($char =~ /[A-Z]/) {
		if ($CharCounts{$char}) { $CharCounts{$char}++; }
		else                    { $CharCounts{$char}=1; }
	    }

	    my $map_pos = $MapMSA[$seq_id][$j];
	    if ($map_pos) {
		if ($MapPosCounts{$map_pos}) { $MapPosCounts{$map_pos}++; }
		else                         { $MapPosCounts{$map_pos}=1; }
	    }

	}

	# For now we won't worry about ranking these -- maybe add in later?

	# Characters
	my @CharList = keys %CharCounts;
	if (@CharList > 0) {
	    foreach my $char (@CharList) {
		my $entry = $char.':'.$CharCounts{$char};
		if ($MSASquish[$j]) { $MSASquish[$j] = $MSASquish[$j].','.$entry; }
		else                { $MSASquish[$j] = $entry;                    }
	    }
	} else {
	    $MSASquish[$j] = 0;
	}

	# Positions
	my @PosList = keys %MapPosCounts;
	if (@PosList > 0) {
	    foreach my $pos (@PosList) {
		my $entry = $pos.':'.$MapPosCounts{$pos};
		if ($MapMSASquish[$j]) { $MapMSASquish[$j] = $MapMSASquish[$j].','.$entry; }
		else                   { $MapMSASquish[$j] = $entry;                       }
	    }
	} else {
	    $MapMSASquish[$j] = 0;
	}

    }

    return (\@MSASquish,\@MapMSASquish);

}







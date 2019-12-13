#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


sub ReadMSA;
sub BuildMapMSA;
sub SquishMSAs;
sub PMParseSPALNOutput;
sub SelenocysteineCheck;
sub UpdateMapMSA;


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


# NOTE: We should be perfectly okay taking '-' if there isn't an actual 'fam.hits.out'
if (@ARGV < 4 || (@ARGV - 4) % 3 != 0) { 
    die "\n  USAGE: PartialMap.pl [fam.afa] [fam.hits.out] [ProteinDB.fa] [Genome.fa] { [seq] [chr] [range] }\n\n"; 
}

my $location = $0;
$location =~ s/PartialMap\.pl//;

# Get those essentials on hand
my $spaln      = $location.'../inc/spaln2.3.3/src/spaln';
my $eslsfetch  = $location.'../inc/easel/miniapps/esl-sfetch';
my $eslseqstat = $location.'../inc/easel/miniapps/esl-seqstat';

my $inmsafname = $ARGV[0];
my $inhitfname = $ARGV[1];
my $protdb     = $ARGV[2];
my $genome     = $ARGV[3];

$inmsafname =~ /\/?([^\/]+)\.afa$/;
my $gene = $1;

# We're (hopefully) going to need to know the lengths of the chromosomes
# in our genome, so let's get our hands on that file Quilter generated.
open(my $chrlenfile,'<',$genome.'.chr_length_index') || die "\n  ERROR:  Failed to locate chromosome length index file for '$genome'\n\n";
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

my ($msaref,$seqnamesref,$num_seqs,$msa_len) = ReadMSA($inmsafname);
my @MSA = @{$msaref};
my @SeqNames = @{$seqnamesref};

# In case we've recorded any ARF info, we'll need the original names
my @OrigNames;
for ($i=0; $i<$num_seqs; $i++) {
    $OrigNames[$i] = $SeqNames[$i];
    $OrigNames[$i] =~ s/\|ARFs?\:[^\|]+\|/\|/;
}

# Use the hitfile to identify where specific positions in the MSA 
# can be attributable to mapping coordinates
my ($mapmsaref,$usesmapref,$num_maps,$map_chr,$min_genome_index,$max_genome_index)
    = BuildMapMSA(\@MSA,\@OrigNames,$num_seqs,$msa_len,$inhitfname);
my @MapMSA = @{$mapmsaref};
my @UsesMap = @{$usesmapref};

# We'll just make ourselves a couple of cozy little arrays
my @MappedSeqIDs;
my @UnmappedSeqIDs;
for (my $i=0; $i<$num_seqs; $i++) {
    if ($UsesMap[$i]) { push(@MappedSeqIDs,$i);   }
    else              { push(@UnmappedSeqIDs,$i); }
}

# If we ever find ourselves latching onto a special chromosome, we'll
# deinitely want to know who the lucky lady / boylady is!
my $chr;

# Alrighty, lads.  Now we have the info. that we need to determine
# which of the paths we're on.
#
# Our overarching goal through this section is to fill in MapMSA
# for any sequences that didn't start off with a mapping.
#
my @NewlyMappedSeqs;
my @MSASquish;
my @MapMSASquish;
if ($num_maps) {

    # RAD! Our MSA is already built around there being some way to map
    # to a particular range of a chromosome, so we'll work off of that.
    
    # The format of the MSA profile is a comma-separated list of character:count pairs,
    # and the format of the mapping profile is the same, but with positions instead of characters.
    # If indels (for non-mappers) created a column with no mapped sequences, it'll have a '0'
    my ($msasquish_ref,$mapmsasquish_ref) = SquishMSAs(\@MSA,\@MapMSA,\@MappedSeqIDs,$num_seqs,$num_maps,$msa_len);
    @MSASquish = @{$msasquish_ref};
    @MapMSASquish = @{$mapmsasquish_ref};

    # Because we're working with a specific sequence range, we can extract the specific region of the
    # genome that we're restricted to considering (we'll pull in an extra 100Kb for good measure, though)
    $chr = $map_chr;
    $chr =~ s/\[revcomp\]//;

    my $low_genome_pos  = $min_genome_index - 100000;
    $low_genome_pos     = 1 if ($low_genome_pos < 1);
    my $high_genome_pos = $max_genome_index + 100000;
    $high_genome_pos    = $ChrLens{$chr} if ($high_genome_pos > $ChrLens{$chr});

    # NOTE: We'll go ahead and flip things if we're in revcomp land
    my $revcomp = 0;
    if ($map_chr =~ /\[revcomp\]/) {
	$revcomp = 1;
	my $temp_genome_pos = $low_genome_pos;
	$low_genome_pos = $high_genome_pos;
	$high_genome_pos = $temp_genome_pos;
    }

    # We should be all set for extraction!
    my $dnafilename = $inmsafname;
    $dnafilename =~ s/\.afa/\.dna\.fa/;
    my $sfetchCmd = $eslsfetch." -c $low_genome_pos\.\.$high_genome_pos \"$genome\" \"$chr\" > \"$dnafilename\"";
    if (system($sfetchCmd)) { die "\n  Failed to pull in DNA region ($sfetchCmd)\n\n"; }

    # This is somewhat unnecessary, but I want to keep my verbage
    # as close as possible to the original ParseSPALNOutput
    my $offset = $low_genome_pos-1; # NOTE: I'm not 100% on why we do this subtraction, but *whatever*

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
	    for (my $i=0; $i<$last_left_exon; $i++) { # Just counting up
		# Advance to the next exon with content for this sequence
		for (my $exon=$break_pos+1; $exon<$num_exons; $exon++) {
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
	    
	    my $Lprotfname = $inmsafname;
	    $Lprotfname =~ s/\.afa$/\./;
	    $Lprotfname = $Lprotfname.$iso.'.left-prot.fa';
	    
	    my $Rprotfname = $inmsafname;
	    $Rprotfname =~ s/\.afa$/\./;
	    $Rprotfname = $Rprotfname.$iso.'.right-prot.fa';
	    
	    # Run along the sequence until we hit the breakpoint, constructing the
	    # sequence that we'll write out to our file.
	    #
	    # More importantly, why the HECK did I make 'j' the sequence index?
	    #
	    my $Lstr = '>Left-'.$gene.'-'.$iso."\n";
	    my $newline = 60;
	    my $Roffset = 0; # This will count what the index of the first Right seq char is
	    for (my $i=0; $i<$break_pos; $i++) {
		if ($MSA[$seqid][$i] =~ /[A-Za-z]/) {
		    $Roffset++;
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
	    
	    # Now we can run our two separate SPALNs and see if they give us *ANY* decent mapping.
	    # This is largely just a copy of the now infamous 'ParseSPALNOutput' from Quilter, but with
	    # some of the finer points shaved down, 'cuz that's the way I do it, baby.
	    # The returned hashes map amino indices from the split sequences to positions in the genome.
	    
	    my $LspalnCmd = $spalnCmd;
	    $LspalnCmd =~ s/\[prot\]/$Lprotfname/;
	    $LspalnCmd =~ s/\[out\]/$Lspalnfname/;
	    my $LAITGPref = PMParseSPALNOutput($LspalnCmd,$Lspalnfname,$revcomp,$offset,$Lprotfname);
	    my %LAminoIndexToGenPos = %{$LAITGPref};
	    
	    my $RspalnCmd = $spalnCmd;
	    $RspalnCmd =~ s/\[prot\]/$Rprotfname/;
	    $RspalnCmd =~ s/\[out\]/$Rspalnfname/;
	    my $RAITGPref = PMParseSPALNOutput($RspalnCmd,$Rspalnfname,$revcomp,$offset,$Rprotfname);
	    my %RAminoIndexToGenPos = %{$RAITGPref};


	    # First -- are these mappings consistent (assuming we did get successful mappings from both)?
	    my @LeftAminoKeys   = sort {$a <=> $b} keys %LAminoIndexToGenPos;
	    my @RightAminoKeys  = sort {$a <=> $b} keys %RAminoIndexToGenPos;
	    next if (scalar(@LeftAminoKeys) + scalar(@RightAminoKeys) == 0);
	    
	    if (scalar(@LeftAminoKeys) && scalar(@RightAminoKeys)) {

		my $consistent = 1; # Let's just assume the best?
		if ($revcomp) {
		    $consistent = 0 if ($LeftAminoKeys[scalar(@LeftAminoKeys)-1] < $RightAminoKeys[0]);
		} else {
		    $consistent = 0 if ($LeftAminoKeys[scalar(@LeftAminoKeys)-1] > $RightAminoKeys[0]);
		}
		
		# If the mappings aren't consistent, I guess we'll pick the longest one?
		if (!$consistent) {
		    if (scalar(@LeftAminoKeys) > scalar(@RightAminoKeys)) {
			%RAminoIndexToGenPos = ();
			@RightAminoKeys = ();
		    } else {
			%LAminoIndexToGenPos = ();
			@LeftAminoKeys = ();			
		    }
		}
	    }

	    # Because our 'Right' indices are all w.r.t. the right part of the sequence
	    # we're going to have to shift them.  These shifts will increase the values,
	    # so we start with the highest and decrease.  We'll remove entries to avoid
	    # conflicts when we integrate w/ LAITGP
	    for (my $i=scalar(@RightAminoKeys)-1; $i>=0; $i--) {
		my $key = $RightAminoKeys[$i];
		$RAminoIndexToGenPos{$key+$Roffset} = $RAminoIndexToGenPos{$key};
		$RAminoIndexToGenPos{$key} = 0;
	    }

	    # Now we combine 'em.  Yeah.  Think about that.
	    my %AminoIndexToGenPos;
	    foreach my $aminoindex (keys %LAminoIndexToGenPos) { $AminoIndexToGenPos{$aminoindex} = $LAminoIndexToGenPos{$aminoindex}; }
	    foreach my $aminoindex (keys %RAminoIndexToGenPos) { $AminoIndexToGenPos{$aminoindex} = $RAminoIndexToGenPos{$aminoindex}; }

	    # UPDATE OUR MAPMSA
	    $mapmsaref = UpdateMapMSA(\@MSA,\@MapMSA,\%AminoIndexToGenPos,$seqid,$msa_len);
	    @MapMSA    = @{$mapmsaref};
	    push(@NewlyMappedSeqs,$seqid);
	    
	} else {

	    #
	    # This is where we find ourselves when there's only one state (good or bad).
	    # In this case we just re-run SPALN and examine whether there's just a little
	    # goof preventing us from calling it a good mapping...
	    #

	    $SeqNames[$seqid] =~ /^[^\|]+\|([^\|]+)\|/;
	    my $iso = $2;

	    my $protfname = $inmsafname;
	    $protfname =~ s/\.afa$/\./;
	    $protfname = $protfname.$iso.'.prot.fa';

	    my $Pstr = '>'.$SeqNames[$seqid]."\n";
	    my $newline = 60;
	    for (my $i=0; $i<$msa_len; $i++) {
		if ($MSA[$seqid][$i] =~ /[A-Za-z]/) {
		    $Pstr = $Pstr.$MSA[$seqid][$i];
		    $newline--;
		    if (!$newline) {
			$Pstr = $Pstr."\n";
			$newline = 60;
		    }
		}
	    }
	    $Pstr = $Pstr."\n" if ($newline != 60);

	    open(my $Protf,'>',$protfname);
	    print $Protf "$Pstr\n";
	    close($Protf);

	    my $spalnfname = $protfname;
	    $spalnfname =~ s/prot\.fa/spaln\.out/;

	    my $spalnCmd = $spaln.' -Q3 -O1 -S3 -ya3 "'.$dnafilename.'" "'.$protfname.'" 1>"'.$spalnfname.'" 2>/dev/null';
	    my $AITGPref = PMParseSPALNOutput($spalnCmd,$spalnfname,$revcomp,$offset,$protfname);
	    my %AminoIndexToGenPos;

	    # UPDATE OUR MAPMSA
	    $mapmsaref = UpdateMapMSA(\@MSA,\@MapMSA,\%AminoIndexToGenPos,$seqid,$msa_len);
	    @MapMSA    = @{$mapmsaref};
	    push(@NewlyMappedSeqs,$seqid);

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
    #
    # Also, we won't do anything to remove [revcomp]s -- DEAL WITH IT
    #
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
	$chr = $top_chr;

	# Now we can remove [revcomp] from chr if appropriate (for pulling purposes)
	my $revcomp = 0;
	if ($chr =~ /\[revcomp\]/) {
	    $revcomp = 1;
	    $chr =~ s/\[revcomp\]//;
	}

	# Now, where on EARTH are we looking on that chromosome?!
	my $low_genome_pos  = 10 ** 10;
	my $high_genome_pos = 0;

	for (my $i=0; $i<scalar(@NearHitChrs); $i++) {
	    if ($NearHitChrs[$i] eq $top_chr) {

		$NearHitRanges[$i] =~ /(\d+)\.\.(\d+)/;
		my $range_start = $1;
		my $range_end   = $2;
		
		# If we're in revcomp land the 'start' index will be larger than the 'end' index
		if ($range_start < $low_genome_pos ) { $low_genome_pos  = $range_start; }
		if ($range_start > $high_genome_pos) { $high_genome_pos = $range_start; }
		if ($range_end   < $low_genome_pos ) { $low_genome_pos  = $range_end;   }
		if ($range_end   > $high_genome_pos) { $high_genome_pos = $range_end;   }

	    }
	}

	# Let's pull in that range, with a spare 1Mb on either side
	$low_genome_pos  -= 1000000;
	$low_genome_pos   = 1 if ($low_genome_pos < 1);
	$high_genome_pos += 1000000;
	$high_genome_pos  = $ChrLens{$chr} if ($high_genome_pos > $ChrLens{$chr});

	# Again, we're going to do this thing I just love doing and will never stop doing EVER
	my $offset = $low_genome_pos-1;

	my $dnafilename = $inmsafname;
	$dnafilename =~ s/\.afa/\.dna\.fa/;
	my $sfetchCmd = $eslsfetch." -c $low_genome_pos\.\.$high_genome_pos \"$genome\" \"$chr\" > \"$dnafilename\"";
	if (system($sfetchCmd)) { die "\n  Failed to pull in DNA region ($sfetchCmd)\n\n"; }


	# As is our natural desire, we shall examine with great care the relationship of
	# our protein sequences to this region of the genome.
	for (my $i=0; $i<$num_seqs; $i++) {
	    
	    # NOTE that this portion is basically identical to how we handle single-state
	    # sequences when we have an established genome region to examine, but that after
	    # we get the mapping things take a different turn.
	    $SeqNames[$i] =~ /^[^\|]+\|([^\|]+)\|/;
	    my $iso = $1;

	    my $protfname = $inmsafname;
	    $protfname =~ s/\.afa$/\./;
	    $protfname = $protfname.$iso.'.prot.fa';
	    
	    my $Pstr = '>'.$SeqNames[$i]."\n";
	    my $newline = 60;
	    for (my $j=0; $j<$msa_len; $j++) {
		if ($MSA[$i][$j] =~ /[A-Za-z]/) {
		    $Pstr = $Pstr.$MSA[$i][$j];
		    $newline--;
		    if (!$newline) {
			$Pstr = $Pstr."\n";
			$newline = 60;
		    }
		}
	    }
	    $Pstr = $Pstr."\n" if ($newline != 60);

	    open(my $Protf,'>',$protfname);
	    print $Protf "$Pstr\n";
	    close($Protf);

	    my $spalnfname = $protfname;
	    $spalnfname =~ s/prot\.fa/spaln\.out/;

	    my $spalnCmd = $spaln.' -Q3 -O1 -S3 -ya3 "'.$dnafilename.'" "'.$protfname.'" 1>"'.$spalnfname.'" 2>/dev/null';
	    my $AITGPref = PMParseSPALNOutput($spalnCmd,$spalnfname,$revcomp,$offset,$protfname);
	    my %AminoIndexToGenPos = %{$AITGPref};

	    next if (!scalar(keys %AminoIndexToGenPos));
	    
	    # UPDATE OUR MAPMSA
	    $mapmsaref = UpdateMapMSA(\@MSA,\@MapMSA,\%AminoIndexToGenPos,$i,$msa_len);
	    @MapMSA    = @{$mapmsaref};
	    push(@NewlyMappedSeqs,$i);

	}

	# Before we call things quits-ies, should we consider marking splice-sites?
	#
	# I THINK NOT! If we only have partially-mapped sequences, then that could
	# misleadingly indicate that the unmapped regions are single exons that
	# failed to map.


    }

} else {

    # No "else" since we require some bit of info. better than just redoing Quilter's work.
    #
    # NOTE: While I sort of hinted at leaving the door open to this, I had forgotten that we
    #       now hold onto the BLAT indices that we've played around with, so there really
    #       doesn't seem to be any reason to think we can try something different with these
    #       sequences and suddenly find ourselves with a decent map.
    #
    exit(0);

}


# If we didn't get anything out of updating our mappings, there's no more work
# to do -- jump ship.
exit(0) if (scalar(@NewlyMappedSeqs) == 0);



#########################################################################################
#########################################################################################
###                                                                                   ###
###  BIG DEAL!!! We now have (1) an updated MapMSA and (2) a list of NewlyMappedSeqs  ###
###                                                                                   ###
#########################################################################################
#########################################################################################


# Sadly, our sense of national unity was short-lived.  We'll need to split off 
# again based on whether or not any sequences mapped initially, since we'll either
# be looking for confirmation of mapped characters (do their NW-based positions
# correspond to characters placed by mapping?) or else... something?

if ($num_maps) {

    # We'll write out the partial mappings to the family's hits.out file.
    # '-' will be used to indicate any unmapped positions.
    my $partialhitstr = "------------------------ PARTIAL MAPPINGS -----------------------\n\n";

    

    # For each of the sequences that we were able to (at least partially) map back to
    # the genome, check whether their newfound mapping positions make sense with the
    # positions established for the other sequences.
    foreach my $seq (@NewlyMappedSeqs) {

	# We'll track all of the positions that are consistent with the established
	# mappings and all of the inconsistent positions, because we can.
	#
	# NOTE that we'll make the values in the arrays [1..seqlen], since that's
	# the way we'll want to record that information for our bio-pals.
	#
	my @ConsistentPos;
	my @InconsistentPos;
	my $consistencies = 0;
	my $inconsistencies = 0;

	# ALSO, since I guess I never got this count, for the sake of measuring
	# success, let's attend to how much of each newly mapped sequence was included
	# in the new mapping.
	my $seqlen = 0;

	for (my $j=0; $j<$msa_len; $j++) {

	    # If this position has a letter, I'll bet it's part of the sequence.
	    # You can't pull the wool over my eyes!
	    $seqlen++ if ($MSA[$seq][$j] =~ /[A-Za-z]/);

	    # If this position has been imbued with a mapping... check it out!
	    if ($MapMSA[$seq][$j]) {
		my $new_pos = $MapMSA[$seq][$j];
		if ($MapMSASquish[$j] =~ /\:$new_pos/) {
		    push(@ConsistentPos,$j+1);
		    $consistencies++;
		} else {
		    push(@InconsistentPos,$j+1);
		    $inconsistencies++;
		}
	    }

	}

	
	# I think it makes sense to hold onto the lists we've built, so I'll
	# just set this to hang around for later debugging output.
	if (1) {

	    # Well, that was easy.  How'd we do? (percent range = [0..100])
	    my $map_accuracy   = int(1000.0 * $consistencies / ($consistencies+$inconsistencies)) / 10.0;
	    my $pct_seq_mapped = int(1000.0 * ($consistencies+$inconsistencies) / $seqlen) / 10.0;
	    
	    # For NOW let's just print these statistics out, since we really haven't
	    # done ANY flippin' verification...
	    my $outstr = "  Percent Seq Mapped: $pct_seq_mapped\%\n  Mapping Accuracy  : $map_accuracy\%\n\n";
	    print "  $gene\n$outstr";

	}

	# Keep in mind that we'll be first just changing the names (no matter
	# what else we do), and then *MAYBE* re-working the MSAs to better
	# represent newly-mapped characters


	##############
	#            #
	#  BIG MOOD  #
	#            #
	##############


	# Regardless of what additional re-tooling we do around finding
	# alternative places to stick characters to ID consistencies,
	# we'll need to identify which characters didn't map (as runs)
	# and add that to the sequence information.
	#
	# Also, this is easy and I'm VERY caffeinated, so I'm blasting
	# some MFing Excision and knocking this basic trash OUT.
	if ($inconsistencies) {

	    # AGAIN, keep in mind that these positions are already '+1'
	    my @UnmappedRuns;
	    my $start_incon_pos = $InconsistentPos[0];
	    my $last_incon_pos  = $start_incon_pos;
	    for (my $i=1; $i<$inconsistencies; $i++) {
		if ($InconsistentPos[$i] != $last_incon_pos+1) {
		    if ($last_incon_pos == $start_incon_pos) { push(@UnmappedRuns,$start_incon_pos);                      }
		    else                                     { push(@UnmappedRuns,$start_incon_pos.'..'.$last_incon_pos); }
		    $start_incon_pos = $InconsistentPos[$i];
		}
		$last_incon_pos = $InconsistentPos[$i];
	    }
	    
	    # Wrap it up with the final run
	    if ($last_incon_pos == $start_incon_pos) { push(@UnmappedRuns,$start_incon_pos);                      }
	    else                                     { push(@UnmappedRuns,$start_incon_pos.'..'.$last_incon_pos); }
	    

	    # Make the string and shove it into the sequence name
	    my $unmapped_str = 'UNMAPPED:'.$UnmappedRuns[0];
	    for (my $i=1; $i<$inconsistencies; $i++) { 
		$unmapped_str = $unmapped_str.','.$UnmappedRuns[$i]; 
	    }
	    
	    my $seqname = $SeqNames[$seq];
	    $seqname =~ /\|([^\|]+)\|[^\|]+$/;
	    my $accession = $1;

	    my $replacement_acc = $unmapped_str.'|'.$accession;
	    $seqname =~ s/\|$accession\|/\|$replacement_acc\|/;

	    $SeqNames[$seq] = $seqname;

	}

	# We'll also want to have a string to concatenate onto this sequence
	# family's hits.out file <-[which might not exist...]
	$partialhitstr = $partialhitstr."Isoform ID : $SeqNames[$seq]\n";
	$partialhitstr = $partialhitstr."Method     : PartialMap\n";
	$partialhitstr = $partialhitstr."Chromosome : $chr\n";
	$partialhitstr = $partialhitstr."Match Pos.s: ";
	for (my $j=0; $j<$msa_len; $j++) {
	    if ($MSA[$seq][$j] =~ /[A-Za-z]/) {
		if ($MapMSA[$seq][$j]) { $partialhitstr = $partialhitstr.$MapMSA[$seq][$j]; }
		else                   { $partialhitstr = $partialhitstr.'-';               }
		$partialhitstr = $partialhitstr.',';
	    }
	}
	$partialhitstr =~ s/\,$//; # We'll have one comma too many at the end of this endeavor
	$partialhitstr = $partialhitstr."\n\n";

    }

    # We are now prepared to cry and scream about where all these partial
    # mappings led us.
    if ($inhitfname !~ /\.hits\.out/) { 
	$inhitfname = $inmsafname;
	$inhitfname =~ s/\.afa$/\.hits\.out/;
    }
    open(my $inhitf,'>>',$inhitfname);
    print $inhitf "$partialhitstr";
    close($inhitf);

    # To keep sequences that we couldn't partially map easily (visually)
    # identifiable, we'll need to make sure that they get written out first,
    # potentially re-ordering sequences in the MSA...
    my $big_out_str = '';

    # First off, make sure we have all of the fully unmapped sequences
    foreach my $seq (@UnmappedSeqIDs) {

	# Scan the list of freshly (partially) mapped sequences to make sure
	# this is a genuinely unmapped sequence.
	my $part_mapped = 0;
	for ($i=0; $i<scalar(@NewlyMappedSeqs); $i++) {
	    if ($seq == $NewlyMappedSeqs[$i]) {
		$part_mapped = 1;
		last;
	    }
	}
	next if ($part_mapped);

	# Uh-oh, looks like somebody didn't even partially map :(
	$big_out_str = $big_out_str.">$SeqNames[$seq]";
	for (my $j=0; $j<$msa_len; $j++) {
	    if ($j % 60 == 0) { $big_out_str = $big_out_str."\n"; }
	    $big_out_str = $big_out_str.$MSA[$seq][$j];
	}
	$big_out_str = $big_out_str."\n\n";

    }

    # Next up, the newly (partially) mapped sequences -- note that they've been
    # upgraded to splice-site-marker status!
    foreach my $seq (@NewlyMappedSeqs) {
	$big_out_str = $big_out_str.">$SeqNames[$seq]";
	for (my $j=0; $j<$msa_len; $j++) {
	    if ($j % 60 == 0) { $big_out_str = $big_out_str."\n"; }
	    if ($MSASquish[$seq][$j] eq '*') { $big_out_str = $big_out_str.'*';            } # All grown up!
	    else                             { $big_out_str = $big_out_str.$MSA[$seq][$j]; }
	}
	$big_out_str = $big_out_str."\n\n";	
    }

    # Finally, tell me what I already knew all over again
    foreach my $seq (@MappedSeqIDs) {
	$big_out_str = $big_out_str.">$SeqNames[$seq]";
	for (my $j=0; $j<$msa_len; $j++) {
	    if ($j % 60 == 0) { $big_out_str = $big_out_str."\n"; }
	    $big_out_str = $big_out_str.$MSA[$seq][$j];
	}
	$big_out_str = $big_out_str."\n\n";
    }

    # Alrighty, time for the big re-write of the book of this gene family
    # -- NOTE: for debugging, let's add a little somethin' somethin'
    open(my $msafile,'>',$inmsafname.'.somethin-somethin');
    print $msafile "$big_out_str";
    close($msafile);

}




1;









#########################################################################
#
#  Function Name: TranslateCodon
#
#  About: Convert a DNA triple to an amino acid. ('x' for 'stop')
#
sub TranslateCodon
{
    
    my $codonref   = shift;
    my @codonArray = @{$codonref};
    
    # In case we didn't get a hold of a full codon (fair game) or
    # something weird slipped in (less fair, but maybe fair),
    # spit an 'X'
    return 'X' if (!(@codonArray && $codonArray[0] && $codonArray[1] && $codonArray[2]));
    my $codon = $codonArray[0].$codonArray[1].$codonArray[2];
    $codon    =  uc($codon);
    $codon    =~ s/U/T/g;
    return 'X' if ($codon =~ /[^ACGT]/);

    
    # Guarantees uppercase
    @codonArray = split('',$codon);

    
    if ($codonArray[0] eq 'A') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "AAA") { return 'K'; }
	    if ($codon eq "AAC") { return 'N'; }
	    if ($codon eq "AAG") { return 'K'; }
	    if ($codon eq "AAT") { return 'N'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "ACA") { return 'T'; }
	    if ($codon eq "ACC") { return 'T'; }
	    if ($codon eq "ACG") { return 'T'; }
	    if ($codon eq "ACT") { return 'T'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "AGA") { return 'R'; }
	    if ($codon eq "AGC") { return 'S'; }
	    if ($codon eq "AGG") { return 'R'; }
	    if ($codon eq "AGT") { return 'S'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "ATA") { return 'I'; }
	    if ($codon eq "ATC") { return 'I'; }
	    if ($codon eq "ATG") { return 'M'; }
	    if ($codon eq "ATT") { return 'I'; }
	}
    }

    if ($codonArray[0] eq 'C') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "CAA") { return 'Q'; }
	    if ($codon eq "CAC") { return 'H'; }
	    if ($codon eq "CAG") { return 'Q'; }
	    if ($codon eq "CAT") { return 'H'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "CCA") { return 'P'; }
	    if ($codon eq "CCC") { return 'P'; }
	    if ($codon eq "CCG") { return 'P'; }
	    if ($codon eq "CCT") { return 'P'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "CGA") { return 'R'; }
	    if ($codon eq "CGC") { return 'R'; }
	    if ($codon eq "CGG") { return 'R'; }
	    if ($codon eq "CGT") { return 'R'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "CTA") { return 'L'; }
	    if ($codon eq "CTC") { return 'L'; }
	    if ($codon eq "CTG") { return 'L'; }
	    if ($codon eq "CTT") { return 'L'; }
	}
    }
    
    if ($codonArray[0] eq 'G') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "GAA") { return 'E'; }
	    if ($codon eq "GAC") { return 'D'; }
	    if ($codon eq "GAG") { return 'E'; }
	    if ($codon eq "GAT") { return 'D'; }
	}
	if ($codonArray[1] eq 'C') {	    
	    if ($codon eq "GCA") { return 'A'; }
	    if ($codon eq "GCC") { return 'A'; }
	    if ($codon eq "GCG") { return 'A'; }
	    if ($codon eq "GCT") { return 'A'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "GGA") { return 'G'; }
	    if ($codon eq "GGC") { return 'G'; }
	    if ($codon eq "GGG") { return 'G'; }
	    if ($codon eq "GGT") { return 'G'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "GTA") { return 'V'; }
	    if ($codon eq "GTC") { return 'V'; }
	    if ($codon eq "GTG") { return 'V'; }
	    if ($codon eq "GTT") { return 'V'; }
	}
    }
    if ($codonArray[0] eq 'T') {
	if ($codonArray[1] eq 'A') {
	    if ($codon eq "TAA") { return 'x'; }
	    if ($codon eq "TAC") { return 'Y'; }
	    if ($codon eq "TAG") { return 'x'; }
	    if ($codon eq "TAT") { return 'Y'; }
	}
	if ($codonArray[1] eq 'C') {
	    if ($codon eq "TCA") { return 'S'; }
	    if ($codon eq "TCC") { return 'S'; }
	    if ($codon eq "TCG") { return 'S'; }
	    if ($codon eq "TCT") { return 'S'; }
	}
	if ($codonArray[1] eq 'G') {
	    if ($codon eq "TGA") { return 'x'; }
	    if ($codon eq "TGC") { return 'C'; }
	    if ($codon eq "TGG") { return 'W'; }
	    if ($codon eq "TGT") { return 'C'; }
	}
	if ($codonArray[1] eq 'T') {
	    if ($codon eq "TTA") { return 'L'; }
	    if ($codon eq "TTC") { return 'F'; }
	    if ($codon eq "TTG") { return 'L'; }
	    if ($codon eq "TTT") { return 'F'; }
	}
    }

    # Weird codon is weird. TO THE BIN WITH YOU!
    return 'X';
    
}






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

    $num_seqs++;

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

    # This probably looks dumb, but for specific cool-dude purposes I'm doing this here
    if (!(-e $hitfilename)) { return(\@MapMSA,\@UsesMap,$num_maps,0,0,0); }

    open(my $hitfile,'<',$hitfilename);
    while (my $line = <$hitfile>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	if ($line =~ /Isoform ID \: (\S+)/) {

	    my $seqname = $1;

	    my $i = 0;
	    while ($SeqNames[$i] ne $seqname) {
		$i++;
		#debugging
		if ($i >= $num_seqs) {
		    print "\n  '$seqname' not in this list:\n";
		    for ($i=0; $i<$num_seqs; $i++) {
			print "   $SeqNames[$i]\n";
		    }
		    die "\n";
		}
	    }

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






#########################################################################
#
#  Function Name: ParseSPALNOutput
#
#  About: The name says it all.  Run a SPALN command and make sense of
#         its wisdom.
#
#  BUT ALSO: This has been modified for PartialMap -- I'll try to just
#         comment-out anything that gets turned off, so in case it ever
#         turns out we need to turn things back on they'll be ready for
#         us.
#
#         I think the direction I'm going to take this is that we do all
#         of the standard error-checking / bug-correction, but don't worry
#         about any issues of length / percent ID -- instead, we'll return
#         a hash where each position in the protein sequence is a key and
#         its mapping on the genome is the value.
#
#  YELLING! I HAVEN'T THOUGHT ABOUT HOW INSERTIONS WILL WORK IN THIS MODEL!!!
#
#    1.  Don't trust percent identities.
#    2.  J = S.
#
sub PMParseSPALNOutput
{
    my ($i,$j,$k);

    my $spalnCmd   = shift;
    my $spalnfname = shift;
    my $revcomp    = shift;
    my $offset     = shift;

    # PM: I'm so bored of these variables
    #my $highscore = shift;
    #my $ChrName = shift;
    #my $seqname = shift;
    #my $prot_len = shift;

    # PM: Need this for 'SelenocysteineCheck'
    my $protfilename = shift;
    
    # NOTE: When we're just looking for SPALN output we can't just
    #       run the command and jump ship, due to how the various
    #       functions depend on score comparison (well, BLAT).
    #       For this reason, we'll set the 'hitstring' to be the full
    #       output AFTER we get the score.
    #
    #my $spalner = shift;

    # Are we timing?
    #my $timing     = shift;
    #my $timingdata = shift;
    #my $timeA;

    my ($line,$hitstring);

    # DEBUGGING
    #my $printspalnCmd = $spalnCmd;
    #$printspalnCmd =~ s/\|$//;
    #system($printspalnCmd);
    # DEBUGGING

    
    #$timeA = [Time::HiRes::gettimeofday()] if ($timing);
    if (system($spalnCmd)) { die "\n  ERROR: Spaln command '$spalnCmd' failed\n\n"; }
    open(my $stdout,'<',$spalnfname) || die "\n  ERROR: Failed to open Spaln output '$spalnfname'\n\n";
    #@{$timingdata}[$timing] += Time::HiRes::tv_interval($timeA) if ($timing);


    # Look for where SPALN has called the start and end of the region.
    my ($range_low,$range_high);
    $line = readline($stdout);
    while (!eof($stdout) && $line !~ /\S+\/(\d+)\-(\d+)/) {
	$line = readline($stdout);
    }

    # If we've hit the end of the file, no point continuing
    if (eof($stdout)) {
	close($stdout);
	return(0,0,1); 
    }


    # Grab the range of positions on the chromosome that we hit in.
    # Note that we want these to be formatted low=start, high=end
    $line =~ /\S+\/(\d+)\-(\d+)/;
    my $start_pos = $1;
    my $end_pos   = $2;


    # Go to the section that lays out the score and identity
    #
    # NOTE: Even though we generally know what's going on with the
    #       complementarity going into our adventures with SPALN,
    #       it at least one case there's been a flip.
    #
    my $saw_complement = 0;
    my ($hitscore,$percent,$true_num_chars);
    while ($line = <$stdout>) {

	# Complementarity marker!
	if ($line =~ /\;C complement\(/) {
	    $saw_complement = 1;
	}

	if ($line =~ /Score \= (\d+)/) {

	    # PM: I don't think we really need any of the information
	    #     stored here -- we'll just keep this while loop to
	    #     make sure we always end up where we want to be.
	    last;
	    
	    # SPALN's score
	    #
	    $hitscore = $1;

	    # It turns out we need to manually confirm that every
	    # character is included, because SPALN loves pranking
	    # me around.
	    #
	    if ($line =~ /Score \= \S+ \S+\, (\d+)\.\d+ \S+\, (\d+)\.\d+ \S+\, (\d+)\.\d+ \S+\, (\d+)\.\d+/) {
		$true_num_chars = $1 + $2 + $3 + $4;
	    } else {

		# It seems that in cases where SPALN gives BAD hits
		# it breaks away from the expected format, so we'll
		# just assume that we want to try again (hence the
		# (0,1) return).
		#
		close($stdout);
		return(0,1,2);
		
	    }

	    last;
	}
    }


    # PM: Shouldn't make much of a difference
    #
    # Check for any unexpected complementarity stuff
    #if ($revcomp) {
    #if (!$saw_complement) {
    #$revcomp = 0;
    #$ChrName =~ s/\[revcomp\]//;
    #}
    #} else {
    #if ($saw_complement) {
    #$revcomp = 1;
    #$ChrName = $ChrName.'[revcomp]';	    
    #}
    #}

    
    # PM: Still don't want an EOF, but we'll take what we can get
    #
    # If we've hit the end of the file or have worse than 97% identity
    # jump ship.
    if (eof($stdout)) {
	close($stdout); 
	return(0,0,3); 
	#} elsif ($true_num_chars < $prot_len) {
	#close($stdout);
	#return(0,1,4); # <- This '1' might give us a second chance (pull in more sequence)
    }
    
    
    # NOW we can go to the real business (the alignment lines)
    $line = readline($stdout); # 'ALIGNMENT'
    while (!eof($stdout) && $line !~ /ALIGNMENT/) {
	$line = readline($stdout);
    }


    # Still not totally stoked on seeing an eof
    if (eof($stdout)) {
	close($stdout);
	return(0,0,5);
    }


    # PM: Not. Relevant. SIR.
    #
    # Here's where we'll gather the SPALN output when we're looking
    # to see the actual SPALN output, for debugging.
    #
    #if ($spalner) {
    #close($stdout);
    #open($stdout,$spalnCmd) || die "\n  ERROR:  Spaln command '$spalnCmd' failed\n\n";
    #$hitstring = 'SEQ: '.$seqname."\n";
    #$hitstring = $hitstring.'CMD: '.$spalnCmd."\n";
    #while ($line = <$stdout>) { $hitstring = $hitstring.$line; }
    #$hitstring = $hitstring."\nEND\n\n";
    #close($stdout);
    #return($hitscore,$hitstring,0);
    #}


    #
    # ALRIGHT, GENTS AND LADIES, it looks like we're actually
    # doing this thing.
    #

    
    # BIG PM NOTE
    # -----------
    #
    # The plan is to let all of this run exactly the way that it would in
    # Quilter, but then we sneak in at the end and re-traverse the output
    # to figure out the specific amino acid coordinates (w.r.t. potentially
    # split-apart sequences) and do some (hopefully) straightforward offset-ery.
    #


    # Scan through the SPALN output constructing 4 arrays:
    #
    #     - An array with the DNA sequence characters
    #     - An array with the protein sequence characters
    #     - An array mapping protein sequence indices to
    #       DNA sequence indices
    #     - An array mapping DNA sequence characters to
    #       positions in the genome
    #
    # The combination of these arrays will allow us to do
    # all of the work we need in terms of generating our
    # final mapping of the protein to the genome, correcting
    # for micro-exons.
    #
    # UPDATE (Nov. 2019): It looks like Spaln will, somewhat
    # regularly, start/end exons with extraneous gaps, so I'm
    # going to try to pay attention to the starts and ends of
    # exons to make sure these don't drag good hits into the
    # mud.
    #
    my $gap_run_len = 0;
    my $starting_exon = 1;
    my $mismatches = 0;
    my @FullNuclSeq;
    my @FullProtSeq;
    my @AAPositions;
    my @NuclPositions;    
    my $trans_line;
    my $full_length = 0;
    my $num_aas     = 0; 
    my $current_pos = $offset;
    my $first_jump  = 0;
    while (!eof($stdout)) {

	# Grab the next line, clean it up, make sure it's
	# meaningful.
	#
	$line = readline($stdout);
	$line =~ s/\n|\r//g;
	next if (!$line);


	# If we're skipping around in the genome, adjust
	# the position variable accordingly.
	#
	if ($line =~ /skip\s+(\d+)\s+nt/) {
	    my $skiplen = int($1);
	    if ($revcomp) { $current_pos -= $skiplen; }
	    else          { $current_pos += $skiplen; }
	    $starting_exon = 1; # Whatever we see next will start an exon
	    $mismatches -= $gap_run_len; # If we had a run of gaps prior to the end of the exon, fugghet about 'em
	    next;
	}


	# If we've made it this far but don't match this
	# terminal format then we're looking at SPALN's
	# translation of the nucleotide sequence
	#
	if ($line !~ /\| \S+\/\d+\-\d+/) {
	    $trans_line = $line;
	    next;
	}

	my $nn_line = $line;

	# EXCELLENT!  This is the next line of DNA characters
	#
	my @NextNucls = split(//,$line);

	# We also split up what would be the translation, so that
	# we can evaluate percent ID for ourselves.
	#
	my @TransLine = split(//,$trans_line);

	# The final check we need to do is to add/subtract
	# the initial 'jump' into the hit (but only if this
	# is the very beginning)
	#
	if ($first_jump == 0) {

	    if ($line =~ /^\s*(\d+)/) {

		$first_jump = $1;

		# Special catch for BLAT parsing
		if ($revcomp==3) { $current_pos  = ($start_pos+1)-$first_jump; }
		else             { $current_pos += $first_jump;                }

		$first_jump = 1;

	    } else {
		return(0,0,6);
	    }

	}

	
	# ... which means that the following line is the next
	# line of protein characters.
	#
	$line = readline($stdout);
	$line =~ s/\n|\r//g;
	my @NextAAs = split('',$line);


	# Advance to the actual DNA sequence.
	# NOTE:  it's possible for SPALN to have a line of
	#        '-' characters, which we'd want to throw
	#        away.  BUT we still need to count up the
	#        position indices when we encounter these. << Do we?
	#
	$i = 0;
	#while ($i < @NextNucls && $NextNucls[$i] !~ /[a-zA-Z]|\|/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] =~ /\s/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] !~ /\s/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] =~ /\s/) { $i++; }
	


	# Some sort of SPALN weirdness...
	#
	next if ($i == @NextNucls || $NextNucls[$i] eq '|');


	# Run along the two lines storing all of the relevant
	# information that we can get out of them.
	#
	while ($NextNucls[$i] !~ /\|/) {
	    
	    push(@FullNuclSeq,$NextNucls[$i]);
	    push(@NuclPositions,$current_pos);

	    # If we're centered on a codon, take note.
	    #
	    if ($NextAAs[$i] =~ /\w/) {		
		push(@FullProtSeq,$NextAAs[$i]);
		push(@AAPositions,$full_length);
		$num_aas++;
	    }
	    
	    # Also, note a *SUBSTITUTION* mismatch (sorry for yelling)
	    if ($NextAAs[$i] =~ /\S/
		&& $NextAAs[$i] ne $TransLine[$i] 
		&& !(uc($TransLine[$i]) eq 'J' && uc($NextAAs[$i]) eq 'S')) { # Sometimes SPALN calls an 'S' a 'J'
		if ($NextAAs[$i] eq '-') {
		    # Gap mismatch (could be a goof...)
		    if (!$starting_exon) {
			$gap_run_len++;
			$mismatches++;
		    }
		} else {
		    $mismatches++; # Genuine mismatch
		}
	    } else {
		$gap_run_len=0;
	    }
	    # I'm using that format because we might see more delightful tricks


	    # Increment the position in the genome
	    #
	    if ($NextNucls[$i] =~ /\w/) {
		if ($revcomp) { $current_pos--; }
		else          { $current_pos++; }
	    }

	    $full_length++;
	    $i++;
	    
	}
	
    }
    

    # Done reading input!  Now we just need to verify that
    # any and all micro-exons get cleaned up!
    #
    close($stdout);


    # If we aren't satisfied with the overall percent identity of the
    # hit go back (possibly to pull in more sequence and try again).
    #
    # * * * THIS IS THE PERCENT IDENTITY CHECK POSITION * * * 
    #
    #return (0,1,7) if ($mismatches/$prot_len > 0.05);
    

    # Selenocysteine does some crazy stuff, man.  I translated some once
    # and I'm still coming down.
    #
    #if ($num_aas != $prot_len) { # <--- Start of original Seleno conditional
	
    # PM: I think we'll just do this automatically, instead of conditionally.
    my ($fps_ref,$aap_ref) = SelenocysteineCheck(\@FullNuclSeq,\@FullProtSeq,\@AAPositions,$protfilename);
    @FullProtSeq = @{$fps_ref};
    @AAPositions = @{$aap_ref};
    $num_aas = scalar(@FullProtSeq);

    # PM: It's still good to check for Selenocysteine, but we no
    #     don't mind so much about our mapping covering the full
    #     protein sequence (well, technically partial, but whatever)
    #return (0,0,8) if ($num_aas != $prot_len);

    #} # <--- End of original Seleno conditional


    # ********************************************************************
    # *** 'IF ZERO' WILL TOGGLE MICRO EXON SEARCHING OFF *****************
    # ********************************************************************
    #
    #
    #if (0) {

    # Run along SPALN's output looking for 'exons' with
    # improbably low amino acid counts.  Currently, we
    # say that an exon should have at least 4 amino acids.
    #
    my @MicroExonRuns;
    my $micro_start = -1;
    my $micro_end   = 0;
    
    $i = 0;
    while ($i < $num_aas) {


	# Figure out the genome position of the current amino
	# acid, so that we can infer exon boundaries (large
	# changes in nucleotide position).
	#
	my $aa_loc   = $NuclPositions[$AAPositions[$i]];
	my $aa_start = $i;
	my $exon_len = 1;

	
	# Run through this exon, counting how many amino acids
	# it has.
	#
	$i++;
	while ($i < $num_aas && abs($NuclPositions[$AAPositions[$i]] - $NuclPositions[$AAPositions[$i-1]]) < 5) {
	    $exon_len++;
	    $i++;
	}


	# If we have too few amino acids, then we're working
	# with a micro exon, baby!  Observe that we chain together
	# runs of contiguous micro-exons.
	#
	if ($exon_len < 4) {
	    
	    if ($micro_start == -1) { $micro_start = $aa_start; }
	    $micro_end = $i;
	    
	} elsif ($micro_start != -1) {
	    
	    push(@MicroExonRuns,$micro_start.'-'.$micro_end);
	    $micro_start = -1;
	    
	}
	
    }
    
    
    # Next we'll run through each of the indicated micro-exon
    # regions ('regions' because they can be chained together)
    # and try to re-align the amino acids to the ends of the
    # 'major' exons that flank them.
    #
    # NOTE: Be careful that we don't step over a breakpoint by
    #       checking whether NuclPositions[i] and NP[i+1] are
    #       1 off from one another.
    #
    foreach my $micro_run_str (@MicroExonRuns) {

	
	# Split the indices out of the string
	#
	my @MicroRun = split(/\-/,$micro_run_str);
	$micro_start = $MicroRun[0];
	$micro_end   = $MicroRun[1]; # One past actual end, so we can do strict '<'


	# Identify the amino acid sequence that constitutes
	# this micro exon (or chain of micro exons)
	#
	my @MicroAAs;
	for ($j = $micro_start; $j < $micro_end; $j++) {
	    push(@MicroAAs,$FullProtSeq[$j]);
	}

	# Figure out how long the micro exon is (keeping in mind that '_end' is
	# one position past the actual end index)
	#
	my $micro_len = $micro_end-$micro_start;

	
	# Look at the preceding exon
	#
	my @RearExt;
	my $rear_ext_len = 0;
	if ($micro_start > 0) {

	    
	    # Get situated so that we can start trying to
	    # transfer amino acids from our micro exons to
	    # the back end of the preceding exon.
	    #
	    # NOTE (for making sense of this):
	    # --------------------------------
	    # NuclPositions are the chromosome indices, AAPositions are the
	    # array indices that SPALN has identified as codon centers -- thus
	    # looking up the NuclPosition of an AAPosition gets you the genome
	    # index for the AAPosition's codon center.
	    #
	    my $rear_codon_center = $NuclPositions[$AAPositions[$micro_start-1]];
	    my $nuclseq_pos = $AAPositions[$micro_start-1]+3; # Which array index are we centered around?
	    while ($rear_ext_len < $micro_len) {
		
		# Make sure we don't overstep (by checking if we've passed over
		# an stretch of genome indices indicating an intron)
		#
		my $step_check = $NuclPositions[$nuclseq_pos];
		if ($revcomp) { $step_check += 3 * ($rear_ext_len+1); }
		else          { $step_check -= 3 * ($rear_ext_len+1); }
		
		last if ($step_check != $rear_codon_center);


		# Figure out what the next codon is
		#
		my @NextCodon = ();
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos-1]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos+1]);

		
		# If we have a match we are very happy.  If
		# we do not have a match we are very sad.
		#
		# Why would we do 'micro_start+(rear_ext_len+1)'?
		#
		if (uc($FullProtSeq[$micro_start+$rear_ext_len]) eq uc(TranslateCodon(\@NextCodon))) {
		    push(@RearExt,$nuclseq_pos);
		    $nuclseq_pos += 3;
		    $rear_ext_len++;
		} else {
		    last;
		}
		
	    }
	    
	}
	
	# If we found a complete run, we're sooooo stoked!
	#
	if ($rear_ext_len == $micro_len) {

	    # Record our triumph
	    #
	    for ($i=0; $i<$micro_len; $i++) { $AAPositions[$micro_start+$i] = $RearExt[$i]; }
	    next;
	    
	}

	
	# Look at the exon ahead
	#
	my @FwdExt;
	my $fwd_ext_len = 0;
	if ($micro_end < $num_aas) {
	    
	    # Try to find a good extension to the front of
	    # the upcoming exon.
	    #
	    my $fwd_codon_center = $NuclPositions[$AAPositions[$micro_end]];
	    my $nuclseq_pos = $AAPositions[$micro_end]-3;
	    while ($fwd_ext_len < $micro_len) {

		
		# Make sure we don't overstep
		#
		my $step_check = $NuclPositions[$nuclseq_pos];
		if ($revcomp) { $step_check -= 3 * ($fwd_ext_len+1); }
		else          { $step_check += 3 * ($fwd_ext_len+1); }
		
		last if ($step_check != $fwd_codon_center);


		# What does the next codon encode?
		#
		my @NextCodon = ();
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos-1]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos+1]);
		

		# If we have a match, in terms of what the codon
		# encodes, then we can continue our considerations.
		# Otherwise we cannot.
		#
		# Note that the first nucleotide position in the
		# 'FwdExt' array corresponds to the last amino acid
		# in the micro-exon region.
		#
		if (uc($FullProtSeq[$micro_end-($fwd_ext_len+1)]) eq uc(TranslateCodon(\@NextCodon))) {
		    push(@FwdExt,$nuclseq_pos);
		    $nuclseq_pos -= 3;
		    $fwd_ext_len++;
		} else {
		    last;
		}
		
	    }
	    
	}

	
	# Did we find a complete run?
	#
	if ($fwd_ext_len == $micro_len) {
	    
	    for ($i=0; $i<$micro_len; $i++) { $AAPositions[$micro_end-($i+1)] = $FwdExt[$i]; }
	    next;
	    
	}
	

	# Append as much as you can to the preceding
	# major exon first.  We don't put any work into
	# evaluating whether amino acids that could go
	# either way should be fitted to one or the other
	# side.
	#
	$i = 0;
	while ($i < $rear_ext_len) {
	    $AAPositions[$micro_start+$i] = $RearExt[$i];
	    $i++;
	}

	
	# Skip ahead to the stuff we're appending to the
	# upcoming major exon, if you need to.
	#
	if ($i < $micro_len - $fwd_ext_len) { $i = $micro_len - $fwd_ext_len; }


	# Append this stuff onto the forward friend.
	#
	$j = $fwd_ext_len - 1;
	while ($i < $micro_len) {
	    $AAPositions[$micro_end-($micro_len-$i)] = $FwdExt[$j];
	    $i++;
	    $j--;
	}
	
    }
    

    #}
    #
    #
    # ********************************************************************
    # *** 'IF ZERO' TOGGLE END *******************************************
    # ********************************************************************

    #
    # PM: This will be our largest departure from the original approach.
    #     instead of a list, we'll pair each position in the amino acid
    #     sequence with its codon center in a hash (since we don't have a
    #     1..seqlen mapping guaranteed).
    #
    #     HOWEVER, we'll benefit from making an initial list like the
    #     original would have outputted, and then rewinding to get the
    #     actual amino indices.
    #
    # *** WHEN THERE ARE INSERTIONS: We'll have to stack them.
    #

    # Radical, dude!  Now all we need to do is just convert
    # our list of amino acid positions into a list of
    # their corresponding middle amino acids (with '*'s to
    # denote splice site boundaries).
    #
    # We'll try to be careful about how we place codons that are
    # called as insertions into the genome.
    #
    my $prev_codon_c = $NuclPositions[$AAPositions[0]];
    my @CodonCenters;
    push(@CodonCenters,$prev_codon_c);
    for ($i = 1; $i < $num_aas; $i++) {

	# A gap of more than 4 is taken to deserve a splice site marker
	#
	#if (abs($NuclPositions[$AAPositions[$i]] - $NuclPositions[$AAPositions[$i-1]]) > 4) {
	#$CodonCenters = $CodonCenters.',*';
	#}

	# Prepping for a fresh character!
	#
	#$CodonCenters = $CodonCenters.',';
	
	# Are we in an insertion (wrt genome)?
	#
	if ($FullNuclSeq[$AAPositions[$i]] !~ /\w/) {
	    #$CodonCenters = $CodonCenters.$prev_codon_c;
	    push(@CodonCenters,$prev_codon_c);
	} else {
	    $prev_codon_c = $NuclPositions[$AAPositions[$i]];
	    #$CodonCenters = $CodonCenters.$prev_codon_c;
	    push(@CodonCenters,$prev_codon_c);
	}
	
    }
    
    # PM: OK, time to convert this list to a hash, using smarts and intelligents.
    #     More specifically, we'll rewind the file and run along the alignment
    #     using Spaln's signposts to tell us which aminos in the sequence have
    #     actually been incorporated into the alignment.
    #
    #     Happily, we know that this file FOLLOWS THE DAMN RULES, so we can move
    #     relatively quickly

    open($stdout,'<',$spalnfname);

    # Scan to the start of the alignment
    while ($line = <$stdout>) {
	last if ($line =~ /ALIGNMENT/);
    }

    my %AminoIndexToGenomePos;
    my $codon_centers_pos = 0;
    while ($line = <$stdout>) {

	# Every time we see a nucleotide line, we have our bearings.
	if ($line =~ /\| \S+\/\d+\-\d+/) { # nucl line

	    $line = <$stdout>; # INPUT PROTEIN LINE!
	    $line =~ /^\s*(\d+)([^\|]+)\|/;

	    my $amino_id = $1 - 1; # We'll work with 0..[seqlen-1] indices
	    my $aminostr = $2;

	    # Break up the line and attribute the appropriate codon center
	    # to each amino index;
	    $aminostr =~ s/[^A-Za-z]//g;
	    foreach my $amino (split(//,$aminostr)) {
		$AminoIndexToGenomePos{$amino_id} = $CodonCenters[$codon_centers_pos];
		$amino_id++;
		$codon_centers_pos++;
	    }

	}

    }

    close($stdout);

    return(\%AminoIndexToGenomePos);

    # Now that we've got our purty little hit, we write it on out and tell
    # 'em how happy that made us feel.
    #my $hitline = "Isoform ID : $seqname\n";
    #$hitline    = $hitline."Method     : SPALN\n";
    #$hitline    = $hitline."Chromosome : $ChrName\n";
    #$hitline    = $hitline."Match Pos.s: $CodonCenters\n\n";
    #return ($hitscore,$hitline,9);

}








#########################################################################
#
#  Function Name: SelenocysteineCheck
#
#  About: This function checks if a selonocyteine (U, encoded by TGA)
#         explains an observed absence of aminos in SPALN output.
#
sub SelenocysteineCheck
{
    my $nucl_seq_ref = shift;
    my $prot_seq_ref = shift;
    my $aa_pos_ref   = shift;
    my $protfilename = shift;

    my @NuclSeq = @{$nucl_seq_ref};
    my @ProtSeq = @{$prot_seq_ref};
    my @AAPos   = @{$aa_pos_ref};

    my @TrueProtSeq;
    open(my $protfile,'<',$protfilename) || die "\n  ERROR:  Failed to open file '$protfilename' during selenocysteine check\n\n";
    while (my $line = <$protfile>) {
	$line =~ s/\n|\r//g;
	next if ($line =~ /\>/);
	foreach my $amino (split(//,$line)) { push(@TrueProtSeq,uc($amino)); }
    }
    close($protfile);


    my $aa_pos = 0;
    while ($aa_pos < scalar(@TrueProtSeq)) {

	# Looks like we've discovered a selenocysteine in our midst!
	if ($TrueProtSeq[$aa_pos] eq 'U' && $aa_pos < scalar(@ProtSeq) && $ProtSeq[$aa_pos] ne 'U') {

	    # Is this consistent with stepping forward 3 from the previous
	    # amino's position?  3 back from forthcoming amino?
	    my $nailed_it = 0;
	    if ($aa_pos) {
		
		# Should be to the right of the last position...
		my $nucl_pos = $AAPos[$aa_pos-1] + 2;
		my $codon = uc($NuclSeq[$nucl_pos].$NuclSeq[$nucl_pos+1].$NuclSeq[$nucl_pos+2]);
		if ($codon eq 'TGA') {

		    # Note:  Splice moves everything at the indicated position to
		    #        the right and sticks in what you want there.
		    splice(@ProtSeq,$aa_pos,0,'U');
		    splice(@AAPos,$aa_pos,0,$nucl_pos+1);
		    $nailed_it = 1;
		} 

	    } 

	    if (!$nailed_it && $aa_pos < scalar(@ProtSeq)-1) {
		
		# ... Or left of the next position
		my $nucl_pos = $AAPos[$aa_pos] - 4;
		my $codon = uc($NuclSeq[$nucl_pos].$NuclSeq[$nucl_pos+1].$NuclSeq[$nucl_pos+2]);
		if ($codon eq 'TGA') {
		    splice(@ProtSeq,$aa_pos-1,0,'U');
		    splice(@AAPos,$aa_pos-1,0,$nucl_pos+1);
		} 

	    }

	}

	$aa_pos++;

    }

    return (\@ProtSeq,\@AAPos);

}






#########################################################################
#
#  FUNCTION: UpdateMapMSA
#
sub UpdateMapMSA
{

    my $msaref     = shift;
    my $mapmsaref  = shift;
    my $maphashref = shift;

    my @MSA     = @{$msaref};
    my @MapMSA  = @{$mapmsaref};
    my %MapHash = %{$maphashref};

    my $seqid   = shift;
    my $msa_len = shift;

    my $seq_char_id = 0;
    for (my $i=0; $i<$msa_len; $i++) {
	if ($MSA[$seqid][$i] =~ /[A-Za-z]/) {
	    if ($MapHash{$seq_char_id}) {
		$MapMSA[$seqid][$i] = $MapHash{$seq_char_id};
	    }
	    $seq_char_id++;
	}
    }

    return \@MapMSA;

}







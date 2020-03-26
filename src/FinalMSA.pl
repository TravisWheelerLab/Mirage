#!/usr/bin/env perl
#
# FinalMSA.pl - Alex Nord - 2016
#
# ABOUT: This script takes an MSA with splice site markers and
#        cleans it up to make mirage's finished product.
#
use warnings;
use strict;
use POSIX;


sub MAX;
sub MIN;
sub RemoveIntronGaps;
sub PostHocCleanup;



# Because this script only makes sense within the context of mirage,
# we don't offer additional options.
if (@ARGV != 2) { die "\n  USAGE:  perl FinalMSA.pl <in> <out>\n\n"; }


my ($i,$j,$k);


# Where we'll store the sequence names and the MSA
my @Headers;
my @BigMSA;

# Length of the seq.s in the MSA
my $MSALength;

# Ripping the sequences from the input file, storing them in BigMSA
open(my $inf,'<',$ARGV[0]) || die "  Failed to open file '$ARGV[0]'\n";
my $numSeqs = 0;
my $line = <$inf>;
while (!eof($inf)) {

    $line =~ s/\n|\r//g;
    while (!eof($inf) && $line !~ /^\>/) {
	$line =  <$inf>;
	$line =~ s/\n|\r//g;
    }

    last if (eof($inf));

    # Make sure we like how the header looks
    my $header = $line;
    if ($header !~ /^>GN\:/) {
	$header =~ s/\>//g;
	$header = '>GN:'.$header;
    }

    $Headers[$numSeqs] = $header;

    $i = 0;
    $line = <$inf>;
    while ($line !~ /^\>/) {
	$line =~ s/\n|\r//g;
	if ($line) {
	    my @seqline = split('',$line);
	    foreach $j (@seqline) {
		$BigMSA[$numSeqs][$i] = $j;
		$i++;
	    }
	}
	last if (eof($inf));
	$line  = <$inf>;
    }
    $MSALength = $i;

    $numSeqs++;

}
close($inf);


# Remove intron gaps from the MSA
#
my ($MSAref,$FinalLength) = RemoveIntronGaps(\@BigMSA,$MSALength,$numSeqs);
my @FinalMSA = @{$MSAref};


# Our post-hoc cleanup -- correcting 'jigsaw' ugliness and other
# obvious imperfections.
#
if ($numSeqs > 1) {
    my ($cleanedMSA,$cleanLength) = PostHocCleanup(\@FinalMSA,$FinalLength,$numSeqs);
    @FinalMSA = @{$cleanedMSA};
    $FinalLength = $cleanLength;
}


# For one final thing, we do a check to see if there's been some 'M'
# misplacement.  This can occur as part of a MultiSeqNW error check that's
# intended to 'cascade' through a series of characters, but can cause
# the first 'M' to get shifted off to the left.
#
# We'll default to moving an initial 'M' followed by gap characters to the
# right, since it would instantiate a bizarre micro-exon or gappy thing
# otherwise <--- on the condition that there's already at least one 'M'
#                at the new start position in the MSA
#
for (my $i=0; $i<$numSeqs; $i++) {

    if (uc($FinalMSA[$i][0]) eq 'M' && $FinalMSA[$i][1] eq '-') {

	my $j=2;
	while ($FinalMSA[$i][$j] eq '-') {
	    $j++;
	}
	$j--;

	my $m_support = 0;
	for (my $k=0; $k<$numSeqs; $k++) {
	    if (uc($FinalMSA[$k][$j]) eq 'M') {
		$m_support = 1;
		last;
	    }
	}

	if ($m_support) {
	    $FinalMSA[$i][$j] = $FinalMSA[$i][0];
	    $FinalMSA[$i][0]  = '-';
	}

    }
}


# Print out our final MSA
open(my $outf,'>',$ARGV[1]);
foreach $i (0..$numSeqs-1) {
    $j = 0;
    print $outf "$Headers[$i]\n";
    while ($j < $FinalLength) {
	print $outf "$FinalMSA[$i][$j]";
	$j++;
	print $outf "\n" if ($j % 50 == 0);
    }
    print $outf "\n" if ($j % 50);
}
close($outf);


# Yeah! Yeah! Yeah! Yeah!
# What a rush!  Let's go again!
1;




#######################
#                     #
#    END OF SCRIPT    #
#                     #
#######################





#
# MAX
#
sub MAX
{
    my $a = shift;
    my $b = shift;
    return $a if ($a > $b);
    return $b;
}






#
# MIN
#
sub MIN
{
    my $a = shift;
    my $b = shift;
    return $a if ($a < $b);
    return $b;
}




########################################################################
#
# Function: RemoveIntronGaps
#
# About: This function removes all splice site columns from an MSA.
#
sub RemoveIntronGaps
{

    my ($i, $j, $k);

    my $MSAref = shift;
    my @BigMSA = @{$MSAref};

    my $MSALength = shift;
    my $num_seqs  = shift;
    
    my @CleanMSA;
    my $CleanLength = 0;
    for ($j=0; $j<$MSALength; $j++) {

	my $splice_site = 0;
	my $gap_count   = 0;

	for ($i=0; $i<$num_seqs; $i++) {

	    # It's a splice boundary marker
	    if ($BigMSA[$i][$j] eq '*') {
		$splice_site = 1;
		last;
	    }

	    # There's a gap
	    if ($BigMSA[$i][$j] eq '-') { $gap_count++; }

	    # Record this character
	    $CleanMSA[$i][$CleanLength] = $BigMSA[$i][$j];
	    
	}

	# If this is a splice site column or all gaps
	# we can overwrite it.
	if ($splice_site == 0 && $gap_count < $num_seqs) { $CleanLength++; }
	
    }

    return (\@CleanMSA,$CleanLength);

}




########################################################################
#
# Function: PostHocCleanup
#
sub PostHocCleanup
{

    my $msa_ref  = shift;
    my $msa_len  = shift;
    my $num_seqs = shift;

    my @MSA = @{$msa_ref};

    my %PositionCharProfile;
    
    #
    # What we're doing is scanning through the MSA looking for
    # gap starts and checking
    #
    #   [1.] if fewer than half of the sequences are being represented
    #
    #   [2.] a strongly matching (<2 mismatches) subsequence exists in one
    #        of the other sequences at that sequence's last gap start or
    #        current gap end.
    #
    #   [3.] confirming that we have a gap to that strongly matching
    #        subsequence.
    #

    #
    #  FORWARD Pass
    #
    my @LastNonGap;
    for (my $i=0; $i<$num_seqs; $i++) { $LastNonGap[$i] = -1; }
    for (my $pos=0; $pos<$msa_len-1; $pos++) {

	my @EndingGap;
	for (my $seq_id = 0; $seq_id < $num_seqs; $seq_id++) {

	    if ($MSA[$seq_id][$pos] eq '-') {

		if ($pos < $msa_len-1 && $MSA[$seq_id][$pos+1] =~ /[A-Za-z]/) {

		    push(@EndingGap,$seq_id);
		    
		    my $hash_key = ($pos+1).uc($MSA[$seq_id][$pos+1]);
		    if ($PositionCharProfile{$hash_key}) { $PositionCharProfile{$hash_key}++; }
		    else                                 { $PositionCharProfile{$hash_key}=1; }

		}

	    } else {

		my $hash_key = $pos.uc($MSA[$seq_id][$pos]);
		if ($PositionCharProfile{$hash_key}) { $PositionCharProfile{$hash_key}++; }
		else                                 { $PositionCharProfile{$hash_key}=1; }

		$LastNonGap[$seq_id] = $pos;

	    }
	}

	# Did fewer than half of the sequences see the end of a gapped region?
	#
        if (@EndingGap && scalar(@EndingGap) <= int($num_seqs/2)) {

	    # Before we consider individual amino acid movements, we'll check whether
	    # this is a "jigsaw" region by checking if [1.] everyone is the same
	    # amino acid and [2.] the same character is what we see at the first gap
	    # index.
	    #
	    my $ref_char   = uc($MSA[$EndingGap[0]][$pos+1]);
	    my $ref_pos    = $LastNonGap[$EndingGap[0]]+1;
	    my $jigsawable = 1;
	    my $all_gap    = 1;
	    for (my $gap_seq_index = 1; $gap_seq_index < scalar(@EndingGap); $gap_seq_index++) {
		my $seq_id = $EndingGap[$gap_seq_index];
		if (uc($MSA[$seq_id][$pos+1]) ne $ref_char) {
		    $jigsawable = 0;
		    last;
		} elsif ($ref_pos != $LastNonGap[$seq_id]+1) {
		    $all_gap = 0;
		    last;
		}
	    }

	    
	    # One last check before the full jigsaw transfer -- are we in agreement that the
	    # reference position suggests this character?
	    #
	    if ($jigsawable && $all_gap && $PositionCharProfile{$ref_pos.$ref_char}
		&& $PositionCharProfile{$ref_pos.$ref_char} >= scalar(@EndingGap)) {
		
		foreach my $seq_id (@EndingGap) {

		    $MSA[$seq_id][$ref_pos] = $MSA[$seq_id][$pos+1];
		    $MSA[$seq_id][$pos+1]   = '-';

		    $LastNonGap[$seq_id] = $ref_pos;

		    $PositionCharProfile{$ref_pos.$ref_char}++;
		    $PositionCharProfile{($pos+1).$ref_char}--;

		}

	    } else {


		# For each sequence that saw the end of a gap, check if the upcoming
		# character matches the "character profile" of anything at the first
		# position of the gap.
		my @DelayedRemovalKeys;
		foreach my $seq_id (@EndingGap) {
		    
		    my $upcoming_char = $MSA[$seq_id][$pos+1];
		    my $first_gap_pos = $LastNonGap[$seq_id]+1;
		    my $new_hash_key  = $first_gap_pos.uc($upcoming_char);
		    my $old_hash_key  = ($pos+1).uc($upcoming_char);
		    if ($first_gap_pos && $PositionCharProfile{$new_hash_key}
			&& $PositionCharProfile{$new_hash_key} > $PositionCharProfile{$old_hash_key}) {
			
			$MSA[$seq_id][$pos+1] = '-';
			$MSA[$seq_id][$first_gap_pos] = $upcoming_char;
			
			$LastNonGap[$seq_id]++;
			
			$PositionCharProfile{$new_hash_key}++;
			$PositionCharProfile{$old_hash_key}--;
			
		    } else {
			
			push(@DelayedRemovalKeys,$old_hash_key);
			
		    }
		    
		}
		
		foreach my $key (@DelayedRemovalKeys) { $PositionCharProfile{$key}--; }
		
	    }

	}
    }


    #
    #  BACKWARD Pass
    #
    for (my $i=0; $i<$num_seqs; $i++) { $LastNonGap[$i] = -1; }
    for (my $pos=$msa_len-1; $pos>0; $pos--) {

	my @StartingGap;
	for (my $seq_id = 0; $seq_id < $num_seqs; $seq_id++) {

	    if ($MSA[$seq_id][$pos] eq '-') {

		if ($MSA[$seq_id][$pos-1] =~ /[A-Za-z]/) {
		    push(@StartingGap,$seq_id);
		}

	    } else {

		$LastNonGap[$seq_id] = $pos;

	    }
	}

	# Did fewer than half of the sequences see the end of a gapped region?
	#
        if (@StartingGap && scalar(@StartingGap) <= int($num_seqs/2)) {

	    # Before we consider individual amino acid movements, we'll check whether
	    # this is a "jigsaw" region by checking if [1.] everyone is the same
	    # amino acid and [2.] the same character is what we see at the first gap
	    # index.
	    #
	    my $ref_char   = uc($MSA[$StartingGap[0]][$pos-1]);
	    my $ref_pos    = $LastNonGap[$StartingGap[0]]-1;
	    my $jigsawable = 1;
	    my $all_gap    = 1;
	    for (my $gap_seq_index = 1; $gap_seq_index < scalar(@StartingGap); $gap_seq_index++) {
		my $seq_id = $StartingGap[$gap_seq_index];
		if (uc($MSA[$seq_id][$pos-1]) ne $ref_char) {
		    $jigsawable = 0;
		    last;
		} elsif ($ref_pos ne $LastNonGap[$seq_id]-1) {
		    $all_gap = 0;
		    last;
		}
	    }

	    
	    # One last check before the full jigsaw transfer -- are we in agreement that the
	    # reference position suggests this character?
	    #
	    if ($jigsawable && $all_gap && $PositionCharProfile{$ref_pos.$ref_char}
		&& $PositionCharProfile{$ref_pos.$ref_char} >= scalar(@StartingGap)) {
		
		foreach my $seq_id (@StartingGap) {

		    $MSA[$seq_id][$ref_pos] = $MSA[$seq_id][$pos-1];
		    $MSA[$seq_id][$pos-1]   = '-';

		    $LastNonGap[$seq_id] = $ref_pos;

		    $PositionCharProfile{$ref_pos.$ref_char}++;
		    $PositionCharProfile{($pos-1).$ref_char}--;

		}

	    } else {

		
		# For each sequence that saw the end of a gap, check if the upcoming
		# character matches the "character profile" of anything at the first
		# position of the gap.
		foreach my $seq_id (@StartingGap) {
		    
		    my $upcoming_char = $MSA[$seq_id][$pos-1];
		    my $first_gap_pos = $LastNonGap[$seq_id]-1;
		    my $new_hash_key  = $first_gap_pos.uc($upcoming_char);
		    my $old_hash_key  = ($pos-1).uc($upcoming_char);
		    if ($first_gap_pos && $PositionCharProfile{$new_hash_key}
			&& $PositionCharProfile{$new_hash_key} > $PositionCharProfile{$old_hash_key}) {
			
			$MSA[$seq_id][$pos-1] = '-';
			$MSA[$seq_id][$first_gap_pos] = $upcoming_char;
			
			$LastNonGap[$seq_id]--;
			
			$PositionCharProfile{$new_hash_key}++;
			$PositionCharProfile{$old_hash_key}--;
			
		    }
		    
		}

	    }

	}

    }

    
    # We may have created some all-gap columns, so let's go ahead and blow 'em
    # to smithereens (sp?).
    #
    my @CleanMSA;
    my $clean_len = 0;
    for (my $pos=0; $pos<$msa_len; $pos++) {

	my $all_gaps = 1;
	for (my $seq_id=0; $seq_id<$num_seqs; $seq_id++) {
	    $CleanMSA[$seq_id][$clean_len] = $MSA[$seq_id][$pos];
	    $all_gaps = 0 if ($MSA[$seq_id][$pos] =~ /[A-Za-z]/);
	}
	
	$clean_len++ if ($all_gaps == 0);

    }

    return (\@CleanMSA,$clean_len);
    
}



#####################################################
##################                 ##################
##################   END OF FILE   ##################
##################                 ##################
#####################################################















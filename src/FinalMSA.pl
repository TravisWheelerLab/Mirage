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
sub CleanMSA; # This is going to be ripped almost directly from MultiMSA...



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

# Do any cleanup that we need to do.
#
# NO CLEANING ALLOWED
#
#my ($cleanMSAref,$clean_len) = CleanMSA(\@BigMSA,$numSeqs,$MSALength);
#@BigMSA    = @{$cleanMSAref};
#$MSALength = $clean_len;

# Remove intron gaps from the MSA, correcting any ugly shifts
my ($MSAref,$FinalLength) = RemoveIntronGaps(\@BigMSA,$MSALength,$numSeqs);
my @FinalMSA = @{$MSAref};


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








################################################################
#
# FUNCTION: CleanMSA <-- Ripped directly from MultiMSA
#
sub CleanMSA
{

    # The MSA that we're interested in cleaning up
    #
    my $InputMSAref = shift;
    my @InputMSA = @{$InputMSAref};

    # The dimensions of the MSA
    #
    my $num_seqs = shift;
    my $input_len = shift;



    ################
    #              #
    #   PART ONE   #
    #              #
    ################



    #
    # The aim of the first part of the cleanup is to locate
    # any segments in the alignment where we have mismatched
    # sequence.  This is done in a forward left-shifting pass
    # and a reverse right-shifting pass.
    #


    # First, the forward pass.  We certainly can't shift anyone
    # in column 0 leftwards.
    #
    for ($j=1; $j<$input_len; $j++) {


	# We don't fool around with splice sites.  Similarly, we
	# can't do a left shift if the last position was a splice
	# site...
	#
	next if ($InputMSA[0][$j] eq '*' || $InputMSA[0][$j-1] eq '*');


	# Is there a consensus character?  Is there disagreement?
	#
	my $sample_char  = 0;
	my $disagreement = 0;
	my @OpenLeft;
	my @ClosedLeft;
	my %ValidChars;
	for ($i=0; $i<$num_seqs; $i++) {

	    # Are we open to the idea of a leftwards shift?
	    #
	    if ($InputMSA[$i][$j-1] eq '-') {
		push(@OpenLeft,$i);
	    } else {
		push(@ClosedLeft,$i);
		$ValidChars{uc($InputMSA[$i][$j-1])} = 1;
	    }

	    # Don't distract me with your dumb old gaps!
	    #
	    next if ($InputMSA[$i][$j] eq '-');

	    # Normalize to upper case
	    #
	    my $this_char = uc($InputMSA[$i][$j]);

	    # Is this the first non-gap we've seen?
	    #
	    if (!$sample_char) {
		$sample_char = $this_char;

	    } elsif ($this_char ne $sample_char) { # COULD IT BE?!
		$disagreement = 1;

	    }

	}

	# Did we find any disagreements within this column?  Is fun on the table?
	#
	if ($disagreement && scalar(@OpenLeft)) { # Here comes the fun!

	    # For each character next to a gap, if we have a match scoot on over!
	    #
	    foreach $k (@OpenLeft) {

		my $this_char = uc($InputMSA[$k][$j]);

		next if ($this_char eq '-');

		# Shift left!
		#
		if ($ValidChars{$this_char}) {
		    $InputMSA[$k][$j-1] = $InputMSA[$k][$j];
		    $InputMSA[$k][$j]   = '-';
		}

	    }

	}

    }


    # Second, the reverse pass.  We certainly can't shift anyone
    # in column input_len-1 rightwards.
    #
    for ($j=$input_len-2; $j>=0; $j--) {


	# We don't fool around with splice sites.  Similarly, we
	# can't do a left shift if the last position was a splice
	# site...
	#
	next if ($InputMSA[0][$j] eq '*' || $InputMSA[0][$j+1] eq '*');


	# Is there a consensus character?  Is there disagreement?
	#
	my $sample_char  = 0;
	my $disagreement = 0;
	my @OpenRight;
	my @ClosedRight;
	my %ValidChars;
	for ($i=0; $i<$num_seqs; $i++) {

	    # Are we open to the idea of a rightwards shift?
	    #
	    if ($InputMSA[$i][$j+1] eq '-') {
		push(@OpenRight,$i);
	    } else {
		push(@ClosedRight,$i);
		$ValidChars{uc($InputMSA[$i][$j-1])} = 1;
	    }

	    # Don't distract me with your dumb old gaps!
	    #
	    next if ($InputMSA[$i][$j] eq '-');

	    # Normalize to upper case
	    #
	    my $this_char = uc($InputMSA[$i][$j]);

	    # Is this the first non-gap we've seen?
	    #
	    if (!$sample_char) {
		$sample_char = $this_char;

	    } elsif ($this_char ne $sample_char) { # COULD IT BE?!
		$disagreement = 1;

	    }

	}

	# Did we find any disagreements within this column?  Is fun on the table?
	#
	if ($disagreement && scalar(@OpenRight)) { # Here comes the fun!

	    # For each character next to a gap, if we have a match scoot on over!
	    #
	    foreach $k (@OpenRight) {

		my $this_char = uc($InputMSA[$k][$j]);

		next if ($this_char eq '-');

		# Shift right!
		#
		if ($ValidChars{$this_char}) {
		    $InputMSA[$k][$j+1] = $InputMSA[$k][$j];
		    $InputMSA[$k][$j]   = '-';
		}

	    }

	}

    }


    # Because we might have introduced some all-gap columns during this
    # first cleanup pass, we need to do some basic cleanup work.  This
    # will involve knowing where splice locations are.
    #
    my @FinalMSA;
    my @SpliceSites;
    my $final_len = 0;
    for ($j=0; $j<$input_len; $j++) {

	# If this is a splice-site column we'll want to know.  This also
	# allows us to skip the other checking work...
	#
	if ($InputMSA[0][$j] eq '*') {

	    # Record that we saw a splice site
	    #
	    push(@SpliceSites,$final_len);

	    # Give everyone an asterisk!
	    #
	    for ($i=0; $i<$num_seqs; $i++) { $FinalMSA[$i][$final_len] = '*'; }

	    # Definitely don't want to overwrite this friend!
	    #
	    $final_len++;

	} else {

	    # Copy over to the FinalMSA.  If it's an all-gap column, we'll
	    # just overwrite it.
	    #
	    my $non_gap_observed = 0;
	    for ($i=0; $i<$num_seqs; $i++) {
		$FinalMSA[$i][$final_len] = $InputMSA[$i][$j];
		$non_gap_observed = 1 if ($InputMSA[$i][$j] ne '-');
	    }
	    $final_len++ if ($non_gap_observed);

	}

    }


    # Copy the 'Final' MSA back over to InputMSA.
    #
    @InputMSA  = @FinalMSA;
    $input_len = $final_len;


    # Reset the 'Final' MSA.
    #
    @FinalMSA  = ();
    $final_len = 0;




    ################
    #              #
    #   PART TWO   #
    #              #
    ################


    #
    #  The aim of the second part of the cleanup is to look
    #  for any 'jigsaw' regions of the MSA, where weird gapping
    #  around a splice site may have caused issues.
    #


    # First, check around each of the splice sites to see
    # if there are any places where it makes sense to flip
    # some peptides around
    #
    for ($k=0; $k<@SpliceSites; $k++) {


	my $site_pos = $SpliceSites[$k];


	# If we need to be careful about edges, adjust how far out
	# we're allowed to go.  By default we only allow flips up
	# to 9 amino acids
	#
	my $ext_dist = 9;
	if ($site_pos < $ext_dist) {
	    $ext_dist = $site_pos-1;
	}
	if ($site_pos + $ext_dist > $input_len) {
	    $ext_dist = $input_len-($site_pos+1);
	}


	# This shouldn't happen, but just in case...
	#
	next if ($ext_dist <= 0);


	# Check those gaps
	#
	my @GapDirs = ();
	my @GapStrs = ();
	for ($i=0; $i<$num_seqs; $i++) {


	    # If there's a gap on only one side, how far can
	    # we go before we hit the sequence on that side?
	    # We don't consider flips larger than 9 amino acids.
	    #
	    my $gap_dir = 0;
	    my $gap_str = '';


	    # 1. Gap on the left, string on the right
	    #
	    if ($InputMSA[$i][$site_pos-1] =~ /\-/ && $InputMSA[$i][$site_pos+1] =~ /\w/) {

		$gap_dir = -1;

		$j=1;
		while ($j<$ext_dist && $InputMSA[$i][$site_pos-$j] =~ /\-/ && $InputMSA[$i][$site_pos+$j] =~ /\w/) {
		    last if ($InputMSA[$i][$site_pos+$j] eq '*');
		    $gap_str = $gap_str.$InputMSA[$i][$site_pos+$j];
		    $gap_dir--;
		    $j++;
		}

		# TOO LONG!
		if ($j == $ext_dist) { $gap_dir = 0; }

	    }
	    #
	    # 2. Gap on the right, string on the left
	    #
	    elsif ($InputMSA[$i][$site_pos-1] =~ /\w/ && $InputMSA[$i][$site_pos+1] =~ /\-/) {

		$gap_dir = 1;

		$j=1;
		while ($j<$ext_dist && $InputMSA[$i][$site_pos-$j] =~ /\w/ && $InputMSA[$i][$site_pos+$j] =~ /\-/) {
		    last if ($InputMSA[$i][$site_pos-$j] eq '*');
		    $gap_str = $gap_str.$InputMSA[$i][$site_pos-$j];
		    $gap_dir++;
		    $j++;
		}

		# TOO LONG!
		if ($j == $ext_dist) { $gap_dir = 0; }

	    }


	    # Regardless of what we saw, add it to the list
	    #
	    push(@GapDirs,$gap_dir);
	    push(@GapStrs,$gap_str);

	}

	# Run through our hits, figuring out who we can swap around
	#
	for ($i=0; $i<$num_seqs; $i++) {

	    next if ($GapDirs[$i] == 0);

	    my @Flippers  = ();
	    my $left_rep  = 0;
	    my $right_rep = 0;

	    if ($GapDirs[$i] < 0) { $right_rep = 1; }
	    else                  { $left_rep  = 1; }

	    for ($j=$i+1; $j<$num_seqs; $j++) {

		# If we have a matching sequence, add it to the flippers
		#
		if ($GapDirs[$j] && $GapStrs[$j] eq $GapStrs[$i]) {

		    push(@Flippers,$j);

		    # Repping a new side?
		    #
		    if ($GapDirs[$i] < 0) { $right_rep = 1; }
		    else                  { $left_rep  = 1; }

		}

	    }

	    # First off, any matches at all?
	    #
	    if (@Flippers) {

		# Cool!  Now we can go through and scoot anybody
		# on the right side over to the left (arbitrary
		# directional choice), assuming both sides are
		# represented.
		#
		if ($left_rep && $right_rep) {

		    # The actual string (as a char array)
		    #
		    my @FlipStr = split('',$GapStrs[$i]);

		    # Make it be so!
		    #
		    for ($j=0; $j<@Flippers; $j++) {
			if ($GapDirs[$Flippers[$j]] < 0) {
			    my $seq        = $Flippers[$j];
			    my $char_start = $site_pos-@FlipStr;
			    my $gap_start  = $site_pos+1;
			    for (my $pos=0; $pos<@FlipStr; $pos++) {
				$InputMSA[$seq][$char_start+$pos] = $FlipStr[$pos];
				$InputMSA[$seq][$gap_start +$pos] = '-';
			    }
			}
		    }

		}

		# These seqs should now be off-limits for flipping
		#
		for ($j=0; $j<@Flippers; $j++) { $GapDirs[$Flippers[$j]] = 0; }

	    }

	}

    }


    # Rad!  All of that flipping may have created some columns that
    # are all gaps, though, so we'll need to clean those up.
    #
    for ($j=0; $j<$input_len; $j++) {

	# If this is a splice-site column we'll want to know.  This also
	# allows us to skip the other checking work...
	#
	if ($InputMSA[0][$j] eq '*') {

	    # Record that we saw a splice site
	    #
	    push(@SpliceSites,$final_len);

	    # Give everyone an asterisk!
	    #
	    for ($i=0; $i<$num_seqs; $i++) { $FinalMSA[$i][$final_len] = '*'; }

	    # Definitely don't want to overwrite this friend!
	    #
	    $final_len++;

	} else {

	    # Copy over to the FinalMSA.  If it's an all-gap column, we'll
	    # just overwrite it.
	    #
	    my $non_gap_observed = 0;
	    for ($i=0; $i<$num_seqs; $i++) {
		$FinalMSA[$i][$final_len] = $InputMSA[$i][$j];
		$non_gap_observed = 1 if ($InputMSA[$i][$j] ne '-');
	    }
	    $final_len++ if ($non_gap_observed);

	}

    }


    # Return our cleaned-up MSA
    #
    return(\@FinalMSA,$final_len);

}








#####################################################
##################                 ##################
##################   END OF FILE   ##################
##################                 ##################
#####################################################















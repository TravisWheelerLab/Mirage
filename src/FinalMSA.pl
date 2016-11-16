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

    # In case we have an extra "chromosome" field, strip it out
    my $header = $line;
    if ($header =~ /([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|\S+$/) {
	$header = $1.'|'.$2.'|'.$3.'|'.$4;
    }
    if ($header !~ /^>GN\:/) {
	$header =~ s/\>//g;
	$header = '>GN:'.$header;
    }

    $Headers[$numSeqs] = $header;

    $i = 0;
    $line = <$inf>;
    while (!eof($inf) && $line !~ /^\>/) {
	$line =~ s/\n|\r//g;
	if ($line) {
	    my @seqline = split('',$line);
	    foreach $j (0..@seqline-1) {
		$BigMSA[$numSeqs][$i] = $seqline[$j];
		$i++;
	    }
	}
	$line  = <$inf>;
    }
    $MSALength = $i;

    $numSeqs++;

}
close($inf);


# Remove intron gaps from the MSA, correcting any ugly shifts
my $MSAref;
($MSAref,$MSALength) = RemoveIntronGaps(\@BigMSA,$MSALength,$numSeqs);
@BigMSA = @{$MSAref};


# Print out our final MSA
open(my $outf,'>',$ARGV[1]);
foreach $i (0..$numSeqs-1) {
    $j = 0;
    print $outf "$Headers[$i]\n";
    while ($j < $MSALength) {
	print $outf "$BigMSA[$i][$j]";
	$j++;
	if ($j % 50 == 0) {
	    print $outf "\n";	   
	}
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
    my ($i,$j,$k);

    my $MSAref = shift;
    my @BigMSA = @{$MSAref};

    my $MSALength = shift;
    my $numSeqs   = shift;

    # Don't mess with the starter intron or ender intron
    $j = 1;
    while ($j < $MSALength-1) {
	
	# Does this position indicate an intron?
	my $intron = 0;
	foreach $i (0..$numSeqs-1) {
	    if ($BigMSA[$i][$j] eq '*') {
		foreach $k (0..$numSeqs-1) {
		    $BigMSA[$k][$j] = '-';
		}
		$intron = 1;
		last;
	    }
	}

	# Indeed it did!  How grand!
	if ($intron) {
	    
	    # Now we check what the tightest gap around the intron boundary
	    # is.
	    my $lowest  = $MSALength;
	    my $highest = 0;
	    foreach $i (0..$numSeqs-1) {

		my $low = $j-1;
		while ($low >= MAX(0,$j-5) && $BigMSA[$i][$low] eq '-') {
		    $low--;
		}
		$low++;

		my $high = $j+1;
		while ($high < MIN($MSALength,$j+5) && $BigMSA[$i][$high] eq '-') {
		    $high++;
		}
		$high--;

		$lowest  = $low  if ($low < $lowest);
		$highest = $high if ($high > $highest);

	    }


	    # In case we don't have a gap on one side or the other (but there is some
	    # gapping) we extend our search just a little
	    if ($lowest != $highest) {
		if ($lowest == $j) {
		    $lowest -= ($highest-$j)+1;
		    $lowest  = 0 if ($lowest < 0);
		} elsif ($highest == $j) {
		    $highest += ($j-$lowest)+1;
		    $highest  = $MSALength-1 if ($highest >= $MSALength);
		}
	    }


	    # So long as everybody has at least some gappage going on,
	    # we should be curious about what's happening here
	    my $gapsize = ($highest-$lowest)+1;
	    if ($gapsize > 1) {

		# Grab all non-gap characters within gapsize of the intron
		# marker, on either side.
		my $pos;
		my @MicroSegs;
		my $microLen = 0;
		foreach $i (0..$numSeqs-1) {
		    $k = 0;
		    foreach $pos ($lowest..$highest) {
			if ($BigMSA[$i][$pos] ne '*' && $BigMSA[$i][$pos] ne '-') {
			    $MicroSegs[$i][$k] = $BigMSA[$i][$pos];
			    $k++;
			}
		    }
		    $microLen = $k if ($k > $microLen);
		}

		# Add in any gap characters to round out
		foreach $i (0..$numSeqs-1) {
		    foreach $k (0..$microLen-1) {
			if (!$MicroSegs[$i][$k]) {
			    $MicroSegs[$i][$k] = '-';
			}
		    }
		}

		# Do these sequences share at least 50% identity?
		my $acceptable = int($gapsize/2);
		my $k = $highest-$lowest;
		foreach $pos (0..$microLen-1) {
		    foreach $i (1..$numSeqs-1) {
			if ($MicroSegs[$i][$pos] ne $MicroSegs[0][$pos]) {
			    $acceptable--;
			    last;
			}
		    }
		    last if ($acceptable <= 0);
		}

		# Do we accept the possibility that these should be aligned?
		if ($acceptable) {

		    foreach $i (0..$numSeqs-1) {
			$pos = 0;
			while ($pos < $microLen) {
			    $BigMSA[$i][$lowest+$pos] = $MicroSegs[$i][$pos];
			    $pos++;
			}
			$pos = $pos+$lowest;
			while ($pos <= $highest) {
			    $BigMSA[$i][$pos] = '-';
			    $pos++;
			}
		    }
		    $j = $highest+1;
		    
		} else { # We didn't move anything
		    $j++;
		}

	    } else { # No gap around the intron	
		$j++;
	    }

	} else { # Not an intron
	    $j++;
	}

    }

    # Finally, we remove all all-gap columns (the only surviving introns should be
    # the first and last, which we just gloss over).
    $k = 0;
    foreach $j (1..$MSALength-2) {

	# Is this an all-gap column?
	my $allGaps = 1;
	foreach $i (0..$numSeqs-1) {
	    if ($BigMSA[$i][$j] ne '-') {
		$allGaps = 0;
		last;
	    }
	}

	# If this was an all-gap column, we don't record it.  Otherwise, we do.
	if (!$allGaps) {
	    foreach $i (0..$numSeqs-1) {
		$BigMSA[$i][$k] = $BigMSA[$i][$j];
	    }
	    $k++;
	}
	
    }
    
    return (\@BigMSA,$k);

}

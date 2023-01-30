#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


if (@ARGV != 1) {
    die "\n  USAGE: Compare-Spaln-Outputs.pl [expected.out] [observed.out]\n\n";
}


my $ObservedFile = OpenInputFile($ARGV[1]);
my %ObservedHits;
while (my $line = <$ObservedFile>) {
    $line =~ s/\n|\r//g;
    if ($line =~ /^\s*Score\s+\=/) {
	$line =~ s/\s//g; # Just in case anything funny can happen with spacing...
	$ObservedHits{$line} = 1;
    }
}
close($ObservedFile);


my $ExpectedFile = OpenInputFile($ARGV[0]);
while (my $line = <$ExpectedFile>) {
    $line =~ s/\n|\r//g;
    if ($line =~ /^\s*Score\s+\=/) {
	my $original_line = $line;
	$line =~ s/\s//g;
	if (!$ObservedHits{$line}) {
	    close($ExpectedFile);
	    die "\n  ERROR: Expected Spaln hit '$original_line' not found in observed test outputs\n\n";
	}
    }
}
close($ExpectedFile);

1;



sub OpenInputFile
{
    my $fname = shift;
    if (!(-e $fname)) {
	die "\n  ERROR:  Failed to locate file '$fname'\n\n";
    }
    open(my $File,'<',$fname)
	|| die "\n  ERROR:  Failed to open input file '$fname'\n\n";
    return $File;
}











# EOF

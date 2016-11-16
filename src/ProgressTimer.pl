#!/usr/bin/env perl

use warnings;
use strict;


if (@ARGV != 3) { die "\n  USAGE:  ./ProgressTimer  <filebase>  <numthreads>  <basecount>\n\n"; }


my $progressbase = $ARGV[0];
my $numthreads   = int($ARGV[1]);
my $overall_prog = int($ARGV[2]);


# Iterate over prog. files, getting a full progress count
for (my $i = 1; $i < $numthreads; $i++) {
    
    # At the very start, some threads might not quite
    # beat 0 to the punch
    if (-e $progressbase.$i) {
	open(my $progfile,'<',$progressbase.$i);
	my $pcount = readline($progfile);
	close($progfile);

	# REALLY special case - file created but nothing
	# to pull out
	if ($pcount && $pcount =~ s/\n|\r//g) {
	    $overall_prog += int($pcount);
	}
    }
    

}

# Print out the overall progress
print "$overall_prog\n";

1;

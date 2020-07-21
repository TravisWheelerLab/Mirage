#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

# We'll definitely want our directory-checking / file IO helpers
sub GetBM { my $lib = $0; $lib =~ s/Quilter2.pl$//; return $lib; }
use lib GetBM();
use BureaucracyMirage;

# So every function doesn't need to loop up for length...
sub PrintProgress;
sub ClearProgress;

# Many functions for specific outs-put
sub ProgressMirageInit;
sub ProgressMirageQuilter;
sub ProgressMirageMultiMSA;
sub ProgressMirageMultiSeqNW;
sub ProgressMirageCleanup;

# Before we get going, let's set a standardized character limit for status messages
my $max_line_len = 80;




###############################################################
#
sub PrintProgress
{
    my $str = shift;
    $str = $str.' ' while (length($str) < $max_line_len);
    print "$str\r";
}
sub ClearProgress { PrintProgress(' '); }



###############################################################
#
sub ProgressMirageInit
{
    # What part of the initialization are we immersed in?
    my $part = shift;
    my $status;

    if ($part eq 'db-speciation') {
	$status = "Dividing protein database according to species";
    }
    
    PrintProgress($status);
}






###############################################################
#
sub ProgressMirageQuilter
{
    # What part of Quilter are we working on?
    my $part = shift;
    my $threadID = shift;
    my $species = shift;
    my $progress_dirname = shift;

    my $status = "  Quilter ($species):";
    if ($part eq 'parsing-gtf') {
	$status = $status."Loading GTF data from file";
    }

    # Only the master gets to report progress (lucky duck!)
    PrintProgress($status) if (!$threadID);
    
}






###############################################################
#
sub ProgressMirageMultiMSA
{

}






###############################################################
#
sub ProgressMirageMultiSeqNW
{

}






###############################################################
#
sub MirageCleanup
{

}





1; # EOF

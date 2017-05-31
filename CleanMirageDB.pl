#!/usr/bin/env perl
#
# CleanMirageDB.pl - Alex Nord - 2016
#
# USAGE: $ perl CleanMirageDB.pl <database>
#
# ABOUT: This program quickly runs through a FASTA-formatted database
#        and makes any changes necessary for mirage to run functionally
#
#
use warnings;
use strict;


# Assist the user
if (@ARGV != 1) {
    print "\n";
    print "  USAGE: perl CleanMirageDB.pl <database>\n\n";
    print "  ABOUT: This program will takes a FASTA-formatted database\n";
    print "         as its input and correct any instances of illegal\n";
    print "         characters (#,&,/) or sequence names without actual\n";
    print "         sequence.  If any issues are found, the original\n";
    print "         sequence names are recorded in a 'Corrections' file\n";
    print "         and any 'ghost' sequences are recorded in a 'NoSequence'\n";
    print "         file.\n";
    print "\n";
    print "         It also looks for the particular structure of fields\n";
    print "         required by mirage, and any sequences that fail to\n";
    print "         conform to this format are removed from the resulting\n";
    print "         database (with names recorded in a 'Removed' file).\n";    
    print "\n";
    print "         If no problems are discovered, the original database\n";
    print "         may be provided as mirage input.  Otherwise, a 'Clean'\n";
    print "         database is generated.\n";
    die   "\n";
}


# Grab the original database name
my $dbfilename = $ARGV[0];
open(my$database,'<',$dbfilename) || die "\n  ERROR: Failed to open input database '$ARGV[0]'\n\n";


# Strip off the '.fa' from the original database name to
# figure out what we'll name our output files.
# Now with '.fasta' support!
my $cleandbname = $dbfilename;
$cleandbname =~ s/\.fa$|\.fasta$//;


# Name all of the files that we might be generating
my $corrections  = $cleandbname;
my $ghostseqfile = $cleandbname;
my $removedfile  = $cleandbname;
$cleandbname  = $cleandbname.'.Clean.fa';
$corrections  = $corrections.'.Corrections.out';
$ghostseqfile = $ghostseqfile.'.NoSequence.out';
$removedfile  = $removedfile.'.Removed.out';


# Open up all of our output files
open(my $clean,'>',$cleandbname)   || die "\n  ERROR: Failed to open output database '$cleandbname'\n\n";
open(my $fixes,'>',$corrections)   || die "\n  ERROR: Failed to open file '$corrections'\n\n";
open(my $ghosts,'>',$ghostseqfile) || die "\n  ERROR: Failed to open file '$ghostseqfile'\n\n";
open(my $gone,'>',$removedfile)    || die "\n  ERROR: Failed to open file '$removedfile'\n\n";


# Record some data about how much we're changing
my $num_fixes  = 0;
my $num_gone   = 0;
my $num_ghosts = 0;
my $comments   = 0;


# If there are any comments starting it off we cut them out
my $line = <$database>;
while ($line && $line !~ /^>/) {
    $comments++;
    $line = <$database>;
}


# If we hit the end of the file, this is bad
if (eof($database)) {
    close($database);
    close($clean);    
    close($fixes);
    close($ghosts);
    close($gone);
    system("rm $clean") if (-e $clean);
    die "\n  ERROR: Database '$dbfilename' does not appear to be FASTA-formatted\n\n";    
}


# Cool! let's get to work
while (!eof($database)) {

    # Grab the next header line.
    $line =~ s/\n|\r//g;
    my $originalhead = $line;
    my $fixedhead    = $line;
    $fixedhead =~ s/\s|\#|\&|\/|\'|\"//g;

    # Is this a non-functional sequence header?
    my $removed = 0;
    if ($fixedhead !~ /(\S+)\|(\S+)\|(\S+)\|(\S+)\|(\S+)/ || $fixedhead !~ /\|[^\|\s]+\s*$/) {
	$removed = 1;
    }
    
    # How many actual lines were there?
    my $seqlines = 0;
    my $sequence = $fixedhead."\n";

    # Scan over the sequence content
    $line = <$database>;
    while ($line && $line !~ /^>/) {
	$line =~ s/\n|\r//g;
	if ($line) {
	    $seqlines++;
	    $sequence = $sequence.$line."\n";
	}
	$line = <$database>;
    }


    # Was this an actual sequence?
    if ($removed) {
	print $gone "$originalhead\n";
	$num_gone++;
    } else {
	if ($seqlines) {
	    print $clean "$sequence";
	    if ($fixedhead ne $originalhead) {
		$num_fixes++;
		print $fixes "$originalhead : $fixedhead\n";
	    }
	} else {
	    $num_ghosts++;
	    print $ghosts "$originalhead\n";
	}
    }

}


# Close up shop
close($database);
close($clean);
close($fixes);
close($gone);


# Let the user know what we discovered
if ($num_ghosts || $num_fixes || $num_gone || $comments) {

    print "\n";

    if ($comments) {
	print "  Removed comment $comments lines from start of file\n";
    }

    if ($num_fixes) {
	print "  Number of sequence names adjusted: $num_fixes (listed in '$corrections')\n";
    } else {
	system("rm $corrections")  if (-e $corrections);
    }

    if ($num_ghosts) {
	print "  Number of names without sequences: $num_ghosts (listed in '$ghostseqfile')\n";
    } else {	
	system("rm $ghostseqfile") if (-e $ghostseqfile);
    }

    if ($num_gone) {
	print "  Number of sequences removed:       $num_gone (listed in '$removedfile')";
    } else {
	system("rm $removedfile") if (-e $removedfile);
    }

    print "\n";
    print "  Clean database filename: '$cleandbname'\n\n";

} else {

    print "\n  Original dataset '$dbfilename' is clean\n\n";

    system("rm $corrections")  if (-e $corrections);
    system("rm $ghostseqfile") if (-e $ghostseqfile);
    system("rm $cleandbname")  if (-e $cleandbname);
    system("rm $removedfile")  if (-e $removedfile);

}


# All done!
1;

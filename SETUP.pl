#!/usr/bin/env perl
#
# Setup.pl - Alex Nord - 2016 
#
# ABOUT: This script is used to set up the directory for running 
#        mirage.pl.  It verifies that all required programs are
#        accessible and correctly compiled.  After running SETUP.pl,
#        the user should only need to add the easel 'miniapps'
#        directory and the spaln 'src' directory to her or his PATH.
#
use warnings;
use strict;
use Cwd;


# The name of the spaln folder (at the top for ease of adjustment)
my $spalnDir = 'inc/spaln2.2.2';
my $easelDir = 'inc/easel';
my $blatDir  = 'inc/blat';
my $spalnTar = $spalnDir.'.tgz';
my $easelTar = $easelDir.'.tgz';
my $blatTar  = $blatDir.'.tgz';


# Before all else, make sure we have autoconf (most likely thing to
# be missing)
open(my $autoconfcheck,"which autoconf |");
my $line = <$autoconfcheck>;
if (!$line) { 
    print "\n  ERROR: Program 'autoconf' is not installed (required for easel installation)";
    die   "\n         Try 'brew install autoconf'\n\n";
}
close($autoconfcheck);


# Where we are RIGHT NOW
my $startDir = getcwd;
$startDir =~ s/\/$//;


# Check that all files are where we need them
my @srcFiles;
push(@srcFiles,'src/makefile');
push(@srcFiles,'src/Diagonals.c');
push(@srcFiles,'src/Diagonals.h');
push(@srcFiles,'src/FindDiagonals.c');
push(@srcFiles,'src/TransSW.c');
push(@srcFiles,'src/MultiSeqNW.c');
push(@srcFiles,'src/MultiSeqNW.h');
push(@srcFiles,'src/DiagonalSets.pm');
push(@srcFiles,'src/Quilter.pl');
push(@srcFiles,'src/MultiMSA.pl');
push(@srcFiles,'src/FinalMSA.pl');
push(@srcFiles,'mirage.pl');
foreach my $file (@srcFiles) {
    if (!(-e $file)) { die "\n  Failed to locate critical file '$file'\n\n"; }
}


# Run the mirage makefile
print "\n  Compiling mirage files\n\n";
chdir("src") || die "\n  Failed to enter directory 'src'\n\n";
if (system("make")) { die "\n  Failed to compile C source files\n\n"; }


# For the other programs included in the package, check to see if they're
# already available before we go through the trouble of setting them up


# Check out spaln
open(my $spalncheck,"which spaln |");
my $spaln = <$spalncheck>;
close($spalncheck);
if (!$spaln) {

    # We need to unpack spaln
    chdir($startDir);
    if (system("tar -xzf $spalnTar -C inc/")) { die "\n  Failed to expand file '$spalnTar'\n\n"; }

    # Configure, make, and install spaln
    chdir("$startDir/$spalnDir/src") || die "\n  Failed to enter directory '$spalnDir/src'\n\n";
    print "\n  Configuring and compiling spaln\n\n";
    if (system("./configure")) { die "\n  Failed to configure 'spaln'\n\n"; }
    if (system("make")) { die "\n  Failed to make 'spaln'\n\n"; } 
    #if (system("make install")) { die "\n  Failed to install 'spaln'\n\n"; }

}    


# Now we need to get easel unpacked and setup, too
open(my $eslsfetchcheck,"which esl-sfetch |");
my $eslsfetch = <$eslsfetchcheck>;
close($eslsfetchcheck);
open(my $eslseqstatcheck,"which esl-seqstat |");
my $eslseqstat = <$eslseqstatcheck>;
close($eslseqstatcheck);
if (!$eslsfetch || !$eslseqstat) {

    chdir($startDir);
    print "\n  Configuring and compiling easel\n\n";
    if (system("tar -xzf $easelTar -C inc/")) { die "\n  Failed to expand '$easelTar'\n\n"; }
    
    # Configure
    chdir($easelDir);
    if (system "autoconf") { die "\n  Failed to configure easel library -- is autoconf installed?\n\n"; }
    if (system "./configure") { die "\n  Failed to configure easel library\n\n"; }

    # Make and make install
    if (system("make")) { die "\n  Failed to compile easel library\n\n"; }
    if (system("make check")) { die "\n  Failure during easel make check\n\n"; }
    #if (system("make install")) { die "\n  Failed to install easel library\n\n"; }

}    


# Finally, BLAT
chdir($startDir);
if (!(-e './inc/blat/blat')) {
    
    if (system("tar -xzf $blatTar -C inc/")) { die "\n  Failed to expand '$blatTar'\n\n"; }    
    chdir($blatDir);
    
    # What sort of OS are we on? (opt.s: Linux, Darwin, other)
    print "\n  Guessing correct BLAT... ";
    open(my $uname, 'uname -a |');
    my $sysline = <$uname>;
    my $BLATlink;
    if (uc($sysline) =~ /^LINUX /) {
	print "Linux\n\n";
	$BLATlink = "ln -s blat.linux.x86_64 blat";
    } elsif (uc($sysline) =~ /^DARWIN /) {
	print "OSX\n\n";
	$BLATlink = "ln -s blat.macOSX.x86_64 blat";
    } else {
	print "OSX i386 (guessing)\n\n";
	$BLATlink = "ln -s blat.macOSX.i386 blat";
    }
    close($uname);

    # Make a symbolic link
    if (system($BLATlink)) { die "\n  Failed to create symbolic link to BLAT binary\n\n"; }
    
}

    
# If everything went well, we can go ahead an clear out the
# spaln and easel tarballs
chdir($startDir);
system("rm $spalnTar") if (-e $spalnTar);
system("rm $easelTar") if (-e $easelTar);
system("rm $blatTar")  if (-e $blatTar);


# Now the only thing the user REALLY needs to do is make
# sure that spaln is on their PATH.
print "\n\n  Setup completed successfully!";
print "\n\n  To complete installation add spaln and easel to your PATH\n";
print "  (see README for assistance)\n\n";


# Happy mirage-ing!
1;


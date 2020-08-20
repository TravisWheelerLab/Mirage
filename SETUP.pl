#!/usr/bin/env perl
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


my $argcheck=0;
my $force=0;
while ($argcheck < scalar(@ARGV)) {
    if (lc($ARGV[$argcheck]) eq '-force') {
	$force=1;
    } else {
	print "\n  Unrecognized option '$ARGV[$argcheck]' ignored.\n\n";
    }
    $argcheck++;
}


# If there's already a link to mirage2, we'll kill it, just in case
# it doesn't look like it's safe to try miragery
if (-e 'mirage2') { system("rm \"mirage2\""); }


# The name of the spaln folder (at the top for ease of adjustment)
my $spalnDir = 'inc/spaln2.3.3';
my $easelDir = 'inc/easel';
my $blatDir  = 'inc/blat';
my $spalnTar = $spalnDir.'.tgz';
my $easelTar = $easelDir.'.tgz';
my $blatTar  = $blatDir.'.tgz';


# Before all else, make sure we have autoconf (most likely thing to
# be missing)
open(my $autoconfcheck,"which autoconf |");
my $line = <$autoconfcheck>;
if (!$line) { die "\n  ERROR: Program 'autoconf' is not installed (required for easel installation)\n         Try 'brew install autoconf'\n\n"; }
close($autoconfcheck);


# Where we are RIGHT NOW
my $startDir = getcwd;
$startDir =~ s/\/$//;


# Check that all files are where we need them
my @srcFiles;
push(@srcFiles,'src/makefile');
push(@srcFiles,'src/BasicBio.c');
push(@srcFiles,'src/BasicBio.h');
push(@srcFiles,'src/BureaucracyMirage.pm');
push(@srcFiles,'src/ExonWeaver.c');
push(@srcFiles,'src/FastMap2.c');
push(@srcFiles,'src/FinalMSA.pl');
push(@srcFiles,'src/Mirage2.pl');
push(@srcFiles,'src/MultiSeqNW.c');
push(@srcFiles,'src/MultiSeqNW.h');
push(@srcFiles,'src/MapsToMSAs.pl');
push(@srcFiles,'src/Quilter2.pl');
push(@srcFiles,'src/run_mirage2.sh');
foreach my $file (@srcFiles) {
    if (!(-e $file)) { die "\n  Failed to locate critical file '$file'\n\n"; }
}


# Run the mirage makefile
print "\n  Compiling mirage files\n\n";
chdir("src") || die "\n  Failed to enter directory 'src'\n\n";
if (system("make")) { die "\n  Failed to compile C source files\n\n"; }


# For the other programs included in the package, check to see if they're
# already available before we go through the trouble of setting them up


# We may need to unpack spaln
chdir($startDir);

# If forcing this might not be necessary
#
if (-e $spalnTar) {

    if (system("tar -xzf $spalnTar -C inc/")) { die "\n  Failed to expand file '$spalnTar'\n\n"; }

    # Configure and make SPALN
    chdir("$startDir/$spalnDir/src") || die "\n  Failed to enter directory '$spalnDir/src'\n\n";
    print "\n  Configuring and compiling spaln\n\n";
    if (system("./configure")) { die "\n  Failed to configure 'spaln'\n\n"; }
    if (system("make")) { die "\n  Failed to make 'spaln'\n\n"; } 

}    

# Now we need to get easel unpacked and setup, too
chdir($startDir);
print "\n  Configuring and compiling easel\n\n";
if (-e $easelTar) {

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


# Unpack BLAT
chdir($startDir);
if (-e $blatTar) {
    if (system("tar -xzf $blatTar -C inc/")) { die "\n  Failed to expand '$blatTar'\n\n"; }    
}
    
# If everything went well, we can go ahead an clear out the
# spaln and easel tarballs
system("rm $spalnTar") if (-e $spalnTar);
system("rm $easelTar") if (-e $easelTar);
system("rm $blatTar")  if (-e $blatTar);

# Create a symbolic link to the mirage-running shell script
my $MirageLink = "ln -s src/run_mirage2.sh mirage2";
if (system($MirageLink)) { die "\n  Failed to create symbolic link to src/run_mirage2.sh\n\n"; }

# Now the only thing the user REALLY needs to do is make
# sure that spaln is on their PATH.
print "\n\n  Setup completed successfully!";
print "\n\n  To complete installation add mirage to your PATH\n";
print "  (see README for assistance)\n\n";


# Happy mirage-ing!
1;


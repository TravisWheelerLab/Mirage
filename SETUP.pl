#!/usr/bin/env perl
#
# ABOUT: This script is used to set up the directory for running 
#        mirage.pl.  It verifies that all required programs are
#        accessible and correctly compiled.
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
my $hsiDir   = 'dependencies/hsi-1.0.0';
my $spalnDir = 'dependencies/spaln2.3.3';
my $blatDir  = 'dependencies/blat';
my $tbnDir   = 'dependencies/tblastn';
my $hsiTar   = $hsiDir.'.tgz';
my $spalnTar = $spalnDir.'.tgz';
my $blatTar  = $blatDir.'.tgz';
my $tbnTar   = $tbnDir.'.tgz';


# Where we are RIGHT NOW
my $startDir = getcwd;
$startDir =~ s/\/$//;


# Check that all files are where we need them
my @srcFiles;
push(@srcFiles,'Makefile');
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
push(@srcFiles,'src/Oasis.pl');
push(@srcFiles,'src/Quilter2.pl');
push(@srcFiles,'src/run_mirage2.sh');
push(@srcFiles,'src/run_oasis.sh');
foreach my $file (@srcFiles) {
    if (!(-e $file)) { die "\n  Failed to locate critical file '$file'\n\n"; }
}


# Run the mirage makefile
print "\n  Compiling mirage files\n\n";
if (system("make")) { die "\n  Failed to compile C source files\n\n"; }


# For the other programs included in the package, check to see if they're
# already available before we go through the trouble of setting them up


# Now we need to get easel unpacked and setup, too
chdir($startDir);
if (-e $hsiTar) {

    if (system("tar -xzf $hsiTar -C dependencies/")) { die "\n  Failed to expand '$hsiTar'\n\n"; }

    # Make and make install
    chdir($hsiDir) || die "\n  Failed to enter directory '$hsiDir'\n\n";
    print "\n  Compiling hsi tools\n\n";
    if (system("cmake .")) { die "\n  Failed to cmake hsi library\n\n"; }
    if (system("make")) { die "\n  Failed to make hsi library\n\n"; }

}


# Unpack Spaln, if necessary
chdir($startDir);
if (-e $spalnTar) {

    if (system("tar -xzf $spalnTar -C dependencies/")) { die "\n  Failed to expand file '$spalnTar'\n\n"; }

    # Configure and make SPALN
    chdir("$startDir/$spalnDir/src") || die "\n  Failed to enter directory '$spalnDir/src'\n\n";
    print "\n  Configuring and compiling spaln\n\n";
    if (system("./configure")) { die "\n  Failed to configure 'spaln'\n\n"; }
    if (system("make")) { die "\n  Failed to make 'spaln'\n\n"; } 

}    


# Unpack BLAT
chdir($startDir);
if (-e $blatTar) {
    if (system("tar -xzf $blatTar -C dependencies/")) { die "\n  Failed to expand '$blatTar'\n\n"; }    
}


# Unpack tblastn
if (-e $tbnTar) {
    if (system("tar -xzf $tbnTar -C dependencies/")) { die "\n  Failed to expand '$tbnTar'\n\n"; }
}

    
# If everything went well, we can go ahead an clear out the
# spaln and easel tarballs
system("rm $hsiTar")   if (-e $hsiTar);
system("rm $spalnTar") if (-e $spalnTar);
system("rm $blatTar")  if (-e $blatTar);
system("rm $tbnTar")   if (-e $tbnTar);

# Create a symbolic link to the mirage-running shell script
my $MirageLink = "ln -s src/run_mirage2.sh mirage2";
if (system($MirageLink)) { die "\n  Failed to create symbolic link to src/run_mirage2.sh\n\n"; }

# Create a symbolic link for running bazaar
# -- Bazaar certainly doesn't need public scrutiny at this stage!
#my $BazaarLink = "ln -s src/run_bazaar.sh bazaar";
#if (system($BazaarLink)) { die "\n  Failed to create symbolic link to src/run_bazaar.sh\n\n"; }

# Create a symbolic link for running oasis
# -- For present purposes, we'll keep Oasis a sort-of-secret
#my $OasisLink = "ln -s src/run_oasis.sh oasis";
#if (system($OasisLink)) { die "\n  Failed to create symbolic link to src/run_oasis.sh\n\n"; }

# Happy mirage-ing!
print "\n  Setup completed successfully!\n\n";
1;


#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


#  FUNCTIONAL SUBROUTINES
#


#  BEUREAUCRATIC SUBROUTINES
#
#  NOTE: Because I think this will be useful outside of Mirage, I'm not going to
#        link it to the src/BureaucracyMirage.pm file, even though that essentially
#        has all of the same functions.
#
sub MIN;
sub MAX;
sub RunSystemCommand;
sub OpenSystemCommand;
sub OpenInputFile;
sub NameOutputFile;
sub OpenOutputFile;
sub OpenFileToAppendTo;
sub DeleteFile;
sub ConfirmDirectory;
sub OpenDirectory;
sub NameDirectory;
sub CreateDirectory;



############
#          #
#  SCRIPT  #
#          #
############



# Because the "default" option is to download a truly insane amount
# of sequence data, we'll require the user to specify that they either
# want to download *everything* or else provide a file listing species
# of interest
if (scalar(@ARGV) != 1){
    print "\n";
    print "  USAGES:\n";
    print "\n";
    print "  % ./DownloadGenomicData.pl [file]   <-- Downloads genome data for species listed in file\n";
    print "  % ./DownloadGenomicData.pl --full   <-- Downloads ALL UCSC genome data (this is a LOT)\n";
    die   "\n";
}


# See if there's a file to parse
my %SpeciesToDownload;
my $download_all = 1;
if (lc($ARGV[0]) !~ /\-full/) {

    $download_all = 0;
    
    my $inf = OpenInputFile($ARGV[0]);
    while (my $line = <$inf>) {

	# Eliminate any leading or trailing whitespace
	$line =~ s/\n|\r//g;
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;

	# Convert any spaces to underscores
	$line =~ s/\s/\_/g;

	$SpeciesToDownload{lc($line)} = 1;

    }
}


# Make an output directory where we can save all our fun files!
my $outdirname = CreateDirectory("Genomic-Data");


# Now it's time to do the terrible work of downloading some GD html
my $ucsc_dir_fname = $outdirname.'ucsc-dir.html';
RunSystemCommand("wget -O $ucsc_dir_fname https://hgdownload.soe.ucsc.edu/downloads.html");

# We'll use a standard filename for temporary html downloads
my $temp_html_fname = $outdirname.'temp.html';

# Scan through the directory file and find the links to all genome / gtf
# subdirectories containing species data we're interested in
my $UCSCf = OpenInputFile($ucsc_dir_fname);
my %SpeciesToGenomes;
my %SpeciesToGTFs;
my %TroubleSpecies;
while (my $line = <$UCSCf>) {

    # Look for the next section break (these had BETTER stay consistent!)
    $line =~ s/\n|\r//g;
    next if ($line !~ /\<\!\-\- (.+) Download/);

    # What would a human being call this species?
    my $plain_species_name = lc($1);
    $plain_species_name =~ s/\s/\_/g;

    # If we have a list of species to look for, see if this dude's made the cut
    unless ($download_all || $SpeciesToDownload{$plain_species_name}) {
	next;
    }

    # The game is afoot!  What's your shorthand name?
    while ($line = <$UCSCf>) {
	last if ($line =~ /\<\/h3\>/);
    }
    $line =~ /\>([^\<]+)\<\/a\>\)\<\/h3\>/;
    my $shorthand_species_name = $1;

    # Assuming everything holds, we don't really need to parse anything else,
    # since there's a standard convention for the URLs containing genome data
    my $genome_dir_fname = 'https://hgdownload.soe.ucsc.edu/goldenPath/'.$shorthand_species_name.'/bigZips/';


    # TIME FOR THE MONEY MOVES!
    
    
    # 1. Identify the genome download link

    # Download the directory
    if ("wget -O $temp_html_fname $genome_dir_fname") {
	my $trouble_key = $plain_species_name.' ('.$shorthand_species_name.')';
	$TroubleSpecies{$trouble_key} = 1;
	next;
    }

    # Luckily, this is (mostly) plain-text rather than html, so we should have an
    # easy time parsing it.
    my $BigZips = OpenInputFile($temp_html_fname);
    my $genome_fname = 0;
    while (my $bz_line = <$BigZips>) {
	if ($bz_line =~ /^\s*(\S+) \- \"Soft-masked\" assembly sequence in one file/) {
	    $genome_fname = $genome_dir_fname.$1;
	    last;
	}
    }
    close($BigZips);

    # If we hit the end of the file and didn't find a genome download link,
    # note this as a "trouble species" and move on
    if (!$genome_fname) {
	my $trouble_key = $plain_species_name.' ('.$shorthand_species_name.')';
	$TroubleSpecies{$trouble_key} = 1;
	next;
    }

    # GENOME AHOY!
    $SpeciesToGenomes{$plain_species_name} = $genome_fname;

    
    # 2. Pull in any gtfs that might be floating around

    # If there's a directory with gtfs it'll be named 'genes' (assuming conventions
    # from time of writing).
    my $gtf_dir_fname = $genome_dir_fname.'genes/';
    next if ("wget -O $temp_html_fname $gtf_dir_fname");

    # Get the names of all the gtf files
    my $Genesf = OpenInputFile($temp_html_fname);
    my $gtf_list_str = '';
    while (my $genes_line = <$Genesf>) {
	if ($genes_line =~ /\<a\s+href\=\"([^\"]+\.gtf\.gz)\"\>/) {
	    my $gtf_fname = $gtf_dir_fname.$1;
	    if ($gtf_list_str) { $gtf_list_str = $gtf_list_str.'|'.$gtf_fname; }
	    else               { $gtf_list_str = $gtf_fname;                   }
	}
    }
    close($Genesf);

    # If there weren't any gtf files in our 'genes' dir, then I have two thoughts:
    # (1.) weird...?, and (2.) whatever
    next if (!$gtf_list_str);

    # GTFs AHOY!
    $SpeciesToGTFs{$plain_species_name} = $gtf_list_str;
    
}
close($UCSCf);


# No need for the main UCSC directory file or the temp html file anymore!
RunSystemCommand("rm $ucsc_dir_fname") if (-e $ucsc_dir_fname);
RunSystemCommand("rm $temp_html_fname") if (-e $temp_html_fname);


# If we didn't get any genomes that isn't ideal...
if (scalar(keys %SpeciesToGenomes) == 0) {
    RunSystemCommand("rmdir $outdirname");
    die "\n  No genomes located...\n\n";
}


# Now we can do our actual BIG downloading!  Hopefully this won't be too ugly...
my $dl_genome_dirname = CreateDirectory($outdirname.'genomes');
my $dl_gtf_dirname = 0;
if (scalar(keys %SpeciesToGTFs)) {
    $dl_gtf_dirname = CreateDirectory($outdirname.'gtfs');
}


# We're going to need to know what directory we're working from
my $pwd_inf = OpenSystemCommand('pwd');
my $pwd = <$pwd_inf>;
close($pwd_inf);
$pwd =~ s/\r|\n//g;
$pwd = $pwd.'/' if ($pwd !~ /\/$/);


# We'll want to position ourselves to really easily build a species guide out of
# what we download, so let's just frickin' do it!
my %SpeciesToSuperNames; # english name + shorthand
my $longest_name_len = 0;
my $longest_genome_len = 0;
my $longest_supername_len = 0;

foreach my $species (sort keys %SpeciesToGenomes) {

    $longest_name_len = length($species) if (length($species) > $longest_name_len);

    my $genome_link = $SpeciesToGenomes{$species};

    # Just in case some are '.fasta.gz'
    $genome_link =~ /\/([^\/]+)(\.[^\.]+)\.gz$/;
    my $species_shorthand = $1;
    my $genome_extension  = $2;

    # What's your "supername"?
    my $supername = $species.' ('.$species_shorthand.')';
    $longest_supername_len = length($supername)
	if (length($supername) > $longest_supername_len);
    $SpeciesToSuperNames{$species} = $supername;

    my $dl_genome_fname = $dl_genome_dirname.$species_shorthand.$genome_extension.'.gz';

    # We'll use rsync, since that's what UCSC recommends
    $genome_link =~ s/https?\:/rsync\:/;
    RunSystemCommand("rsync -a -P $genome_link $dl_genome_fname");

    # Unpack that sucker!
    RunSystemCommand("gunzip $dl_genome_fname");

    # Record the full path to the file
    $dl_genome_fname =~ s/\.gz$//g;
    $dl_genome_fname = $pwd.$dl_genome_fname;
    $SpeciesToGenomes{$species} = $dl_genome_fname;

    $longest_genome_len = length($dl_genome_fname)
	if (length($dl_genome_fname) > $longest_genome_len);

    # If there isn't at least one gtf for this species, then we're moving on!
    next if (!$SpeciesToGTFs{$species});

    my @GTFList = split(/\|/,$SpeciesToGTFs{$species});
    for (my $i=0; $i<scalar(@GTFList); $i++) {

	my $gtf_link = $GTFList[$i];
	
	$gtf_link =~ /\/([^\/]+)$/;
	my $dl_gtf_fname = $dl_gtf_dirname.$1;

	$gtf_link =~ s/https?\:/rsync\:/;

	RunSystemCommand("rsync -a -P $gtf_link $dl_gtf_fname");
	RunSystemCommand("gunzip $dl_gtf_fname");

	$dl_gtf_fname =~ s/\.gz$//;
	$GTFList[$i] = $dl_gtf_fname;
	
    }

    # Now we'll concatenate all of our GTFs into one big honkin' nasty boi!
    my $cat_cmd = '> '.$dl_gtf_dirname.$species_shorthand.'.gtf';
    foreach my $dl_gtf_fname (@GTFList) {
	$cat_cmd = $dl_gtf_fname.' '.$cat_cmd;
    }
    $cat_cmd = 'cat '.$cat_cmd;
    RunSystemCommand($cat_cmd);

    # Record the loooooooong name
    $SpeciesToGTFs{$species} = $pwd.$dl_gtf_dirname.$species_shorthand.'.gtf';

    # We got our big honkin' nasty boi, so let's clean up our lil' beepin' chill gals!
    foreach my $dl_gtf_fname (@GTFList) {
	RunSystemCommand("rm $dl_gtf_fname");
    }
    
}


# WOWEE! How's about that -- we've done all of our downloading!
# Now it's time to build a SpeciesGuide file and let the user know
# how we did overall.


my $SpeciesGuide = OpenOutputFile($outdirname.'Species-Guide');
foreach my $species (sort keys %SpeciesToGenomes) {

    my $genome = $SpeciesToGenomes{$species};

    my $gtf = '-';
    if ($SpeciesToGTFs{$species}) {
	$gtf = $SpeciesToGTFs{$species};
    }

    while (length($species) < $longest_name_len) {
	$species = $species.' ';
    }

    while (length($genome) < $longest_genome_len) {
	$genome = $genome.' ';
    }

    print $SpeciesGuide "$species $genome $gtf\n";

}
close($SpeciesGuide);

my %AllSpecies;
foreach my $species (keys %SpeciesToGenomes) { $AllSpecies{$species} = 1; }
foreach my $species (keys %TroubleSpecies  ) { $AllSpecies{$species} = 1; }

# Oh, let's give even the "trouble species" supernames -- it's a fundamental right!
foreach my $species (keys %TroubleSpecies) {
    $SpeciesToSuperNames{$species} = $species;
    $longest_supername_len = length($species)
	if (length($species) > $longest_supername_len);
}

my $SummaryFile = OpenOutputFile($outdirname.'download-summary.out');

my $species_header = 'Species';
my $species_underl = '-------';
while (length($species_header) < $longest_supername_len) {
    $species_header = $species_header.' ';
    $species_underl = $species_underl.' ';
}

# Three spaces between columns
print $SummaryFile "$species_header   Genome?   GTF(s)?\n";
print $SummaryFile "$species_underl   -------   -------\n";

foreach my $species (sort keys %AllSpecies) {

    my $genome = 'Yes    ';
    my $gtf    = 'Yes    ';
    if ($TroubleSpecies{$species}) {
	$genome = 'No     ';
	$gtf    = 'No     ';
    } elsif (!$SpeciesToGTFs{$species}) {
	$gtf    = 'No     ';
    }

    my $supername = $SpeciesToSuperNames{$species};

    # Let's get pretty with it!
    $supername =~ s/\_/ /g;
    my @NameComponents = split(/ /,$supername);
    $supername = '';
    for (my $i=0; $i<scalar(@NameComponents); $i++) {
	if ($NameComponents[$i] =~ /^([a-z])/) {
	    my $first_letter = uc($1);
	    $NameComponents[$i] =~ s/^[a-z]//;
	    $NameComponents[$i] = $first_letter.$NameComponents[$i];
	}
	$supername = $supername.' ' if ($supername);
	$supername = $supername.$NameComponents[$i];
    }

    while (length($supername) < $longest_supername_len) {
	$supername = $supername.' ';
    }

    print $SummaryFile "$supername   $genome   $gtf\n";
    
}
close($SummaryFile);


# DONE!
print "\n";
print "  Genomic data assembled!  Results in '$outdirname'\n";
print "\n";


1;  ##   END OF SCRIPT   ##





############################
#                          #
#  FUNCTIONAL SUBROUTINES  #
#                          #
############################





#################################################################
#
#  FUNCTION:  
#







###############################
#                             #
#  BEUREAUCRATIC SUBROUTINES  #
#                             #
###############################



#################################################################
#
#  FUNCTION:  MIN
#
sub MIN
{
    my $a = shift;
    my $b = shift;
    return $a if ($a < $b);
    return $b;
}



#################################################################
#
#  FUNCTION:  MAX
#
sub MAX
{
    my $a = shift;
    my $b = shift;
    return $a if ($a > $b);
    return $b;
}



#################################################################
#
#  FUNCTION:  RunSystemCommand
#
sub RunSystemCommand
{
    my $command = shift;
    if (system($command)) { die "\n  ERROR:  System command '$command' failed during execution\n\n"; }
}



#################################################################
#
#  FUNCTION:  OpenSystemCommand
#
sub OpenSystemCommand
{
    my $command = shift;
    if ($command !~ /\s+\|\s*$/) { $command = $command.' |'; }
    open(my $command_output,$command) || die "\n  ERROR:  Failed to open output from system command '$command'\n\n";
    return $command_output;
}



#################################################################
#
#  FUNCTION:  OpenInputFile
#
sub OpenInputFile
{
    my $filename = shift;
    if (!(-e $filename)) { die "\n  ERROR:  Failed to locate input file '$filename'\n\n"; }
    open(my $filehandle,'<',$filename) || die "\n  ERROR:  Failed to open input file '$filename'\n\n";
    return $filehandle;
}



#################################################################
#
#  FUNCTION:  NameOutputFile
#
sub NameOutputFile
{
    my $intended_name = shift;
    my $basename;
    my $extension;
    if ($intended_name =~ /(\S+)(\.[^\.]+)$/) {
	$basename = $1;
	$extension = $2;
    } else {
	$basename = $intended_name;
	$extension = '';
    }
    my $filename = $basename.$extension;
    my $attempt = 1;
    while (-e $filename) {
	$attempt++;
	$filename = $basename.'_'.$attempt.$extension;
    }
    return $filename;
}



#################################################################
#
#  FUNCTION:  OpenOutputFile
#
sub OpenOutputFile
{
    my $filename = shift;
    $filename = NameOutputFile($filename);
    open(my $filehandle,'>',$filename) || die "\n  ERROR:  Failed to open output file '$filename'\n\n";
    return $filehandle;
}



#################################################################
#
#  FUNCTION:  OpenFileToAppendTo
#
sub OpenFileToAppendTo
{
    my $filename = shift;
    open(my $filehandle,'>>',$filename) || die "\n  ERROR:  Failed to open output file '$filename' (for appending)\n\n";
    return $filehandle;
}



#################################################################
#
#  FUNCTION:  ConfirmDirectory
#
sub ConfirmDirectory
{
    my $dirname = shift;
    if (!(-d $dirname)) { die "\n  ERROR:  Failed to locate directory '$dirname'\n\n"; }
    if ($dirname !~ /\/$/) { $dirname = $dirname.'/'; }
    return $dirname;
}



#################################################################
#
#  FUNCTION:  DeleteFile
#
sub DeleteFile
{
    my $filename;
    if (-e $filename) {
	my $rm_command = 'rm '.$filename;
	RunSystemCommand($rm_command);
    }
}



#################################################################
#
#  FUNCTION:  OpenDirectory
#
sub OpenDirectory
{
    my $dirname = shift;
    $dirname = ConfirmDirectory($dirname);
    opendir(my $dirhandle,$dirname) || die "\n  ERROR:  Failed to open directory '$dirname'\n\n";
    return $dirhandle;
}



#################################################################
#
#  FUNCTION:  NameDirectory
#
sub NameDirectory
{
    my $intended_name = shift;
    $intended_name =~ s/\/$//;
    my $dirname = $intended_name;
    my $attempt = 1;
    while (-d $dirname) {
	$attempt++;
	$dirname = $intended_name.'_'.$attempt;
    }
    $dirname = $dirname.'/';
    return $dirname;
}



#################################################################
#
#  FUNCTION:  CreateDirectory
#
sub CreateDirectory
{
    my $dirname = shift;
    $dirname = NameDirectory($dirname);
    RunSystemCommand("mkdir $dirname");
    if (!(-d $dirname)) { die "\n  ERROR:  Creation of directory '$dirname' failed\n\n"; }
    return $dirname;
}


##   END OF FILE   ##














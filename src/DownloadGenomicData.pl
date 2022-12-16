#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


#  FUNCTIONAL SUBROUTINES
#
sub ParseSpeciesTreeFile;
sub RecursiveTreeCheck;
sub RecursiveTreeReduce;
sub GenNameSwapHashes;
sub GetSpeciesNamePairs;
sub GetStaticLatinToEnglish;
sub GetStaticLatinTree;
sub GetStaticEnglishTree;


#  BEUREAUCRATIC SUBROUTINES
#
#  NOTE: Because I think this will be useful outside of Mirage, I'm not going to
#        link it to the src/Bureaucracy.pm file, even though that essentially
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
if (scalar(@ARGV) < 1){
    print "\n";
    print "  USAGES:\n";
    print "\n";
    print "  % ./DownloadGenomicData.pl [OPT.s] [file]   <-- Downloads genome data for species listed in file\n";
    print "  % ./DownloadGenomicData.pl [OPT.s] --full   <-- Downloads ALL UCSC genome data (this is a LOT)\n";
    print "\n";
    print "  OPT.s:\n";
    print "\n";
    print "    -outdirname [string] : Name the output directory (default:Genomic-Data)\n";
    print "    --parse-species-tree : Instead of downloading genomic data, parse phylip tree from NCBI (see README)\n";
    print "\n";
    print "  Help :\n";
    print "    --list-known-species : List species that this script has existing knowledge of\n";
    die   "\n";
}


# See if the user wants to try their hand at providing some commandline
# arguments.
my $species_list_fname = 0;
my $outdirname = 'Genomic-Data';
my $relative_outdir = 1;
for (my $opt_num = 0; $opt_num < scalar(@ARGV); $opt_num++) {

    if (lc($ARGV[$opt_num]) =~ /\-outdirname/) {

	$outdirname = $ARGV[++$opt_num];
	if ($outdirname =~ /^\// || $outdirname =~ /^\~/) {
	    $relative_outdir = 0;
	}

    } elsif ($opt_num == scalar(@ARGV)-1) {

	if (-e $ARGV[$opt_num]) {
	    $species_list_fname = $ARGV[$opt_num];
	} elsif (lc($ARGV[$opt_num]) !~ /\-full$/) {
	    die "\n  ERROR:  Either a species list file or '--full' must be provided as the final argument\n\n";
	}
	
    } elsif (lc($ARGV[$opt_num]) =~ /\-parse\-species\-tree/) {

	if (@ARGV < 3) {
	    die "\n  USAGE:  ./DownloadGenomicData.pl --parse-species-tree [tree.phy] [scientific-to-English-names.txt OR -]\n\n";
	}
	
	my ($latin_tree_str,$english_tree_str)
	    = ParseSpeciesTreeFile($ARGV[scalar(@ARGV)-2],$ARGV[scalar(@ARGV)-1]);
	
	print "\n\n";
	print "Latin Species Tree\n";
	print "------------------\n";
	print "$latin_tree_str\n";
	print "\n\n";
	print "English Species Tree (if known)\n";
	print "-------------------\n";
	print "$english_tree_str\n";
	print "\n\n";
    
	exit(0);

    } elsif (lc($ARGV[$opt_num]) =~ /\-list\-known\-species/) {

	my @EnglishNames;
	my @LatinNames;
	my $longest_english_name_len = 0;
	foreach my $species_name_pair (@{GetSpeciesNamePairs()}) {

	    $species_name_pair =~ s/\_/ /g;
	    $species_name_pair =~ /^([^\|]+)\|([^\|]+)$/;
	    my $latin_name = $1;
	    my $english_name = $2;
	    
	    push(@LatinNames,$latin_name);
	    push(@EnglishNames,$english_name);

	    if (length($english_name) > $longest_english_name_len) {
		$longest_english_name_len = length($english_name);
	    }

	}

	print "\n";
	for (my $i=0; $i<scalar(@EnglishNames); $i++) {
	    my $english_name = $EnglishNames[$i];
	    while (length($english_name) < $longest_english_name_len) {
		$english_name = $english_name.' ';
	    }
	    my $latin_name = $LatinNames[$i];
	    print "  $english_name ($latin_name)\n";
	}
	print "\n";
	print "  NOTE:  When making a list of species to download, include either\n";
	print "         the Latin OR the English name (not both).\n";
	print "\n";

	exit(0);

    } else {
	print "  Unrecognized option '$ARGV[$opt_num]' ignored\n";
    }
    
}


# We'll want to have these on-hand ASAP!
my ($latin_to_english_ref,$english_to_latin_ref) = GenNameSwapHashes('-');
my %LatinToEnglish = %{$latin_to_english_ref};
my %EnglishToLatin = %{$english_to_latin_ref};


# See if there's a file to parse
my %SpeciesToDownload;
if ($species_list_fname) {

    my $inf = OpenInputFile($species_list_fname);
    while (my $line = <$inf>) {

	# Eliminate any leading or trailing whitespace
	$line =~ s/\n|\r//g;
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;

	# Convert any spaces to underscores
	$line =~ s/\s/\_/g;

	my $species_name = lc($line);
	if ($EnglishToLatin{$species_name}) {
	    $species_name = $EnglishToLatin{$species_name};
	}
	
	$SpeciesToDownload{$species_name} = 1;

    }
    
} else {

    foreach my $species_name (keys %LatinToEnglish) {
	$SpeciesToDownload{$species_name} = 1;
    }

}


# Make an output directory where we can save all our fun files!
$outdirname = CreateDirectory($outdirname);


# Now it's time to do the terrible work of downloading some GD html
my $ucsc_dir_fname = $outdirname.'ucsc-dir.html';
my $genome_index_url = "https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/";
RunSystemCommand("wget --quiet -O $ucsc_dir_fname $genome_index_url");

# We'll use a standard filename for temporary html downloads
my $temp_html_fname = $outdirname.'temp.html';

# Scan through the directory file and find the links to all genome / gtf
# subdirectories containing species data we're interested in
my $UCSCf = OpenInputFile($ucsc_dir_fname);
my %SpeciesToGenomes;
my %SpeciesToGTFs;
my %TroubleSpecies;

# Let the user know what we're doing, since there'll be a bit of time
# spent performing quiet downloads
print "\n  Scanning UCSC Genome Browser for available genomes...\n";

# Get to the list of current genomes
while (my $line = <$UCSCf>) {
    last if ($line =~ /\<pre\>/);
}

if (eof($UCSCf)) {
    die "\n  ERROR: Failed to make sense of the HTML index file '$ucsc_dir_fname'\n\n";
}

while (my $line = <$UCSCf>) {

    $line =~ s/\n|\r//g;
    last if ($line =~ /\<\/pre\>/);
    
    next if ($line !~ /\<a href\=\"([^\"]+)\"\>([^\/|\<]+)\/?\</);
    my $species_dir_url = $1;
    my $latin_species_name = lc($2);

    next if (!$SpeciesToDownload{$latin_species_name});
    
    # What would an English-speaking person normally call this species?
    my $english_species_name = $LatinToEnglish{$latin_species_name};


    # 1. Identify the genome download link
    
    my $genome_dir_fname = $genome_index_url.$species_dir_url.'bigZips/';

    # Download the directory
    if (system("wget --quiet -O $temp_html_fname $genome_dir_fname")) {
	$TroubleSpecies{$latin_species_name} = 1;
	next;
    }

    # We're going to be a little funny -- look for 'xyz.2bit' and then swap
    # the '.2bit' with '.fa.gz'
    my $BigZips = OpenInputFile($temp_html_fname);
    my $genome_fname = 0;
    my $reduced_name = 0;
    while (my $bz_line = <$BigZips>) {
	if ($bz_line =~ /(\S+)\.2bit/) {
	    $reduced_name = $1;
	    $genome_fname = $genome_dir_fname.$reduced_name.'.fa.gz';
	    last;
	}
    }
    close($BigZips);

    # If we hit the end of the file and didn't find a genome download link,
    # note this as a "trouble species" and move on
    if (!$genome_fname) {
	$TroubleSpecies{$latin_species_name} = 1;
	next;
    }

    
    # GENOME AHOY!
    $SpeciesToGenomes{$latin_species_name} = $genome_fname;

    
    # 2. Pull in any gtfs that might be floating around

    # If there's a directory with gtfs it'll be named 'genes' (assuming conventions
    # from time of writing).
    my $gtf_dir_fname = $genome_dir_fname.'genes/';
    next if (system("wget --quiet -O $temp_html_fname $gtf_dir_fname"));
    
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

    # GTFs ahoy?
    if ($gtf_list_str) {
	# GTFs AHOY!
	$SpeciesToGTFs{$latin_species_name} = $gtf_list_str;
    }

}
close($UCSCf);


# No need for the main UCSC directory file or the temp html file anymore!
RunSystemCommand("rm $ucsc_dir_fname") if (-e $ucsc_dir_fname);
RunSystemCommand("rm $temp_html_fname") if (-e $temp_html_fname);


# If we didn't get any genomes that isn't ideal...
if (scalar(keys %SpeciesToGenomes) == 0) {
    RunSystemCommand("rmdir $outdirname");
    print "\n";
    print "  No genomes downloaded (likely reason: failed to match listed species names)\n";
    print "  Program exiting\n";
    die "\n";
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
my %LatinToSuper; # "English name (Latin name/shorthand)"
my $longest_english_name_len = 0;
my $longest_latin_name_len = 0;
my $longest_supername_len = 0;
my $longest_genome_fname_len = 0;

foreach my $latin_species_name (sort keys %SpeciesToGenomes) {

    if (length($latin_species_name) > $longest_latin_name_len) {
	$longest_latin_name_len = length($latin_species_name);
    }

    my $genome_link = $SpeciesToGenomes{$latin_species_name};

    # Just in case some are '.fasta.gz'
    $genome_link =~ /\/([^\/]+)(\.[^\.]+)\.gz$/;
    my $species_shorthand = $1;
    my $genome_extension = $2;

    my $english_species_name = $LatinToEnglish{$latin_species_name};
    if (length($english_species_name) > $longest_english_name_len) {
	$longest_english_name_len = length($english_species_name);
    }

    # What's your "supername"?
    my $supername = $english_species_name.' ('.$latin_species_name.'/'.$species_shorthand.')';
    if (length($supername) > $longest_supername_len) {
	$longest_supername_len = length($supername)
    }
    $LatinToSuper{$latin_species_name} = $supername;

    my $dl_genome_fname = $dl_genome_dirname.$species_shorthand.$genome_extension.'.gz';

    # We'll use rsync, since that's what UCSC recommends
    $genome_link =~ s/https?\:/rsync\:/;
    RunSystemCommand("rsync -a -P $genome_link $dl_genome_fname");

    # Unpack that sucker!
    RunSystemCommand("gunzip $dl_genome_fname");

    # Record the full path to the file
    $dl_genome_fname =~ s/\.gz$//g;
    $dl_genome_fname = $pwd.$dl_genome_fname if ($relative_outdir);
    $SpeciesToGenomes{$latin_species_name} = $dl_genome_fname;

    if (length($dl_genome_fname) > $longest_genome_fname_len) {
	$longest_genome_fname_len = length($dl_genome_fname);
    }

    # If there isn't at least one gtf for this species, then we're moving on!
    next if (!$SpeciesToGTFs{$latin_species_name});

    my @GTFList = split(/\|/,$SpeciesToGTFs{$latin_species_name});
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
    if ($relative_outdir) {
	$SpeciesToGTFs{$latin_species_name} =
	    $pwd.$dl_gtf_dirname.$species_shorthand.'.gtf';
    } else {
	$SpeciesToGTFs{$latin_species_name} =
	    $dl_gtf_dirname.$species_shorthand.'.gtf';
    }
    
    # We got our big honkin' nasty boi, so let's clean up our lil' beepin' chill gals!
    foreach my $dl_gtf_fname (@GTFList) {
	RunSystemCommand("rm $dl_gtf_fname");
    }
    
}


# WOWEE! How's about that -- we've done all of our downloading!
# Now it's time to build a SpeciesGuide file and let the user know
# how we did overall.


# ... But first, how about we decide what our species trees should look like?
#
my ($english_tree,$latin_tree)
    = RecursiveTreeReduce(GetStaticLatinTree());

# Species guides comin' atchya!
my $LatinGuide = OpenOutputFile($outdirname.'Latin-Species-Guide');
my $EnglishGuide = OpenOutputFile($outdirname.'English-Species-Guide');

if (scalar(split(/\,/,$latin_tree)) > 2) {
    print $LatinGuide "$latin_tree\n";
    print $EnglishGuide "$english_tree\n";
}

foreach my $latin_species_name (sort keys %SpeciesToGenomes) {

    my $genome = $SpeciesToGenomes{$latin_species_name};

    my $gtf = '-';
    if ($SpeciesToGTFs{$latin_species_name}) {
	$gtf = $SpeciesToGTFs{$latin_species_name};
    }

    my $english_species_name = $LatinToEnglish{$latin_species_name};
    while (length($english_species_name) < $longest_english_name_len) {
	$english_species_name = $english_species_name.' ';
    }

    while (length($latin_species_name) < $longest_latin_name_len) {
	$latin_species_name = $latin_species_name.' ';
    }

    while (length($genome) < $longest_genome_fname_len) {
	$genome = $genome.' ';
    }

    print $LatinGuide "$latin_species_name $genome $gtf\n";
    print $EnglishGuide "$english_species_name $genome $gtf\n";

}
close($LatinGuide);
close($EnglishGuide);


# If any of the requested species weren't in UCSC, we'll want to be sure to
# add them to the 'TroubleSpecies' hash
foreach my $species (keys %SpeciesToDownload) {
    if (!$SpeciesToGenomes{$species}) {
	$TroubleSpecies{$species} = 1;
    }
}

my %AllSpecies;
foreach my $species (keys %SpeciesToGenomes) { $AllSpecies{$species} = 1; }
foreach my $species (keys %TroubleSpecies  ) { $AllSpecies{$species} = 1; }

# Oh, let's give even the "trouble species" supernames -- it's a fundamental right!
foreach my $species (keys %TroubleSpecies) {
    $LatinToSuper{$species} = $species;
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

    my $supername = $LatinToSuper{$species};

    # Let's get pretty with it!
    $supername =~ s/\_/ /g;
    my @NameComponents = split(/ /,$supername);
    $supername = '';
    for (my $i=0; $i<scalar(@NameComponents); $i++) {
	if ($NameComponents[$i] !~ /\// && $NameComponents[$i] =~ /^([a-z])/) {
	    my $first_letter = uc($1);
	    $NameComponents[$i] =~ s/^[a-z]/$first_letter/;
	} elsif ($NameComponents[$i] =~ /^\(([a-z])/) {
	    my $first_letter = uc($1);
	    $NameComponents[$i] =~ s/^\([a-z]/\($first_letter/;
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
#  FUNCTION:  ParseSpeciesTreeFile
#
sub ParseSpeciesTreeFile
{
    
    my $treefname = shift;
    my $mapfname  = shift;

    my $treef = OpenInputFile($treefname);
    my $giant_str = '';
    while (my $line = <$treef>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	if ($line =~ /\)/) {
	    $line = '),'; # We only want species, not families etc.
	}
	$line =~ s/\'//g;
	$line =~ s/\:\d+//g;
	$line =~ s/\s/\_/g;

	$giant_str = $giant_str.$line;

    }
    close($treef);

    $giant_str =~ s/\,$//;
    return RecursiveTreeCheck(lc($giant_str));

}






#################################################################
#
#  FUNCTION:  RecursiveTreeCheck
#
sub RecursiveTreeCheck
{
    my $in_tree_str = shift;
    my @TreeChars = split(//,$in_tree_str);

    # What we'll do here is sort all characters into either being
    # at this sub-tree's node, or belonging to a subtree of their own
    my $tree_depth = 0;
    my $node_str = '';
    my $subtree_str;
    my @SubTreeStrs;
    for (my $char_pos=1; $char_pos < scalar(@TreeChars)-1; $char_pos++) {

	my $char = $TreeChars[$char_pos];

	if ($char eq '(') {
	    $subtree_str = '' if (!$tree_depth);
	    $tree_depth++;
	}

	if ($tree_depth) { $subtree_str = $subtree_str.$char; }
	else             { $node_str    =    $node_str.$char; }

	if ($char eq ')') {
	    $tree_depth--;
	    push(@SubTreeStrs,$subtree_str) if (!$tree_depth);
	}
	
    }

    while ($node_str =~ /\,\,/) {
	$node_str =~ s/\,\,/\,/;
    }

    my $latin_out_tree_str = '';
    my $english_out_tree_str = '';
    foreach my $latin_species (split(/\,/,$node_str)) {

	my $english_species = $LatinToEnglish{$latin_species};
	
	if ($latin_out_tree_str) {
	    $latin_out_tree_str = '('.$latin_out_tree_str.','.$latin_species.')';
	    $english_out_tree_str = '('.$english_out_tree_str.','.$english_species.')';
	} else {
	    $latin_out_tree_str = $latin_species;
	    $english_out_tree_str = $english_species;
	}

    }

    my $latin_out_subtree_str = '';
    my $english_out_subtree_str = '';
    for (my $i=0; $i<scalar(@SubTreeStrs); $i++) {

	my ($latin_subtree_str,$english_subtree_str)
	    = RecursiveTreeCheck($SubTreeStrs[$i]);

	if ($latin_out_subtree_str) {
	    $latin_out_subtree_str = '('.$latin_out_subtree_str.','.$latin_subtree_str.')';
	    $english_out_subtree_str = '('.$english_out_subtree_str.','.$english_subtree_str.')';
	} else {
	    $latin_out_subtree_str = $latin_subtree_str;
	    $english_out_subtree_str = $english_subtree_str;
	}

    }

    if ($latin_out_subtree_str) {
	$latin_out_tree_str = '('.$latin_out_tree_str.','.$latin_out_subtree_str.')';
	$english_out_tree_str = '('.$english_out_tree_str.','.$english_out_subtree_str.')';
    }

    return ($latin_out_tree_str,$english_out_tree_str);
    
}





#################################################################
#
#  FUNCTION:  RecursiveTreeReduce
#
sub RecursiveTreeReduce
{
    my $in_tree_str = shift;
    my @TreeChars = split(//,$in_tree_str);
    
    # Sadly, I'm essentially repeating what's in 'RecursiveTreeCheck' but with
    # the additional bit of computation involved in reacting to absent species
    my $tree_depth = 0;
    my $node_str = '';
    my $subtree_str;
    my @SubTreeStrs;
    for (my $char_pos=1; $char_pos < scalar(@TreeChars)-1; $char_pos++) {

	my $char = $TreeChars[$char_pos];

	if ($char eq '(') {
	    $subtree_str = '' if (!$tree_depth);
	    $tree_depth++;
	}

	if ($tree_depth) { $subtree_str = $subtree_str.$char; }
	else             { $node_str    =    $node_str.$char; }

	if ($char eq ')') {
	    $tree_depth--;
	    push(@SubTreeStrs,$subtree_str) if (!$tree_depth);
	}
	
    }

    # Because we're going off our pre-formatted string, we have a guarantee that
    # no node will have more than two species
    my $english_node_str = '';
    my $latin_node_str = '';
    foreach my $latin_species (split(/\,/,$node_str)) {

	next if (!$SpeciesToGenomes{$latin_species});
	
	my $english_species = $LatinToEnglish{$latin_species};
	if ($latin_node_str) {
	    $latin_node_str = '('.$latin_node_str.','.$latin_species.')';
	    $english_node_str = '('.$english_node_str.','.$english_species.')';
	} else {
	    $latin_node_str = $latin_species;
	    $english_node_str = $english_species;
	}
	
    }

    # Get the subtree string(s) (if there are any)
    my $english_subtree_str = '';
    my $latin_subtree_str = '';
    foreach $subtree_str (@SubTreeStrs) {

	my ($english_sub_str,$latin_sub_str)
	    = RecursiveTreeReduce($subtree_str);
	next if (!$english_sub_str);
	
	if ($english_subtree_str) {
	    $english_subtree_str = '('.$english_subtree_str.','.$english_sub_str.')';
	    $latin_subtree_str = '('.$latin_subtree_str.','.$latin_sub_str.')';
	} else {
	    $english_subtree_str = $english_sub_str;
	    $latin_subtree_str = $latin_sub_str;
	}
	
    }

    # Put it all together!
    my $english_out_str = '';
    my $latin_out_str = '';
    if ($english_node_str) {
	if ($english_subtree_str) {
	    $english_out_str = '('.$english_node_str.','.$english_subtree_str.')';
	    $latin_out_str = '('.$latin_node_str.','.$latin_subtree_str.')';
	} else {
	    $english_out_str = $english_node_str;
	    $latin_out_str = $latin_node_str;
	}
    } elsif ($english_subtree_str) {
	$english_out_str = $english_subtree_str;
	$latin_out_str = $latin_subtree_str;
    }

    return ($english_out_str,$latin_out_str);

}






#################################################################
#
#  FUNCTION:  GenNameSwapHashes
#
sub GenNameSwapHashes
{
    my $fname = shift;

    my %LatinToEnglish;
    my %EnglishToLatin;
    if ($fname eq '-') {
    
	my $species_name_pairs_ref = GetSpeciesNamePairs();
	foreach my $species_name_pair (@{$species_name_pairs_ref}) {
	    
	    $species_name_pair =~ /^([^\|]+)\|([^\|]+)$/;
	    my $latin_name = $1;
	    my $english_name = $2;

	    $LatinToEnglish{$latin_name} = $english_name;
	    $EnglishToLatin{$english_name} = $latin_name;
	
	}

    } else {

	my $mapf = OpenInputFile($fname);
	my %LatinToEnglish;
	while (my $line = <$mapf>) {
	    
	    if ($line =~ /^\s*([^\=]+)\s*\=\s*([^\=]+)\s*$/) {
		
		my $latin_name   = lc($1);
		my $english_name = lc($2);
		
		$latin_name =~ s/\s/\_/g;
		$english_name =~ s/\s/\_/g;
		
		$LatinToEnglish{$latin_name} = $english_name;
		$EnglishToLatin{$english_name} = $latin_name;
		
	    }
	    
	}
	close($mapf);

    }

    return(\%LatinToEnglish,\%EnglishToLatin);

}







#################################################################
#
#  FUNCTION:  GetSpeciesNamePairs
#
sub GetSpeciesNamePairs
{
    my @SpeciesNamePairs;
    push(@SpeciesNamePairs,'anopheles_gambiae|a._gambiae');
    push(@SpeciesNamePairs,'vicugna_pacos|alpaca');
    push(@SpeciesNamePairs,'alligator_mississippiensis|american_alligator');
    push(@SpeciesNamePairs,'dasypus_novemcinctus|armadillo');
    push(@SpeciesNamePairs,'gadus_morhua|atlantic_cod');
    push(@SpeciesNamePairs,'papio_anubis|baboon');
    push(@SpeciesNamePairs,'bison_bison_bison|bison');
    push(@SpeciesNamePairs,'pan_paniscus|bonobo');
    push(@SpeciesNamePairs,'apteryx_australis|brown_kiwi');
    push(@SpeciesNamePairs,'melopsittacus_undulatus|budgerigar');
    push(@SpeciesNamePairs,'otolemur_garnettii|bushbaby');
    push(@SpeciesNamePairs,'felis_catus|cat');
    push(@SpeciesNamePairs,'gallus_gallus|chicken');
    push(@SpeciesNamePairs,'pan_troglodytes|chimp');
    push(@SpeciesNamePairs,'cricetulus_griseus|chinese_hamster');
    push(@SpeciesNamePairs,'manis_pentadactyla|chinese_pangolin');
    push(@SpeciesNamePairs,'ciona_intestinalis|ciona');
    push(@SpeciesNamePairs,'latimeria_chalumnae|coelacanth');
    push(@SpeciesNamePairs,'bos_taurus|cow');
    push(@SpeciesNamePairs,'macaca_fascicularis|crab-eating_macaque');
    push(@SpeciesNamePairs,'drosophila_melanogaster|d._melanogaster');
    push(@SpeciesNamePairs,'canis_lupus_familiaris|dog');
    push(@SpeciesNamePairs,'tursiops_truncatus|dolphin');
    push(@SpeciesNamePairs,'loxodonta_africana|elephant');
    push(@SpeciesNamePairs,'callorhinchus_milii|elephant_shark');
    push(@SpeciesNamePairs,'mustela_putorius_furo|ferret');
    push(@SpeciesNamePairs,'takifugu_rubripes|fugu');
    push(@SpeciesNamePairs,'thamnophis_sirtalis|garter_snake');
    push(@SpeciesNamePairs,'nomascus_leucogenys|gibbon');
    push(@SpeciesNamePairs,'aquila_chrysaetos_canadensis|golden_eagle');
    push(@SpeciesNamePairs,'rhinopithecus_roxellana|golden_snub-nosed_monkey');
    push(@SpeciesNamePairs,'gorilla_gorilla_gorilla|gorilla');
    push(@SpeciesNamePairs,'chlorocebus_sabaeus|green_monkey');
    push(@SpeciesNamePairs,'cavia_porcellus|guinea_pig');
    push(@SpeciesNamePairs,'neomonachus_schauinslandi|hawaiian_monk_seal');
    push(@SpeciesNamePairs,'erinaceus_europaeus|hedgehog');
    push(@SpeciesNamePairs,'equus_caballus|horse');
    push(@SpeciesNamePairs,'homo_sapiens|human');
    push(@SpeciesNamePairs,'dipodomys_ordii|kangaroo_rat');
    push(@SpeciesNamePairs,'petromyzon_marinus|lamprey');
    push(@SpeciesNamePairs,'anolis_carolinensis|lizard');
    push(@SpeciesNamePairs,'galeopterus_variegatus|malayan_flying_lemur');
    push(@SpeciesNamePairs,'trichechus_manatus_latirostris|manatee');
    push(@SpeciesNamePairs,'callithrix_jacchus|marmoset');
    push(@SpeciesNamePairs,'oryzias_latipes|medaka');
    push(@SpeciesNamePairs,'geospiza_fortis|medium_ground_finch');
    push(@SpeciesNamePairs,'pteropus_vampyrus|megabat');
    push(@SpeciesNamePairs,'myotis_lucifugus|microbat');
    push(@SpeciesNamePairs,'balaenoptera_acutorostrata_scammoni|minke_whale');
    push(@SpeciesNamePairs,'mus_musculus|mouse');
    push(@SpeciesNamePairs,'microcebus_murinus|mouse_lemur');
    push(@SpeciesNamePairs,'heterocephalus_glaber|naked_mole-rat');
    push(@SpeciesNamePairs,'oreochromis_niloticus|nile_tilapia');
    push(@SpeciesNamePairs,'pongo_pygmaeus_abelii|orangutan');
    push(@SpeciesNamePairs,'chrysemys_picta_bellii|painted_turtle');
    push(@SpeciesNamePairs,'ailuropoda_melanoleuca|panda');
    push(@SpeciesNamePairs,'sus_scrofa|pig');
    push(@SpeciesNamePairs,'ochotona_princeps|pika');
    push(@SpeciesNamePairs,'ornithorhynchus_anatinus|platypus');
    push(@SpeciesNamePairs,'nasalis_larvatus|proboscis_monkey');
    push(@SpeciesNamePairs,'oryctolagus_cuniculus|rabbit');
    push(@SpeciesNamePairs,'rattus_norvegicus|rat');
    push(@SpeciesNamePairs,'macaca_mulatta|rhesus');
    push(@SpeciesNamePairs,'procavia_capensis|rock_hyrax');
    push(@SpeciesNamePairs,'aplysia_californica|sea_hare');
    push(@SpeciesNamePairs,'ovis_aries|sheep');
    push(@SpeciesNamePairs,'sorex_araneus|shrew');
    push(@SpeciesNamePairs,'choloepus_hoffmanni|sloth');
    push(@SpeciesNamePairs,'enhydra_lutris_nereis|southern_sea_otter');
    push(@SpeciesNamePairs,'spermophilus_tridecemlineatus|squirrel');
    push(@SpeciesNamePairs,'saimiri_boliviensis|squirrel_monkey');
    push(@SpeciesNamePairs,'tarsius_syrichta|tarsier');
    push(@SpeciesNamePairs,'sarcophilus_harrisii|tasmanian_devil');
    push(@SpeciesNamePairs,'echinops_telfairi|tenrec');
    push(@SpeciesNamePairs,'nanorana_parkeri|tibetan_frog');
    push(@SpeciesNamePairs,'tupaia_belangeri|tree_shrew');
    push(@SpeciesNamePairs,'meleagris_gallopavo|turkey');
    push(@SpeciesNamePairs,'strongylocentrotus_purpuratus|urchin');
    push(@SpeciesNamePairs,'notamacropus_eugenii|wallaby');
    push(@SpeciesNamePairs,'ceratotherium_simum|white_rhinoceros');
    push(@SpeciesNamePairs,'xenopus_tropicalis|x._tropicalis');
    push(@SpeciesNamePairs,'taeniopygia_guttata|zebra_finch');
    push(@SpeciesNamePairs,'danio_rerio|zebrafish');
    
    return \@SpeciesNamePairs;

}







#################################################################
#
#  FUNCTION:  GetStaticLatinTree / GetStaticEnglishTree
#
sub GetStaticLatinTree
{
    return '((strongylocentrotus_purpuratus,aplysia_californica),(((((((alligator_mississippiensis,chrysemys_picta_bellii),latimeria_chalumnae),callorhinchus_milii),petromyzon_marinus),ciona_intestinalis),(((((thamnophis_sirtalis,anolis_carolinensis),((((((((((((galeopterus_variegatus,manis_pentadactyla),trichechus_manatus_latirostris),tupaia_belangeri),procavia_capensis),loxodonta_africana),echinops_telfairi),dasypus_novemcinctus),choloepus_hoffmanni),notamacropus_eugenii),sarcophilus_harrisii),ornithorhynchus_anatinus),((((((((pteropus_vampyrus,myotis_lucifugus),(((((ictidomys_tridecemlineatus,heterocephalus_glaber),cavia_porcellus),cricetulus_griseus),dipodomys_ordii),(rattus_norvegicus,mus_musculus))),(oryctolagus_cuniculus,ochotona_princeps)),(ceratotherium_simum,equus_caballus)),((((balaenoptera_acutorostrata_scammoni,vicugna_pacos),sus_scrofa),tursiops_truncatus),((bison_bison_bison,ovis_aries),bos_taurus))),((((neomonachus_schauinslandi,felis_catus),ailuropoda_melanoleuca),canis_lupus_familiaris),(enhydra_lutris_nereis,mustela_putorius_furo))),((((carlito_syrichta,nomascus_leucogenys),otolemur_garnettii),microcebus_murinus),(((((homo_sapiens,pongo_abelii),gorilla_gorilla_gorilla),(pan_troglodytes,pan_paniscus)),((((rhinopithecus_roxellana,chlorocebus_sabaeus),nasalis_larvatus),papio_anubis),(macaca_mulatta,macaca_fascicularis))),(saimiri_boliviensis,callithrix_jacchus)))),(sorex_araneus,erinaceus_europaeus)))),(((aquila_chrysaetos_canadensis,melopsittacus_undulatus),apteryx_australis),((taeniopygia_guttata,geospiza_fortis),(meleagris_gallopavo,gallus_gallus)))),(nanorana_parkeri,xenopus_tropicalis)),((((takifugu_rubripes,oreochromis_niloticus),oryzias_latipes),gadus_morhua),danio_rerio))),(drosophila_melanogaster,anopheles_gambiae)))';
}
sub GetStaticEnglishTree
{
    return '((urchin,sea_hare),(((((((american_alligator,painted_turtle),coelacanth),elephant_shark),lamprey),ciona),(((((garter_snake,lizard),((((((((((((malayan_flying_lemur,chinese_pangolin),manatee),tree_shrew),rock_hyrax),elephant),tenrec),armadillo),sloth),wallaby),tasmanian_devil),platypus),((((((((megabat,microbat),(((((ictidomys_tridecemlineatus,naked_mole-rat),guinea_pig),chinese_hamster),kangaroo_rat),(rat,mouse))),(rabbit,pika)),(white_rhinoceros,horse)),((((minke_whale,alpaca),pig),dolphin),((bison,sheep),cow))),((((hawaiian_monk_seal,cat),panda),dog),(southern_sea_otter,ferret))),((((carlito_syrichta,gibbon),bushbaby),mouse_lemur),(((((human,pongo_abelii),gorilla),(chimp,bonobo)),((((golden_snub-nosed_monkey,green_monkey),proboscis_monkey),baboon),(rhesus,crab-eating_macaque))),(squirrel_monkey,marmoset)))),(shrew,hedgehog)))),(((golden_eagle,budgerigar),brown_kiwi),((zebra_finch,medium_ground_finch),(turkey,chicken)))),(tibetan_frog,x._tropicalis)),((((fugu,nile_tilapia),medaka),atlantic_cod),zebrafish))),(d._melanogaster,a._gambiae)))';
}





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














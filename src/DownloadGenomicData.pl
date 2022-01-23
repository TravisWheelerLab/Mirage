#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


#  FUNCTIONAL SUBROUTINES
#
sub ParseSpeciesTreeFile;
sub RecursiveTreeCheck;
sub RecursiveTreeReduce;
sub GetLatinToCommon;
sub GetStaticLatinToCommon;
sub GetStaticLatinTree;
sub GetStaticCommonTree;


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
    die   "\n";
}


# See if the user wants to try their hand at providing some commandline
# arguments.
my $outdirname = 'Genomic-Data';
my $relative_outdir = 1;
for (my $opt_num = 0; $opt_num < scalar(@ARGV)-1; $opt_num++) {
    if (lc($ARGV[$opt_num]) =~ /\-outdirname/) {
	$outdirname = $ARGV[++$opt_num];
	if ($outdirname =~ /^\// || $outdirname =~ /^\~/) {
	    $relative_outdir = 0;
	}
    } elsif (lc($ARGV[$opt_num]) =~ /\-parse\-species\-tree/) {

	if (@ARGV < 3) {
	    die "\n  USAGE:  ./DownloadGenomicData.pl --parse-species-tree [tree.phy] [scientific-to-common-names.txt OR -]\n\n";
	}
	my ($latin_tree_str,$common_tree_str)
	    = ParseSpeciesTreeFile($ARGV[scalar(@ARGV)-2],$ARGV[scalar(@ARGV)-1]);
	
	print "\n\n";
	print "Latin Species Tree\n";
	print "------------------\n";
	print "$latin_tree_str\n";
	print "\n\n";
	print "Common Species Tree (if known)\n";
	print "-------------------\n";
	print "$common_tree_str\n";
	print "\n\n";
    
	exit(0);

    } else {
	print "  Unrecognized option '$ARGV[$opt_num]' ignored\n";
    }
}


# See if there's a file to parse
my %SpeciesToDownload;
my $download_all = 1;
if (lc($ARGV[scalar(@ARGV)-1]) !~ /\-full/) {

    $download_all = 0;
    
    my $inf = OpenInputFile($ARGV[scalar(@ARGV)-1]);
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
$outdirname = CreateDirectory($outdirname);


# Now it's time to do the terrible work of downloading some GD html
my $ucsc_dir_fname = $outdirname.'ucsc-dir.html';
RunSystemCommand("wget --quiet -O $ucsc_dir_fname https://hgdownload.soe.ucsc.edu/downloads.html");

# We'll use a standard filename for temporary html downloads
my $temp_html_fname = $outdirname.'temp.html';

# Scan through the directory file and find the links to all genome / gtf
# subdirectories containing species data we're interested in
my $UCSCf = OpenInputFile($ucsc_dir_fname);
my %SpeciesToGenomes;
my %SpeciesToSciNames;
my %SpeciesToGTFs;
my %TroubleSpecies;

my $ucscf_line = <$UCSCf>; # prime our reader

# Let the user know what we're doing, since there'll be a bit of time
# spent performing quiet downloads
print "\n  Scanning UCSC Genome Browser for available genomes...\n";

while (!eof($UCSCf)) {

    # Look for the next section break (these had BETTER stay consistent!)
    $ucscf_line =~ s/\n|\r//g;
    if ($ucscf_line !~ /\<\!\-\- (.+) Download/) {
	$ucscf_line = <$UCSCf>;
	next;
    }

    # What would a human being call this species?
    my $plain_species_name = lc($1);
    $plain_species_name =~ s/\s/\_/g;

    # What's your shorthand name?
    while ($ucscf_line = <$UCSCf>) {
	if ($ucscf_line =~ /\<\/h3\>/ || $ucscf_line =~ /\<\!\-\- .+ Download/) {
	    last;
	}
    }
    last if (eof($UCSCf));

    # There are some utility sections that use the same general formatting,
    # so we need to make sure that if we aren't looking at a genome then
    # we can jump ship at an appropriate time.
    my $shorthand_species_name;
    if ($ucscf_line =~ /\>([^\<]+)\<\/a\>\)\<\/h3\>/) {
	$shorthand_species_name = $1;
    } else {
	# NOTE: We don't pull in a new line in case we 'last'-ed to a '<!--'
	next;
    }

    # We'll want to grab the species' scientific name, for smartness
    my $desc_fname = 'https://hgdownload.soe.ucsc.edu/gbdb/'.$shorthand_species_name.'/html/description.html';
    my $scientific_name = $shorthand_species_name; # Placeholder
    if (!system("wget --quiet -O $temp_html_fname $desc_fname")) {

	# We'll look for the scientific name in either the genome ID (preferred)
	# or under the photo description.
	my $descf = OpenInputFile($temp_html_fname);
	my $photo_name = 0;
	my $gen_id_name = 0;
	while (my $desc_line = <$descf>) {
	    if (lc($desc_line) =~ /\<font size\=\-1\>\s*\<em\>([^\<]+)\<\/em\>\s*\<br\>/) {
		$photo_name = $1;
		$photo_name =~ s/^\s+//g;
		$photo_name =~ s/\s+$//g;
		$photo_name =~ s/\s/\_/g;
	    } elsif ($desc_line =~ /NCBI Genome ID\:/) {
		while (!eof($descf) && lc($desc_line) !~ /\<\/a\>\s*\(/) {
		    $desc_line = <$descf>;
		}
		last if (eof($descf));
		if (lc($desc_line) =~ /\<\/a\>\s*\(([^\)|\(]+)/) {
		    $gen_id_name = $1;
		    $gen_id_name =~ s/^\s+//;
		    $gen_id_name =~ s/\s+$//;
		    $gen_id_name =~ s/\s/\_/g;
		}
	    }
	}
	close($descf);

	# In case there are both gen_id and scientific names, we'll want
	# to make sure that we grab whichever one looks the most Latin
	# (i.e., has spaces).
	# NOTE that even though there seems to be some redundancies in how
	# this is coded it's reflecting a series of preferences, so we need
	# all this conditioning.
	if ($gen_id_name && $photo_name) {
	    if ($gen_id_name =~ /\_/) {
		$scientific_name = $gen_id_name;
	    } elsif ($photo_name =~ /\_/) {
		$scientific_name = $photo_name;
	    } else {
		$scientific_name = $gen_id_name;
	    }
	} elsif ($gen_id_name) {
	    $scientific_name = $gen_id_name;
	} elsif ($photo_name) {
	    $scientific_name = $photo_name;
	}
	
    }

    $SpeciesToSciNames{$plain_species_name} = $scientific_name;

    # If we have a list of species to look for, see if this dude's made the cut
    unless ($download_all || $SpeciesToDownload{$plain_species_name} || $SpeciesToDownload{$scientific_name}) {
	$ucscf_line = <$UCSCf>;
	next;
    }

    # Assuming everything holds, we don't really need to parse anything else,
    # since there's a standard convention for the URLs containing genome data
    my $genome_dir_fname = 'https://hgdownload.soe.ucsc.edu/goldenPath/'.$shorthand_species_name.'/bigZips/';


    # TIME FOR THE MONEY MOVES!
    
    
    # 1. Identify the genome download link

    # Download the directory
    if (system("wget --quiet -O $temp_html_fname $genome_dir_fname")) {
	my $trouble_key = $plain_species_name.' ('.$shorthand_species_name.')';
	$TroubleSpecies{$trouble_key} = 1;
	$ucscf_line = <$UCSCf>;
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
	$ucscf_line = <$UCSCf>;
	next;
    }

    # GENOME AHOY!
    $SpeciesToGenomes{$plain_species_name} = $genome_fname;

    
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

    # If there weren't any gtf files in our 'genes' dir, then I have two thoughts:
    # (1.) weird...?, and (2.) whatever
    if (!$gtf_list_str) {
	$ucscf_line = <$UCSCf>;
	next;
    }

    # GTFs AHOY!
    $SpeciesToGTFs{$plain_species_name} = $gtf_list_str;

    
    # Moving right along
    $ucscf_line = <$UCSCf>;

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
my %SpeciesToSuperNames; # english name + scientific name + shorthand
my $longest_name_len = 0;
my $longest_sci_name_len = 0;
my $longest_genome_len = 0;
my $longest_supername_len = 0;

foreach my $species (sort keys %SpeciesToGenomes) {

    $longest_name_len = length($species) if (length($species) > $longest_name_len);

    my $genome_link = $SpeciesToGenomes{$species};

    # Just in case some are '.fasta.gz'
    $genome_link =~ /\/([^\/]+)(\.[^\.]+)\.gz$/;
    my $species_shorthand = $1;
    my $genome_extension  = $2;

    my $scientific_name = $SpeciesToSciNames{$species};
    $longest_sci_name_len = length($scientific_name)
	if (length($scientific_name) > $longest_sci_name_len);

    # What's your "supername"?
    my $supername = $species.' ('.$scientific_name.'/'.$species_shorthand.')';
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
    $dl_genome_fname = $pwd.$dl_genome_fname if ($relative_outdir);
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
    if ($relative_outdir) {
	$SpeciesToGTFs{$species} = $pwd.$dl_gtf_dirname.$species_shorthand.'.gtf';
    } else {
	$SpeciesToGTFs{$species} = $dl_gtf_dirname.$species_shorthand.'.gtf';
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
my %ReducedSpeciesToSciNames;
foreach my $species (keys %SpeciesToGenomes) {
    $ReducedSpeciesToSciNames{$species} = $SpeciesToSciNames{$species};
}
my ($common_tree,$latin_tree)
    = RecursiveTreeReduce(GetStaticCommonTree(),\%ReducedSpeciesToSciNames);

# We'll count the number of species in our trees -- if there are less than
# three, then there really isn't any point in having a species tree
my $num_species_in_tree = scalar(split(/\,/,$common_tree));


# 1. Common name species guide
my $SpeciesGuide = OpenOutputFile($outdirname.'Common-Species-Guide');

print $SpeciesGuide "$common_tree\n" if ($num_species_in_tree > 2);
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

# 2. Scientific name species guide
$SpeciesGuide = OpenOutputFile($outdirname.'Latin-Species-Guide');

print $SpeciesGuide "$latin_tree\n" if ($num_species_in_tree > 2);
my %AllSciSpeciesNames;
foreach my $species (sort keys %SpeciesToGenomes) {

    my $scientific_name = $SpeciesToSciNames{$species};
    $AllSciSpeciesNames{$scientific_name} = 1;
    
    my $genome = $SpeciesToGenomes{$species};

    my $gtf = '-';
    if ($SpeciesToGTFs{$species}) {
	$gtf = $SpeciesToGTFs{$species};
    }

    while (length($scientific_name) < $longest_sci_name_len) {
	$scientific_name = $scientific_name.' ';
    }

    while (length($genome) < $longest_genome_len) {
	$genome = $genome.' ';
    }

    print $SpeciesGuide "$scientific_name $genome $gtf\n";

}
close($SpeciesGuide);

# If any of the requested species weren't in UCSC, we'll want to be sure to
# add them to the 'TroubleSpecies' hash
foreach my $species (keys %SpeciesToDownload) {
    if (!$SpeciesToGenomes{$species} && !$AllSciSpeciesNames{$species}) {
	$TroubleSpecies{$species} = 1;
    }
}

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

    my $latin_to_common_ref = GetLatinToCommon($mapfname);
    my %LatinToCommon = %{$latin_to_common_ref};

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
    return RecursiveTreeCheck(lc($giant_str),\%LatinToCommon);

}






#################################################################
#
#  FUNCTION:  RecursiveTreeCheck
#
sub RecursiveTreeCheck
{
    my $in_tree_str = shift;
    my @TreeChars = split(//,$in_tree_str);

    my $latin_to_common_ref = shift;
    my %LatinToCommon = %{$latin_to_common_ref};

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
    my $common_out_tree_str = '';
    foreach my $species (split(/\,/,$node_str)) {

	my $common_species;
	if ($LatinToCommon{$species}) {
	    $common_species = $LatinToCommon{$species};
	} else {
	    $common_species = $species;
	}
	
	if ($latin_out_tree_str) {
	    $latin_out_tree_str = '('.$latin_out_tree_str.','.$species.')';
	    $common_out_tree_str = '('.$common_out_tree_str.','.$common_species.')';
	} else {
	    $latin_out_tree_str = $species;
	    $common_out_tree_str = $common_species;
	}

    }

    my $latin_out_subtree_str = '';
    my $common_out_subtree_str = '';
    for (my $i=0; $i<scalar(@SubTreeStrs); $i++) {

	my ($latin_subtree_str,$common_subtree_str)
	    = RecursiveTreeCheck($SubTreeStrs[$i],\%LatinToCommon);

	if ($latin_out_subtree_str) {
	    $latin_out_subtree_str = '('.$latin_out_subtree_str.','.$latin_subtree_str.')';
	    $common_out_subtree_str = '('.$common_out_subtree_str.','.$common_subtree_str.')';
	} else {
	    $latin_out_subtree_str = $latin_subtree_str;
	    $common_out_subtree_str = $common_subtree_str;
	}

    }

    if ($latin_out_subtree_str) {
	$latin_out_tree_str = '('.$latin_out_tree_str.','.$latin_out_subtree_str.')';
	$common_out_tree_str = '('.$common_out_tree_str.','.$common_out_subtree_str.')';
    }

    return ($latin_out_tree_str,$common_out_tree_str);
    
}





#################################################################
#
#  FUNCTION:  RecursiveTreeReduce
#
sub RecursiveTreeReduce
{
    my $in_tree_str = shift;
    my @TreeChars = split(//,$in_tree_str);
    
    my $common_to_latin_ref = shift;
    my %CommonToLatin = %{$common_to_latin_ref};

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
    my $common_node_str = '';
    my $latin_node_str = '';
    foreach my $species (split(/\,/,$node_str)) {
	if ($CommonToLatin{$species}) {
	    if ($common_node_str) {
		$common_node_str = '('.$common_node_str.','.$species.')';
		$latin_node_str = '('.$latin_node_str.','.$CommonToLatin{$species}.')';
	    } else {
		$common_node_str = $species;
		$latin_node_str = $CommonToLatin{$species};
	    }
	}
    }

    # Get the subtree string(s) (if there are any)
    my $common_subtree_str = '';
    my $latin_subtree_str = '';
    foreach $subtree_str (@SubTreeStrs) {
	my ($common_sub_str,$latin_sub_str)
	    = RecursiveTreeReduce($subtree_str,\%CommonToLatin);
	next if (!$common_sub_str);
	if ($common_subtree_str) {
	    $common_subtree_str = '('.$common_subtree_str.','.$common_sub_str.')';
	    $latin_subtree_str = '('.$latin_subtree_str.','.$latin_sub_str.')';
	} else {
	    $common_subtree_str = $common_sub_str;
	    $latin_subtree_str = $latin_sub_str;
	}
    }

    # Put it all together!
    my $common_out_str = '';
    my $latin_out_str = '';
    if ($common_node_str) {
	if ($common_subtree_str) {
	    $common_out_str = '('.$common_node_str.','.$common_subtree_str.')';
	    $latin_out_str = '('.$latin_node_str.','.$latin_subtree_str.')';
	} else {
	    $common_out_str = $common_node_str;
	    $latin_out_str = $latin_node_str;
	}
    } elsif ($common_subtree_str) {
	$common_out_str = $common_subtree_str;
	$latin_out_str = $latin_subtree_str;
    }

    return ($common_out_str,$latin_out_str);

}






#################################################################
#
#  FUNCTION:  GetLatinToCommon
#
sub GetLatinToCommon
{
    my $fname = shift;
    return GetStaticLatinToCommon() if ($fname eq '-');

    my $mapf = OpenInputFile($fname);
    my %LatinToCommon;
    while (my $line = <$mapf>) {

	if ($line =~ /^([^\=]+)\=([^\=]+)$/) {

	    my $latin_name  = lc($1);
	    my $common_name = lc($2);

	    $latin_name =~ s/^\s+//;
	    $latin_name =~ s/\s+$//;
	    $latin_name =~ s/\s/\_/g;

	    $common_name =~ s/^\s+//;
	    $common_name =~ s/\s+$//;
	    $common_name =~ s/\s/\_/g;

	    $LatinToCommon{$latin_name} = $common_name;

	}

    }
    close($mapf);

    return \%LatinToCommon;
    
}






#################################################################
#
#  FUNCTION:  GetStaticLatinToCommon
#
sub GetStaticLatinToCommon
{
    my %LatinToCommon;

    $LatinToCommon{'anopheles_gambiae'} = 'a._gambiae';
    $LatinToCommon{'vicugna_pacos'} = 'alpaca';
    $LatinToCommon{'alligator_mississippiensis'} = 'american_alligator';
    $LatinToCommon{'dasypus_novemcinctus'} = 'armadillo';
    $LatinToCommon{'gadus_morhua'} = 'atlantic_cod';
    $LatinToCommon{'papio_anubis'} = 'baboon';
    $LatinToCommon{'bison_bison_bison'} = 'bison';
    $LatinToCommon{'pan_paniscus'} = 'bonobo';
    $LatinToCommon{'apteryx_australis'} = 'brown_kiwi';
    $LatinToCommon{'melopsittacus_undulatus'} = 'budgerigar';
    $LatinToCommon{'otolemur_garnettii'} = 'bushbaby';
    $LatinToCommon{'felis_catus'} = 'cat';
    $LatinToCommon{'gallus_gallus'} = 'chicken';
    $LatinToCommon{'pan_troglodytes'} = 'chimp';
    $LatinToCommon{'cricetulus_griseus'} = 'chinese_hamster';
    $LatinToCommon{'manis_pentadactyla'} = 'chinese_pangolin';
    $LatinToCommon{'ciona_intestinalis'} = 'ciona';
    $LatinToCommon{'latimeria_chalumnae'} = 'coelacanth';
    $LatinToCommon{'bos_taurus'} = 'cow';
    $LatinToCommon{'macaca_fascicularis'} = 'crab-eating_macaque';
    $LatinToCommon{'drosophila_melanogaster'} = 'd._melanogaster';
    $LatinToCommon{'canis_lupus_familiaris'} = 'dog';
    $LatinToCommon{'tursiops_truncatus'} = 'dolphin';
    $LatinToCommon{'loxodonta_africana'} = 'elephant';
    $LatinToCommon{'callorhinchus_milii'} = 'elephant_shark';
    $LatinToCommon{'mustela_putorius_furo'} = 'ferret';
    $LatinToCommon{'takifugu_rubripes'} = 'fugu';
    $LatinToCommon{'thamnophis_sirtalis'} = 'garter_snake';
    $LatinToCommon{'nomascus_leucogenys'} = 'gibbon';
    $LatinToCommon{'aquila_chrysaetos_canadensis'} = 'golden_eagle';
    $LatinToCommon{'rhinopithecus_roxellana'} = 'golden_snub-nosed_monkey';
    $LatinToCommon{'gorilla_gorilla_gorilla'} = 'gorilla';
    $LatinToCommon{'chlorocebus_sabaeus'} = 'green_monkey';
    $LatinToCommon{'cavia_porcellus'} = 'guinea_pig';
    $LatinToCommon{'neomonachus_schauinslandi'} = 'hawaiian_monk_seal';
    $LatinToCommon{'erinaceus_europaeus'} = 'hedgehog';
    $LatinToCommon{'equus_caballus'} = 'horse';
    $LatinToCommon{'homo_sapiens'} = 'human';
    $LatinToCommon{'dipodomys_ordii'} = 'kangaroo_rat';
    $LatinToCommon{'petromyzon_marinus'} = 'lamprey';
    $LatinToCommon{'anolis_carolinensis'} = 'lizard';
    $LatinToCommon{'galeopterus_variegatus'} = 'malayan_flying_lemur';
    $LatinToCommon{'trichechus_manatus_latirostris'} = 'manatee';
    $LatinToCommon{'callithrix_jacchus'} = 'marmoset';
    $LatinToCommon{'oryzias_latipes'} = 'medaka';
    $LatinToCommon{'geospiza_fortis'} = 'medium_ground_finch';
    $LatinToCommon{'pteropus_vampyrus'} = 'megabat';
    $LatinToCommon{'myotis_lucifugus'} = 'microbat';
    $LatinToCommon{'balaenoptera_acutorostrata_scammoni'} = 'minke_whale';
    $LatinToCommon{'mus_musculus'} = 'mouse';
    $LatinToCommon{'microcebus_murinus'} = 'mouse_lemur';
    $LatinToCommon{'heterocephalus_glaber'} = 'naked_mole-rat';
    $LatinToCommon{'oreochromis_niloticus'} = 'nile_tilapia';
    $LatinToCommon{'pongo_pygmaeus_abelii'} = 'orangutan';
    $LatinToCommon{'chrysemys_picta_bellii'} = 'painted_turtle';
    $LatinToCommon{'ailuropoda_melanoleuca'} = 'panda';
    $LatinToCommon{'sus_scrofa'} = 'pig';
    $LatinToCommon{'ochotona_princeps'} = 'pika';
    $LatinToCommon{'ornithorhynchus_anatinus'} = 'platypus';
    $LatinToCommon{'nasalis_larvatus'} = 'proboscis_monkey';
    $LatinToCommon{'oryctolagus_cuniculus'} = 'rabbit';
    $LatinToCommon{'rattus_norvegicus'} = 'rat';
    $LatinToCommon{'macaca_mulatta'} = 'rhesus';
    $LatinToCommon{'procavia_capensis'} = 'rock_hyrax';
    $LatinToCommon{'aplysia_californica'} = 'sea_hare';
    $LatinToCommon{'ovis_aries'} = 'sheep';
    $LatinToCommon{'sorex_araneus'} = 'shrew';
    $LatinToCommon{'choloepus_hoffmanni'} = 'sloth';
    $LatinToCommon{'enhydra_lutris_nereis'} = 'southern_sea_otter';
    $LatinToCommon{'spermophilus_tridecemlineatus'} = 'squirrel';
    $LatinToCommon{'saimiri_boliviensis'} = 'squirrel_monkey';
    $LatinToCommon{'tarsius_syrichta'} = 'tarsier';
    $LatinToCommon{'sarcophilus_harrisii'} = 'tasmanian_devil';
    $LatinToCommon{'echinops_telfairi'} = 'tenrec';
    $LatinToCommon{'nanorana_parkeri'} = 'tibetan_frog';
    $LatinToCommon{'tupaia_belangeri'} = 'tree_shrew';
    $LatinToCommon{'meleagris_gallopavo'} = 'turkey';
    $LatinToCommon{'strongylocentrotus_purpuratus'} = 'urchin';
    $LatinToCommon{'notamacropus_eugenii'} = 'wallaby';
    $LatinToCommon{'ceratotherium_simum'} = 'white_rhinoceros';
    $LatinToCommon{'xenopus_tropicalis'} = 'x._tropicalis';
    $LatinToCommon{'taeniopygia_guttata'} = 'zebra_finch';
    $LatinToCommon{'danio_rerio'} = 'zebrafish';
    
    return \%LatinToCommon;

}







#################################################################
#
#  FUNCTION:  GetStaticLatinTree / GetStaticCommonTree
#
sub GetStaticLatinTree
{
    return '((strongylocentrotus_purpuratus,aplysia_californica),(((((((alligator_mississippiensis,chrysemys_picta_bellii),latimeria_chalumnae),callorhinchus_milii),petromyzon_marinus),ciona_intestinalis),(((((thamnophis_sirtalis,anolis_carolinensis),((((((((((((galeopterus_variegatus,manis_pentadactyla),trichechus_manatus_latirostris),tupaia_belangeri),procavia_capensis),loxodonta_africana),echinops_telfairi),dasypus_novemcinctus),choloepus_hoffmanni),notamacropus_eugenii),sarcophilus_harrisii),ornithorhynchus_anatinus),((((((((pteropus_vampyrus,myotis_lucifugus),(((((ictidomys_tridecemlineatus,heterocephalus_glaber),cavia_porcellus),cricetulus_griseus),dipodomys_ordii),(rattus_norvegicus,mus_musculus))),(oryctolagus_cuniculus,ochotona_princeps)),(ceratotherium_simum,equus_caballus)),((((balaenoptera_acutorostrata_scammoni,vicugna_pacos),sus_scrofa),tursiops_truncatus),((bison_bison_bison,ovis_aries),bos_taurus))),((((neomonachus_schauinslandi,felis_catus),ailuropoda_melanoleuca),canis_lupus_familiaris),(enhydra_lutris_nereis,mustela_putorius_furo))),((((carlito_syrichta,nomascus_leucogenys),otolemur_garnettii),microcebus_murinus),(((((homo_sapiens,pongo_abelii),gorilla_gorilla_gorilla),(pan_troglodytes,pan_paniscus)),((((rhinopithecus_roxellana,chlorocebus_sabaeus),nasalis_larvatus),papio_anubis),(macaca_mulatta,macaca_fascicularis))),(saimiri_boliviensis,callithrix_jacchus)))),(sorex_araneus,erinaceus_europaeus)))),(((aquila_chrysaetos_canadensis,melopsittacus_undulatus),apteryx_australis),((taeniopygia_guttata,geospiza_fortis),(meleagris_gallopavo,gallus_gallus)))),(nanorana_parkeri,xenopus_tropicalis)),((((takifugu_rubripes,oreochromis_niloticus),oryzias_latipes),gadus_morhua),danio_rerio))),(drosophila_melanogaster,anopheles_gambiae)))';
}
sub GetStaticCommonTree
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














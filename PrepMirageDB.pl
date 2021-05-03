#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;


# We'll need bureaucratic powers
sub GetSrcDir { my $lib = $0; $lib =~ s/PrepMirageDB.pl$/src/; return $lib; }
use lib GetSrcDir();
use BureaucracyMirage;


# Additional subroutines
sub ReplaceIllegalChars;
sub ParseUniProt;


# Give usage instructions if it looks like we can be helpful
if (@ARGV != 1) { die "\n  USAGE:  ./PrepMirageDB.pl [protein-database.fa]\n\n"; }


# Grab the name of the input database, since that'll guide all of
# our naming of other sequence files.
my $infname = $ARGV[0];

# Split the input file name into the pre-extension and the extension
$infname =~ /^(\S+)\.([^\.]+)$/;
my $infname_base = $1;
my $infname_ext  = $2;

# Before we create an output file, make sure that we can actually
# access the input file...
my $inf = OpenInputFile($infname);

# Create an output file
my $outfname = $infname_base.'.mirage-ready.'.$infname_ext;
if (-e $outfname) { die "\n  Mirage-ready database '$outfname' already exists\n\n"; }
my $outf = OpenOutputFile($outfname);

# We'll also want a file to track name changes
my $namechangefname = $infname_base.'.name-changes.out';
open(my $namechangef,'>',$namechangefname);

# If we're getting non-Miragey sequence names, we'll want to be able to generate our
# own ID numbers
my %NameUniquenessCheck;

# Generate our Mirage-ready version of the protein database
while (my $line = <$inf>) {

    # Clean it up
    $line =~ s/\n|\r//g;

    # Is this a header line?
    if ($line =~ /^\>(.*)$/) {

	my $full_header = $1;

	# If there are any comments (a '#' followed by whatever) clear 'em out
	my $orig_name = $full_header;
	my $comments  = '';
	if ($full_header =~ /([^\#]+)(\#.*)/) {
	    $orig_name = $1;
	    $comments  = $2;
	    $orig_name =~ s/\s+$//g;
	}

	# Next up, check the format of the sequence name.
	# We'll allow for three formats:
	#
	#   1. species|gene(s)|id
	#   2. More English-y version of (1.)
	#   3. UniProt
	#
	if ($orig_name =~ /^([^\|]+)\|([^\|]+)\|([^\|]+)$/) { ################# 1

	    # The recommended format is being followed, so let's see if they nailed it!
	    my $species = $1;
	    my $input_genes = $2;
	    my $id = $3;

	    $species = ReplaceIllegalChars($species);

	    my $genes = '';
	    foreach my $listed_gene (split(/\//,$input_genes)) {
		my $adjusted_gene = ReplaceIllegalChars($listed_gene);
		if ($genes) {
		    $genes = $genes.'/'.$adjusted_gene;
		} else {
		    $genes = $adjusted_gene;
		}
	    }

	    # Last check -- do we have a unique id?
	    my $attempt_num = 0;
	    my $adj_id = $id;
	    while ($NameUniquenessCheck{$species.'|'.$genes.'|'.$adj_id}) {
		$attempt_num++;
		$adj_id = $id.'.'.$attempt_num;
	    }
	    $NameUniquenessCheck{$species.'|'.$genes.'|'.$adj_id} = 1;

	    # If we changed the name, report it -- otherwise, cruise on along!
	    my $new_name = $species.'|'.$genes.'|'.$adj_id;
	    if ($new_name ne $orig_name) {
		print $namechangef "$orig_name ===[changed-to]==> $new_name\n";
	    }
	    
	    print $outf "\>$new_name";
	    if ($comments) {
		print $outf " $comments";
	    }
	    print $outf "\n";
	    
	} elsif (lc($orig_name) =~ /species\:/) { ################################## 2

	    # Looks like it's the more verbose version of the Mirage name format,
	    # so let's see what we can get
	    my $lc_line = lc($line);
	    $lc_line =~ /species\:\s*(\S+)/;
	    my $species = $1;
	    $species = ReplaceIllegalChars($species);

	    my $input_genes;
	    if ($lc_line =~ /genes?\:\s*(\S+)/) {
		$input_genes = $1;
	    } else {
		close($outf); RunSystemCommand("rm \"$outfname\"");
		close($namechangef); RunSystemCommand("rm \"$namechangefname\"");
		die "\n  Failed to determine gene (gene:[name]) in '$line'\n\n";
	    }

	    my $genes = '';
	    foreach my $listed_gene (split(/\//,$input_genes)) {
		my $adjusted_gene = ReplaceIllegalChars($listed_gene);
		if ($genes) {
		    $genes = $genes.'/'.$adjusted_gene;
		} else {
		    $genes = $adjusted_gene;
		}
	    }

	    my $id;
	    if ($lc_line =~ /id\:\s*(\S+)/) {
		$id = $1;
	    } else {
		close($outf); RunSystemCommand("rm \"$outfname\"");
		close($namechangef); RunSystemCommand("rm \"$namechangefname\"");
		die "\n  Failed to determine id (id:[name]) in '$line'\n\n";
	    }

	    my $attempt_num = 0;
	    my $adj_id = $id;
	    while ($NameUniquenessCheck{$species.'|'.$genes.'|'.$adj_id}) {
		$attempt_num++;
		$adj_id = $id.'.'.$attempt_num;
		$adj_id =~ s/\s//g;
	    }
	    $NameUniquenessCheck{$species.'|'.$genes.'|'.$adj_id} = 1;

	    # Nailed it!
	    my $new_name = $species.'|'.$genes.'|'.$adj_id;

	    # Print it out!
	    print $namechangef "$orig_name ===[changed-to]==> $new_name\n";
	    print $outf "\>$new_name";
	    if ($comments) {
		print $outf " $comments";
	    }
	    print $outf "\n";

	} else { ################################################################### 3

	    my ($species,$genes,$id,$new_comments) = ParseUniProt($orig_name);

	    # If we already have comments, add the 'new' comments to the end
	    if ($comments && $new_comments) {
		$new_comments =~ s/^\#/ /;
		$comments = $comments.$new_comments;
	    }

	    # Can we use the recommended id? (if there is one...)
	    my $adj_id;
	    if ($id) {
		my $attempt_num = 0;
		my $adj_id = $id;
		while ($NameUniquenessCheck{$species.'|'.$genes.'|'.$adj_id}) {
		    $attempt_num++;
		    $adj_id = $id.'.'.$attempt_num;
		}
		$NameUniquenessCheck{$species.'|'.$genes.'|'.$adj_id} = 1;
	    } else {
		my $adj_id = 1;
		while ($NameUniquenessCheck{$species.'|'.$genes.'|'.$adj_id}) {
		    $adj_id++;
		}
	    }
	    
	    my $new_name = $species.'|'.$genes.'|'.$adj_id;

	    print $namechangef "$orig_name ===[changed-to]==> $new_name\n";
	    print $outf "\>$new_name";
	    if ($comments) {
		print $outf " $comments";
	    }
	    print $outf "\n";		
	    
	}

    } else {

	# Sequence line: copy-paste
	print $outf "$line\n";

    }

}

# Close the output files
close($inf);
close($outf);
close($namechangef);

# If there's no difference between the input file and the output file,
# then our work was unnecessary (but it's always empowering to know that
# you started off with a beautifully-formatted database!)
my $diff = OpenSystemCommand("diff \"$infname\" \"$outfname\" \| wc -l");
my $num_diff_lines = <$diff>;
close($diff);
if ($num_diff_lines =~ /^\s*0\s*$/) {

    RunSystemCommand("rm \"$outfname\"");
    RunSystemCommand("rm \"namechangefname\"");
    print "\n  Input sequence file is formatted beautifully and ready for Mirage!\n\n";

} else {

    # We must have seen some differences, so let the user know what's up
    print "\n";
    print "  Database preparation complete!\n";
    print "\n";
    print "  Mirage-ready sequence database: $outfname\n";
    print "  List of sequence name changes : $namechangefname\n";
    print "\n";

}


1; # END OF SCRIPT


###################
#                 #
#   SUBROUTINES   #
#                 #
###################



#############################################################################
#
#  Function: ReplaceIllegalChars
#
#  NOTE: The focus of this function is to remove characters that would
#        interfere with either shell commands or regex, so we're looking
#        to swap out:
#
#             + whitespace
#             + [ ]
#             + | & * $
#
#        with '_'
#
sub ReplaceIllegalChars
{
    my $string = shift;
    $string =~ s/\s|\||\&|\*|\$/\_/g;
    $string =~ s/\[/\(/g;
    $string =~ s/\]/\)/g;
    return $string;
}


#############################################################################
#
#  Function: ParseUniProt
#
sub ParseUniProt
{

    my $orig_name = shift;

    my @EvidenceMeanings;
    $EvidenceMeanings[0] = 'Protein-level-evidence(1)';
    $EvidenceMeanings[1] = 'Transcript-level-evidence(2)';
    $EvidenceMeanings[2] = 'Homology-inferred(3)';
    $EvidenceMeanings[3] = 'Predicted(4)';
    $EvidenceMeanings[4] = 'Uncertain(5)';
    $EvidenceMeanings[5] = 'No-evidence-indicator';

    
    my $species;
    if ($orig_name =~ /OS\=([^\=]+\=?)/) {
	$species = lc($1);
	$species =~ s/\s+\S\S\=$//;
	$species = ReplaceIllegalChars($species);
    } else {
	close($outf); RunSystemCommand("rm \"$outfname\"");
	close($namechangef); RunSystemCommand("rm \"$namechangefname\"");
	die "\n  Failed to determine species (OS=[name]) in '$orig_name'\n\n";
    }
    
    # Let's see if it's pulled straight from UniProt (and thus full of juicy tidbits)
    my $source_db;
    my $accession;
    my $long_gene;
    my $recommended_id = 0;
    if ($orig_name =~ /\>([^\|]+)\|([^\|]+)\|([^\|]+)/) {
	$source_db = 'Database:'.$1;
	$accession = 'Accession:'.$2;
	$long_gene = lc($3);
	$long_gene = ReplaceIllegalChars($long_gene);
	if ($accession =~ /\:[^\-]+\-(\d+)/) {
	    $recommended_id = $1;
	}
    }

    my $gene = '';
    if ($orig_name =~ /GN\=([^\=]+\=?)/) {
	$gene = lc($1);
	$gene =~ s/\s+\S\S\=$//;
	$gene = ReplaceIllegalChars($gene);
	if ($long_gene && $long_gene ne $gene) {
	    $gene = $gene.'/'.$long_gene;
	}
    } elsif ($long_gene) {
	$gene = $long_gene;
    } else {
	close($outf); RunSystemCommand("rm \"$outfname\"");
	close($namechangef); RunSystemCommand("rm \"$namechangefname\"");
	die "\n  Failed to determine gene (GN=[name]) in '$orig_name'\n\n";
    }
    
    my $evidence = 6; # No Entry
    if ($orig_name =~ /PE\=(\d)/) {
	$evidence = $1;
    }

    my $new_comments = '';
    if ($source_db) { $new_comments = $new_comments.' '.$source_db; }
    if ($accession) { $new_comments = $new_comments.' '.$accession; }
    if ($new_comments || $evidence != 6) {
	$evidence = $EvidenceMeanings[$evidence-1];
	$new_comments = $new_comments.' '.$evidence;
    }
    if ($new_comments) { $new_comments = '#'.$new_comments; }

    return ($species,$gene,$recommended_id,$new_comments);
    
}


# EOF



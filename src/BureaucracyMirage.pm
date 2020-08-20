#
#  BureaucracyMirage.pm : File I/O-type stuff for Mirage 'n' Pals
#
use warnings;
use strict;
use POSIX;

# Very general functions
sub Max;
sub Min;
sub RunSystemCommand;
sub OpenSystemCommand;
sub SpawnProcesses;
sub ConfirmFile;
sub OpenInputFile;
sub NameOutputFile;
sub OpenOutputFile;
sub ConfirmDirectory;
sub OpenDirectory;
sub CreateDirectory;

# More bio-specific functions
sub ParseMirageSeqName;
sub TranslateCodon;

# Global codon table
#my @CODONTABLE =
#('K','N','K','N',  # AA-
#'T','T','T','T',  # AC-
#'R','S','R','S',  # AG-
#'I','I','M','I',  # AT-
#'Q','H','Q','H',  # CA-
#'P','P','P','P',  # CC-
#'R','R','R','R',  # CG-
#'L','L','L','L',  # CT-
#'E','D','E','D',  # GA-
#'A','A','A','A',  # GC-
#'G','G','G','G',  # GG-
#'V','V','V','V',  # GT-
#'X','Y','X','Y',  # TA-
#'S','S','S','S',  # TC-
#'S','C','W','C',  # TG-
#'L','F','L','F'); # TT-


# TRANSLATION HASH
my %CODONHASH = (
    'AAA','K',  'AAC','N',  'AAG','K',  'AAT','N',
    'ACA','T',  'ACC','T',  'ACG','T',  'ACT','T',
    'AGA','R',  'AGC','S',  'AGG','R',  'AGT','S',
    'ATA','I',  'ATC','I',  'ATG','M',  'ATT','I',
    'CAA','Q',  'CAC','H',  'CAG','Q',  'CAT','H',
    'CCA','P',  'CCC','P',  'CCG','P',  'CCT','P',
    'CGA','R',  'CGC','R',  'CGG','R',  'CGT','R',
    'CTA','L',  'CTC','L',  'CTG','L',  'CTT','L',
    'GAA','E',  'GAC','D',  'GAG','E',  'GAT','D',
    'GCA','A',  'GCC','A',  'GCG','A',  'GCT','A',
    'GGA','G',  'GGC','G',  'GGG','G',  'GGT','G',
    'GTA','V',  'GTC','V',  'GTG','V',  'GTT','V',
    'TAA','X',  'TAC','Y',  'TAG','X',  'TAT','Y',
    'TCA','S',  'TCC','S',  'TCG','S',  'TCT','S',
    'TGA','S',  'TGC','C',  'TGG','W',  'TGT','C',
    'TTA','L',  'TTC','F',  'TTG','L',  'TTT','F');


###################################################################
#
#  FUNCTION:  Max
#
sub Max
{
    my $a = shift;
    my $b = shift;
    return $a if ($a > $b);
    return $b;
}


###################################################################
#
#  FUNCTION:  Min
#
sub Min
{
    my $a = shift;
    my $b = shift;
    return $a if ($a < $b);
    return $b;
}


###################################################################
#
#  FUNCTION:  RunSystemCommand
#
sub RunSystemCommand
{
    my $command = shift;
    if (system($command)) { die "\n  ERROR:  System command '$command' failed\n\n"; }
}


###################################################################
#
#  FUNCTION:  OpenSystemCommand
#
sub OpenSystemCommand
{
    my $command = shift;
    if ($command !~ / \|$/) { $command = $command.' |'; }
    open(my $commandf,$command) || die "\n  ERROR:  Failed to open results from system command '$command'\n\n";
    return $commandf;
}


###################################################################
#
#  FUNCTION:  SpawnProcesses
#
sub SpawnProcesses
{
    my $requested_processes = shift;

    my $num_processes = 1;
    my $threadID = 0;
    my $pid;
    while ($num_processes < $requested_processes) {
	if ($pid = fork) {
	    if (not defined $pid) { die "\n  ERROR:  Fork failed\n\n"; }
	    $num_processes++;
	} else {
	    $threadID = $num_processes;
	    last;
	}
    }
    return $threadID;
}




###################################################################
#
#  FUNCTION:  ConfirmFile
#
sub ConfirmFile
{
    my $fname = shift;
    if (!(-e $fname)) { die "\n  ERROR:  Failed to locate file '$fname'\n\n"; }
    return $fname;
}


###################################################################
#
#  FUNCTION:  OpenInputFile
#
sub OpenInputFile
{
    my $fname = shift;
    ConfirmFile($fname);
    open(my $inf,'<',$fname) || die "\n  ERROR:  Failed to open input file '$fname'\n\n";
    return $inf;
}


###################################################################
#
#  FUNCTION:  NameOutputFile
#
sub NameOutputFile
{
    my $fname = shift;

    my $extension = '';
    if ($fname =~ /(\.[^\.]+)$/) {
	$extension = $1;
	$fname =~ s/\.[^\.]+$//;
    }

    my $outfname = $fname.$extension;
    my $fcount = 1;
    while (-e $outfname) {
	$fcount++;
	$outfname = $fname.'-'.$fcount.$extension;
    }

    return $outfname;
}


###################################################################
#
#  FUNCTION:  OpenOutputFile
#
sub OpenOutputFile
{
    my $fname = shift;
    $fname = NameOutputFile($fname);
    open(my $outf,'>',$fname) || die "\n  ERROR:  Failed to open output file '$fname'\n\n";
    return $outf;
}


###################################################################
#
#  FUNCTION:  ConfirmDirectory
#
sub ConfirmDirectory
{
    my $dirname = shift;
    if ($dirname !~ /\/$/) { $dirname = $dirname.'/'; }
    if (!(-d $dirname)) { die "\n  ERROR:  Failed to locate directory '$dirname'\n\n"; }
    return $dirname;
}


###################################################################
#
#  FUNCTION:  OpenDirectory
#
sub OpenDirectory
{
    my $dirname = shift;
    $dirname = ConfirmDirectory($dirname);
    opendir(my $dir,$dirname) || die "\n  ERROR:  Failed to open directory '$dirname'\n\n";
    return $dir;
}


###################################################################
#
#  FUNCTION:  CreateDirectory
#
sub CreateDirectory
{
    my $dirname = shift;

    my $dirbase = $dirname;
    if ($dirbase =~ /\/$/) {
	$dirbase =~ s/\/$//;
    }

    my $dircount = 1;
    while (-d $dirname) {
	$dircount++;
	$dirname = $dirbase.'-'.$dircount;
    }
    $dirname = $dirname.'/';

    if (system("mkdir \"$dirname\"")) { die "\n  ERROR:  Failed to create directory '$dirname'\n\n"; }
    return ConfirmDirectory($dirname);
}


###################################################################
#
#  BIO. FUNCTION:  ParseMirageSeqName
#
#  INPUT  : seqname
#  OUTPUT : $species,$gene,$iso_id,$mirage_id,$extra_info
#
sub ParseMirageSeqname
{
    my $seqname = shift;
    $seqname =~ /^\>?([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)/;
    my $species   = $1;
    my $gene      = $2;
    my $iso_id    = $3;
    my $mirage_id = $4;

    my $extra_info = 0;
    if ($seqname =~ /\#(.+)$/) {
	$extra_info = $1;
    }

    return($species,$gene,$iso_id,$mirage_id,$extra_info);
    
}



###################################################################
#
#  BIO. FUNCTION:  TranslateCodon
#
#  INPUT : $codon_str
#  OUTPUT: $amino
#
#  NOTE  :  We'll return '-' for abnormal inputs and 'X' for stop codons
#
sub TranslateCodon
{
    my $codon_str = shift;
    $codon_str = uc($codon_str);
    $codon_str =~ s/U/T/g; # DNA
    if (length($codon_str) != 3 || $codon_str =~ /[^ACGT]/) { return '-'; }
    return $CODONHASH{$codon_str};
}



1; # EOF

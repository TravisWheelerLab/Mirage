#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;
use Getopt::Long;
use Time::HiRes;
use File::Basename;
use lib dirname (__FILE__);
use Cwd;

# YUCKITY YUCK YUCK
sub GetThisDir { my $lib = $0; $lib =~ s/\/Quilter2.pl$//; return $lib; }
use lib GetThisDir();
use BureaucracyMirage;
use DisplayProgress;

sub PrintUsage;
sub ParseArgs;
sub ParseChromSizes;
sub UseGTF;
sub ParseGTF;
sub CheckGeneAliases;
sub UseFastMap;
sub ParseFastMapOutput;
sub GenSpliceGraph;
sub PrintHitToWeaverInf;
sub GetSpliceStrengths;
sub CheckForFullMaps;
sub ExtractCanonXWMap;
sub AttemptSpalnFill;
sub AttemptChimericXWMap;
sub ExtractXWInputExons;
sub ExtractChimericXWMap;
sub FindMaximalXWHitSet;
sub SaveTopXWInputs;
sub GetAminoHitGaps;
sub RecordMaximalHits;
sub RunBlatOnFileSet;
sub GenBlatMaps;
sub AttemptBlatFill;
sub BlatToSpalnSearch;
sub GetSpalnSfetchRanges;
sub SpalnSearch;
sub ClusterNuclRanges;
sub AdjustNuclCoord;
sub ParseBlatLine;
sub ParseSpalnOutput;
sub FinalFileCheck;

# Original Spaln parsing code
sub OPSO_ParseSpalnOutput;
sub OPSO_SelenocysteineCheck;


# Help?
if (@ARGV < 3) { PrintUsage(0); }



# As you'd expect, we'll parse the commandline arguments right up top
my $opts_ref = ParseArgs();
my %Opts     = %{$opts_ref};

# As you'd expect, we want to know what there is to know about our proteins
my $species_dirname = ConfirmDirectory($ARGV[0]);
my $seq_dirname = ConfirmDirectory($species_dirname.'seqs');

# As you'd expect, next up is confirming the genome
my $genome = ConfirmFile($ARGV[1]);

# As you'd expect, we'll want to parse the gtf, if one's been provided
my $gtfname = $ARGV[2];

# As long as everything's going as expected, figure out where we are
my $srcdir = $0;
$srcdir =~ s/Quilter2.pl$//;
my $sfetch = $srcdir.'../inc/hsi/sfetch';
my $spaln  = $srcdir.'../inc/spaln2.3.3/src/spaln';

# Spaln requires certain environment variables to be set, so why don'tcha set 'em?
$spaln =~ /^(.*)spaln$/;
$ENV{'ALN_TAB'} = $1.'../table';
$ENV{'ALN_DBS'} = $1.'../seqdb';

# Since we use the same options for all spaln searches, let's just set those now
# NOTE: I'm going to test out forcing forward-strand only (S1)
#       as opposed to both strands (S3).  I'm also boosting the weight
#       for coding potential (-yz=2 -> -yz=4) and lowering the weight
#       for splice site signal (-yy=8 -> -yy=4).
my $spaln_opts = ' -Q3 -O1 -S1 -ya3 -yz4 -yy4 ';
$spaln = $spaln.$spaln_opts;

# How many CPUs do we intend to use?
my $ThreadGuide = OpenInputFile($seq_dirname.'Thread-Guide');
my $num_cpus = <$ThreadGuide>;
close($ThreadGuide);
$num_cpus =~ /Num CPUs\: (\d+)/;
$num_cpus = $1;

# We'll also read in our chromosome sizes early, since that info can
# be useful in various places
my $chrsizes_ref = ParseChromSizes($genome);
my %ChrSizes = %{$chrsizes_ref};

# Initialize the progress displaying variables
InitQuilterProgressVars($seq_dirname,$num_cpus);

# If we're using a GTF for our first round of searching, hop to it!
# Otherwise, jump straight to building up one big file for BLAT search.
my $blat_naming;
my @BlatFileNames;
if ($gtfname ne '-') {

    UseGTF($seq_dirname,$genome,$gtfname,$num_cpus,\%Opts);

    # Each process MIGHT have left us with a blat file if it had
    # any unfinishable business, so we'll want to aggregate those
    for (my $i=0; $i<$num_cpus; $i++) {
	my $blatfname = $seq_dirname.$i.'.blat.fa';
	next if (!(-e $blatfname));
	push(@BlatFileNames,$blatfname);
    }

    # We're running our BLAT search using the sequence names as the way
    # to derive gene family membership, and you're powerless to stop it.
    $blat_naming = 'seqname';

} else {

    # We'll need to put all of the files that we want to use in our
    # BLAT search into a list, because that's what I've decided we
    # need to do.
    my $SeqDir = OpenDirectory($seq_dirname);
    while (my $fname = readdir($SeqDir)) {
	next if ($fname !~ /\.fa$/);
	push(@BlatFileNames,$seq_dirname.$fname);
    }
    closedir($SeqDir);

    # We're running our BLAT search using file names as the way to derive
    # the gene families that sequences belong to, and we're doing it now!
    $blat_naming = 'filename';

}


# If we have any sequences to run BLAT on, we'd better get to it!
if (scalar(@BlatFileNames)) {
    
    my ($blat_outfname,$blat_nameguide_ref,$blat_genes_ref)
	= RunBlatOnFileSet(\@BlatFileNames,$genome,$blat_naming);
    my %BlatNameGuide = %{$blat_nameguide_ref};
    my %BlatGenes = %{$blat_genes_ref};

    # Now we've run BLAT -- time to make sense of all the cool things it's
    # telling us!
    GenBlatMaps($blat_outfname,\%BlatNameGuide,\%BlatGenes,$num_cpus);

    # Cleanup on line 160! << CRITICAL: KEEP THIS LINE NUMBER CORRECT, FOR JOKE
    RunSystemCommand("rm \"$blat_outfname\"");

}


# We'll depart on a final endeavor to make sure our ".quilter.out" files are
# well-organized (i.e., communicate canonical chromosome and have the sequences
# sorted)
FinalFileCheck($seq_dirname,$num_cpus);



1;


###################
#  END OF SCRIPT  #
###################





############################################################
#
#  Function: PrintUsage
#
sub PrintUsage
{
    my $help_requested = shift;
    print "\n  USAGE:  ./Quilter2.pl [species-dir] [genome.fa] [.gtf] {OPT.s}\n";
    die "\n" if (!$help_requested);
    die "\n"; # Eventually we'll give more info if help is requestied...
}




############################################################
#
#  Function: ParseArgs
#
sub ParseArgs
{
    my %Options;
    &GetOptions(
	\%Options,
	"help",
	"v",
	"time",
	) || die "\n  ERROR:  Failed to parse Quilter2 commandline arguments\n\n";

    if ($Options{help}) { PrintUsage(1); }
    
    return \%Options;
}





############################################################
#
#  Function: ParseChromSizes
#
sub ParseChromSizes
{
    my $genome = shift;

    my $chrsize_fname = $genome;
    $chrsize_fname =~ s/\.fa$/\.chrom\.sizes/;

    my $ChrSizeFile = OpenInputFile($chrsize_fname);
    my %ChrSizes;
    while (my $line = <$ChrSizeFile>) {
	if ($line =~ /^(\S+)\s+(\d+)$/) {
	    my $chr = $1;
	    my $size = $2;
	    $ChrSizes{$chr} = $size;
	    $ChrSizes{$chr.'[revcomp]'} = $size;
	}
    }
    close($ChrSizeFile);

    return \%ChrSizes;
}





############################################################
#
#  Function: UseGTF
#
sub UseGTF
{
    my $seq_dirname = shift;
    my $genome      = shift;
    my $gtfname     = shift;
    my $num_cpus    = shift;
    my $opts_ref    = shift;
    my %Opts = %{$opts_ref};

    # Get the name of the species from the directory name
    $seq_dirname =~ /\/([^\/]+)\/seqs\/$/;
    my $species = $1;

    # UHHHHHH, how exactly did you plan on using a GTF without PARSING IT FIRST?!
    my $gtf_ref = ParseGTF($gtfname,$species,$seq_dirname.'../../../gene-aliases');
    my %GTF = %{$gtf_ref};

    # Go ahead and spin up them there CPUs!
    my $threadID = SpawnProcesses($num_cpus);
    
    # Re-open the thread guide and scan to your assignment region
    my $ThreadGuide = OpenInputFile($seq_dirname.'Thread-Guide');
    my $gene;
    while (my $line = <$ThreadGuide>) {
	if ($line =~ /^$threadID (\S+)$/) {
	    $gene = $1;
	    last;
	}
    }

    # A file where we can write out sequences that we need to BLAT.
    # Because these files end up being appended together by the main
    # thread, we can leave setting up naming conventions until that time.
    my $blatfname = $seq_dirname.$threadID.'.blat.fa';
    my $BlatFile  = OpenOutputFile($blatfname);

    # Work it, thready!
    my $genes_completed = 0;
    while ($gene) {

	# Let's name the protein file we'll be working with
	my $gene_fname = $seq_dirname.$gene.'.fa';

	# Try to build a splicegraph-like thing using FastMap2 and our GTF data
	UseFastMap($gene_fname,$genome,\%GTF,$BlatFile);

	# Grab the next family
	last if (eof($ThreadGuide));
	$gene = <$ThreadGuide>;
	if ($gene =~ /^$threadID (\S+)$/) { $gene = $1; }
	else                              { $gene = 0;  }

	# Report your progress
	$genes_completed++;
	DispProgQuilter('fm2|'.$threadID.'|'.$genes_completed);

	
    }

    # Close up the blat file, and (if we're so blessed) clear it if it's empty
    close($BlatFile);
    if (!(-s $blatfname)) {
	RunSystemCommand("rm \"$blatfname\"") if (-e $blatfname);
    }
    
    # Wait for everyone to finish
    if ($threadID) { exit(0); }
    while (wait() != -1) {}

    # GTF used!
    
}





############################################################
#
#  Function: ParseGTF
#
sub ParseGTF
{
    my $gtfname    = shift;
    my $species    = shift;
    my $aliasfname = shift;

    DispProgQuilter('parsing-gtf');

    # First off, if there are any aliases, let's take notice!
    my %GeneAliases;
    if (-e $aliasfname) {
	my $gene_aliases_ref = CheckGeneAliases($aliasfname);
	%GeneAliases = %{$gene_aliases_ref};
    }
    
    my $GTFile = OpenInputFile($gtfname);

    # Essentially, what we want to do is pull all of the GTF entries for our
    # species into a hash
    my %GTF;
    my %DuplicateCheck;
    while (my $line = <$GTFile>) {

	my $entry_type = '';
	if ($line =~ /^\S+\s+\S+\s+(\S+)/) {
	    $entry_type = lc($1);
	}
	
        if ($line =~ /^\#/ || ($entry_type ne 'exon' && $entry_type ne 'cds')) {
	    next;
	}

	# Note that startPos is always less than endPos, even when we're
        # indexing into the reverse complement.
        $line =~ /^(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+([\+\-])/;
        my $chr     = $1;
        my $start   = $2;
        my $end     = $3;
        my $revcomp = $4;

	# If this chromosome doesn't appear in our genome, we can't use
	# this exon data
	#
	# NOTE that a common mismatch between genomes and GTFs is that
	# genomes will name their chromosomes like 'chr18', where the
	# GTF would just use '18' so we'll see if this is info. worth
	# noting before nexting.
	#
	if (!$ChrSizes{$chr}) {
	    if ($ChrSizes{'chr'.$chr}) {
		$chr = 'chr'.$chr;
	    } else {
		next;
	    }
	}
	
        if ($revcomp eq '-') {
	    $chr = $chr.'[revcomp]';
            my $temp = $start;
            $start = $end;
            $end = $temp;
	}

        $line =~ /gene\_name \"([^\"]+)\"\;/;
        my $gene  = lc($1);
	my $entry = $chr.':'.$start.'..'.$end;

	next if ($DuplicateCheck{$entry});
	$DuplicateCheck{$entry} = 1;
	
	if ($GTF{$gene}) { $GTF{$gene} = $GTF{$gene}.'#'.$entry; }
	else             { $GTF{$gene} = $entry;                 }

	# If this gene has any aliases, record to them, too!
	next if (!$GeneAliases{$gene});
	foreach my $alias (split(/\//,$GeneAliases{$gene})) {
	    if ($GTF{$alias}) { $GTF{$alias} = $GTF{$alias}.'#'.$entry; }
	    else              { $GTF{$alias} = $entry;                  }
	}

    }
    close($GTFile);

    # We'll do a quick confirmation that we found anything of value in the
    # GTF -- if we didn't, that's a warning that the genome and GTF don't
    # agree on sequence names
    my @AllGTFGenes = keys %GTF;
    if (scalar(@AllGTFGenes) == 0) {
	print "\n";
	print "  WARNING:  No usable exon entries found in GTF.\n";
	print "            Recommend checking that chromosome names match genome.\n";
	print "\n";
    }

    return \%GTF;
    
}





############################################################
#
#  Function: CheckGeneAliases
#
sub CheckGeneAliases
{
    my $fname = shift;

    my %Aliases;
    my $inf = OpenInputFile($fname);
    while (my $line = <$inf>) {

	my @Genes = split(/\s/,$line);
	my $gene  = $Genes[0];
	for (my $i=1; $i<scalar(@Genes); $i++) {
	    my $alias = $Genes[$i];
	    if ($Aliases{$alias}) { $Aliases{$alias} = $Aliases{$alias}.'/'.$gene; }
	    else                  { $Aliases{$alias} = $gene;                      }
	}

    }
    close($inf);

    return \%Aliases;
}





############################################################
#
#  Function: UseFastMap
#
sub UseFastMap
{
    my $gene_fname = shift;
    my $genome     = shift;
    my $gtf_ref    = shift;
    my $BlatFile   = shift;

    my %GTF = %{$gtf_ref};

    $gene_fname =~ /\/([^\/]+)\.fa$/;
    my $gene = $1;

    # Before we get into the mapping business, let's load in the protein sequences,
    # since they might have something to say about how we build our graph...
    my $GeneFile = OpenInputFile($gene_fname);
    my @Seqs;
    my @SeqNames;
    my $num_seqs = -1;
    while (my $line = <$GeneFile>) {
	$line =~ s/\n|\r//g;
	if ($line =~ /\>(\S+)/) {
	    $SeqNames[++$num_seqs] = $1;
	    $Seqs[$num_seqs] = '';
	} else {
	    $Seqs[$num_seqs] = $Seqs[$num_seqs].uc($line);
	}
    }
    $num_seqs++;
    close($GeneFile);

    # If there isn't an entry for this gene, fuhget about it!
    # And, obviously, by 'fuhget about it' I mean 'add it to the Blat file'
    if (!$GTF{$gene}) {
	for (my $i=0; $i<$num_seqs; $i++) {
	    print $BlatFile ">$gene\|$SeqNames[$i]\n";
	    my @Seq = split(//,$Seqs[$i]);
	    for (my $j=0; $j<scalar(@Seq); $j++) {
		print $BlatFile "$Seq[$j]";
		print $BlatFile "\n" if (($j+1) % 60 == 0);
	    }
	    print $BlatFile "\n" if (scalar(@Seq) % 60);
	    print $BlatFile "\n";
	}
	return;
    }

    # Break it down by chromosome (in case this gene is associated with multiple)    
    my %RangesByChr;
    foreach my $chr_range (split(/\#/,$GTF{$gene})) {

	$chr_range =~ /^(\S+)\:(\d+\.\.\d+)$/;
	my $chr = $1;
	my $range = $2;

	if ($RangesByChr{$chr}) { $RangesByChr{$chr} = $RangesByChr{$chr}.'/'.$range; }
	else                    { $RangesByChr{$chr} = $range;                        }
	
    }

    # Name the file we'll write our chromosomal sequences out to
    my $nucl_fname = $gene_fname;
    $nucl_fname =~ s/\.fa$/\.nucl\.fa/;

    # Our goal is to build a splice graph (possibly per chromosome...) for each
    # sequence in this family, so we'll first build up a list of filenames of
    # ExonWeaver outputs.
    my @SpliceGraphs;
    my @XWInputNames;
    for (my $i=0; $i<$num_seqs; $i++) {
	$XWInputNames[$i] = 0;
	$SpliceGraphs[$i] = 0;
    }

    # Now we can go chromosome-by-chromosome and try to build the
    # best splice graph for our sequences we can on a single chromosome
    #
    # While we're doing this, we'll also build up a list of ranges to use
    # in a potential Spaln search.
    my @SpalnStarts;
    my @SpalnEnds;
    my @SpalnChrs;
    foreach my $chr (keys %RangesByChr) {

	my @Ranges = split(/\//,$RangesByChr{$chr});

	# We'll need to put our values into these arrays so we can make sure
	# our Spaln ranges don't end up being too large
	my @RangeStarts;
	my @RangeEnds;

	# Is this your *real* name, MR(S). DR. CHROMOSOME?!
	my $revcomp = 0;
	my $search_start;
	my $search_end;
	my $exon_list_str = '';
	if ($chr =~ /\[revcomp\]/) {

	    $chr =~ s/\[revcomp\]//;
	    $revcomp = 1;

	    $search_start = 0;
	    $search_end   = $ChrSizes{$chr};
	    foreach my $range (@Ranges) {
		
		$range =~ /^(\d+)\.\.(\d+)$/;
		my $high = $1;
		my $low  = $2;
		
		if ($high > $search_start) { $search_start = $high; } 
		if ($low  < $search_end)   { $search_end   = $low;  }
		$exon_list_str = $exon_list_str.' '.$high.' '.$low;

		push(@RangeStarts,$high);
		push(@RangeEnds,$low);
		
	    }

	    # Buff it up! Until you can feel it!
	    $search_start += 1000;
	    $search_end   -= 1000;
	    
	} else {

	    $search_start = $ChrSizes{$chr};
	    $search_end   = 0;
	    foreach my $range (@Ranges) {
		
		$range =~ /^(\d+)\.\.(\d+)$/;
		my $low  = $1;
		my $high = $2;
		
		if ($low  < $search_start) { $search_start = $low;  }
		if ($high > $search_end)   { $search_end   = $high; }
		$exon_list_str = $exon_list_str.' '.$low.' '.$high;

		push(@RangeStarts,$low);
		push(@RangeEnds,$high);
		
	    }

	    # Buff it up! And you don't even need it!
	    $search_start -= 1000;
	    $search_end   += 1000;
	    
	}

	# Get the ranges we'd use for sfetchery if push comes to Spaln
	my ($starts_ref,$ends_ref,$num_spaln_ranges,$sum_len)
	    = GetSpalnSfetchRanges(\@RangeStarts,\@RangeEnds,$ChrSizes{$chr},0);
	@RangeStarts = @{$starts_ref};
	@RangeEnds   = @{$ends_ref};

	# Record the range for this chromosome
	for (my $i=0; $i<$num_spaln_ranges; $i++) {
	    push(@SpalnStarts,$RangeStarts[$i]);
	    push(@SpalnEnds,$RangeEnds[$i]);
	    if ($revcomp) { push(@SpalnChrs,$chr.'-'); }
	    else          { push(@SpalnChrs,$chr.'+'); }
	}

	# Just in case we're pushing up against the edge of the
	# genomic sequence (this is primarily a concern with 'alt'
	# chromosomes), make sure we don't overstep.
	$search_end   = Min($search_end,  $ChrSizes{$chr});
	$search_start = Max($search_start,1);

	# Extract the appropriate genomic region for our search
	RunSystemCommand($sfetch." -range $search_start\.\.$search_end \"$genome\" \"$chr\" > \"$nucl_fname\"");

	# It's... FASTMAP2 TIME!
	my $mapcmd = $srcdir."FastMap2 \"$gene_fname\" $num_seqs \"$nucl_fname\" $exon_list_str";

	
	################################################################
	#
	# It's... JUST GATHERING TIMING DATA!
	#
	#$mapcmd = $mapcmd.' 1>/dev/null';
	#my $starttime = [Time::HiRes::gettimeofday()];
	#RunSystemCommand($mapcmd);
	#my $runtime = Time::HiRes::tv_interval($starttime);
	#return $runtime;
	#
	# It's... A TEST TO SEE IF BLOCK-Y FASTMAP IS TOO SLOW!
	#
	################################################################

	my $hitref = ParseFastMapOutput($mapcmd,$num_seqs);

	# We now have an array with the hits for each sequence stored in
	# an &-separated list.  These hits are in |-separated format,
	# and ordered as follows:
	#
	#   0. amino start position
	#   1. amino end   position
	#   2. nucl  start position
	#   3. nucl  end   position
	#   4. exon number
	#   5. reading frame (0-2)
	#   6. score (BLOSUM 62, half-bit)
	#   7. nucleotide sequence (w/ 17 lowercase buffer nucls on each side)
	#
	# Moreover, they've been run through Perl's sort, so we should be in
	# order of starting position in the amino acid sequence.
	#
	# Our next task is to ID the best splice graph that we can for each of
	# the sequences, individually.

	my @HitsBySeq = @{$hitref};
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($HitsBySeq[$i]) {

		# Get the name of the files provided to and produced by ExonWeaver
		my ($xw_infname,$xw_outfname) =
		    GenSpliceGraph($HitsBySeq[$i],$Seqs[$i],$SeqNames[$i],$chr,$revcomp,$gene_fname);

		# If we didn't get any spliced alignments out of ol' XW, move along
		next if (!$xw_outfname);

		# That's lookin' like a graph to me!  Toss it on the pile!
		if ($SpliceGraphs[$i]) {
		    $XWInputNames[$i] = $XWInputNames[$i].';'.$xw_infname;
		    $SpliceGraphs[$i] = $SpliceGraphs[$i].';'.$xw_outfname;
		} else {
		    $XWInputNames[$i] = $xw_infname;
		    $SpliceGraphs[$i] = $xw_outfname;
		}
	    }
	}
    }

    # We'll save removing the nucleotide file until the end, in case we end up
    # reusing the file name for SPALN-assisted search
    
    
    # SpliceGraphs is now an array of file names (or zeros) for ExonWeaver output,
    # which we can read in to check for full-protein splicings (or else to see if
    # there's at least agreement on some canonical exons...)

    
    # First off, we'll see which sequences have full-protein mappings and
    # record the mapped chromosomes, so we can determine a canonical
    # chromosome
    my @FullMaps;
    my %FullMapsByChr;
    my $num_unmapped = 0;
    for (my $i=0; $i<$num_seqs; $i++) {
	if ($SpliceGraphs[$i]) {

	    $FullMaps[$i] = CheckForFullMaps($SpliceGraphs[$i],length($Seqs[$i]));

	    # If there weren't any full mappings, skip to the next sequence.
	    # Otherwise, get the names of the chromosomes (from the file names)
	    # that had full mappings
	    if (!$FullMaps[$i]) {
		$num_unmapped++;
		next;
	    }

	    foreach my $xw_fname (split(/\;/,$FullMaps[$i])) {

		$xw_fname =~ /\.([^\.]+)\.weaver\.out/;
		my $chr = $1;

		if ($FullMapsByChr{$chr}) { $FullMapsByChr{$chr}++; }
		else                      { $FullMapsByChr{$chr}=1; }

	    }
	    
	} else {	    
	    $FullMaps[$i] = 0;
	    $num_unmapped++;
	}
    }

    # Which chromosome scored the most full maps?
    my $top_chr = 0;
    my $top_chr_mapcount = 0;
    foreach my $chr (keys %FullMapsByChr) {
	if ($FullMapsByChr{$chr} > $top_chr_mapcount) {
	    $top_chr = $chr;
	    $top_chr_mapcount = $FullMapsByChr{$chr};
	}
    }

    # If we don't have any full maps, we'll end up punting to BLAT.
    # Otherwise, let's kick off the mapping parsery!
    if ($top_chr) {

	# Because we're using a filename-friendly way of transcribing
	# the chromosome name, 

	# Sweet!  We've got a canonical chromosome!
	#
	# Our plan now is to separate the sequences between those with
	# full maps to the canonical chromosome and those without.
	#
	# We'll repurpose 'FullMaps' to indicate whether a full map was
	# found to the canonical chromosome, and use a separate array
	# to track the names of those output files.
	#
	for (my $i=0; $i<$num_seqs; $i++) {

	    # Determine which camp we fall in
	    my $top_chr_fname = 0;
	    if ($SpliceGraphs[$i] =~ /([^\;]+\.$top_chr\.weaver\.out)/) {

		$top_chr_fname = $1;

		# Because gene names can have parentheses, we need
		# to do some very dumb explicitification
		$top_chr_fname =~ s/\(/\\\(/g;
		$top_chr_fname =~ s/\)/\\\)/g;
		if ($FullMaps[$i] =~ /$top_chr_fname/) {
		    $FullMaps[$i] = 1;
		} else {
		    $FullMaps[$i] = 0;
		}
		
	    } else {
		$FullMaps[$i] = 0;
	    }


	    # If we have a full mapping, extract it!
	    # Otherwise, we'll check if we can produce a chimeric mapping...
	    my $seq_len = length($Seqs[$i]);
	    if ($FullMaps[$i]) {
		# TOO EASY!!!
		$top_chr_fname =~ s/\\\(/\(/g;
		$top_chr_fname =~ s/\\\)/\)/g;
		$FullMaps[$i] = ExtractCanonXWMap($top_chr_fname,$seq_len,
						  $SeqNames[$i],$top_chr);
	    }
	    next;
	    #
	    #
	    #
	    # * * * * * * * * START UNUSED CODE BLOCK * * * * * * * * *
	    #
	    # NOTE: The above 'next' is essentially a way to block out the below code,
	    #       at least for now.  I might re-activate it at some point, but I
	    #       want to switch over to a different approach (without losing the
	    #       currently available infrastructure)
	    #
	    #} elsif ($top_chr_fname) {
	    if (!$FullMaps[$i] && $top_chr_fname) {
		
		$top_chr_fname =~ s/\\\(/\(/g;
		$top_chr_fname =~ s/\\\)/\)/g;

		# We know this sequence is at least partially hitting to the
		# expected chromosome.  We'll see if we can use SPALN to find
		# putative exons to fill in the gaps in the 
		my ($max_exon_str,$max_amino_str,$max_hits_str) =
		    FindMaximalXWHitSet($SpliceGraphs[$i],$top_chr_fname);

		# We can easily build up a set of ranges of that we'll need to fill in,
		# knowing where we have solid hits to the top chromosome.
		my ($gap_starts_ref,$gap_ends_ref)
		    = GetAminoHitGaps($max_amino_str,length($Seqs[$i]));
		my @GapStarts = @{$gap_starts_ref};
		my @GapEnds   = @{$gap_ends_ref};

		# Before we get too excited about the possibility of trans splicing,
		# let's consider whether the gaps we're seeing are due to unindexed
		# exons that we could detect with Spaln.
		#
		# The biggest nuisance with this is that we're going to need to rip
		# SPALN putative exon ranges, but then plug them into a new ExonWeaver
		# run, since we want clean splice boundary delineation between our
		# exons.
		#
		$FullMaps[$i] =
		    AttemptSpalnFill($top_chr_fname,$max_hits_str,\@GapStarts,\@GapEnds,
				     $Seqs[$i],$SeqNames[$i],$genome,$gene_fname);

		# At this point, if we still don't have a full mapping, it seems like
		# something nonstandard could be going on, so we'll quickly check if
		# the full set of hits (across chromosomes) could be useful...
		if (!$FullMaps[$i]) {
		    $FullMaps[$i] =
			AttemptChimericXWMap($XWInputNames[$i],$SpliceGraphs[$i],
					     $top_chr_fname,$max_exon_str,\@GapStarts,
					     \@GapEnds,$Seqs[$i],$SeqNames[$i],
					     $gene_fname);
		}
		
		# At this point, we either have a full hit or we don't, as is the way
		# of logic.
		#
		# If we don't have a full mapping, we'll do some prep work for BLAT.
		# In addition to writing out the full sequence, we'll also put out
		# the parts of the sequence that correspond to the gaps left by our
		# partially mapping.
		#
		if (!$FullMaps[$i]) {

		    # First, write out the full sequence...
		    print $BlatFile ">$gene\|$SeqNames[$i]\n";
		    my @Seq = split(//,$Seqs[$i]);
		    for (my $j=0; $j<scalar(@Seq); $j++) {
			print $BlatFile "$Seq[$j]";
			print $BlatFile "\n" if (($j+1) % 60 == 0);
		    }
		    print $BlatFile "\n" if (scalar(@Seq) % 60);
		    print $BlatFile "\n";

		    # Next up, write out our Sequence McNuggets
		    for (my $gap_id=0; $gap_id<scalar(@GapStarts)-1; $gap_id++) {
			my $start_amino = $GapStarts[$gap_id];
			my $end_amino = $GapEnds[$gap_id];
			my $part_name = $SeqNames[$i].'s'.$start_amino.'e'.$end_amino;
			print $BlatFile ">$gene\|$part_name\n";
			for (my $j=0; $j<=$end_amino-$start_amino; $j++) {
			    print $BlatFile "$Seq[$j+$start_amino]";
			    print $BlatFile "\n" if (($j+1) % 60 == 0);
			}
			print $BlatFile "\n" if (($end_amino-$start_amino+1) % 60);
			print $BlatFile "\n";
		    }

		    # We'll want to hold onto our maximal hits in a special file
		    # (or, if not special, at least one that won't be destroyed)
		    # so let's go ahead and make that.
		    #my $max_hits_fname = $top_chr_fname;
		    #$max_hits_fname =~ s/\.weaver\.out$/\-partial\.weaver\.out/;
		    #RecordMaximalHits($top_chr_fname,$max_hits_fname,$max_hits_str);
		    #
		    # ACTUALLY, let's only hold onto the sequences we'd feed into
		    # ExonWeaver -- this could be optimized, but I'm pretty sure we'd
		    # see diminishing returns...
		    my $top_xwinfname = $top_chr_fname;
		    $top_xwinfname =~ s/\.out$/\.in/;
		    SaveTopXWInputs($top_xwinfname);
		    
		}

		# As a last little step, we'll create a separate file with just the
		# maximal hits to reference after BLAT.


	    } else {

		# Alas, this string appears to be a misfit :'(
		# Let's see if BLAT+SPALN can help out at all
		print $BlatFile ">$gene\|$SeqNames[$i]\n";
		my @Seq = split(//,$Seqs[$i]);
		for (my $j=0; $j<scalar(@Seq); $j++) {
		    print $BlatFile "$Seq[$j]";
		    print $BlatFile "\n" if (($j+1) % 60 == 0);
		}
		print $BlatFile "\n" if (scalar(@Seq) % 60);
		print $BlatFile "\n";

	    }
	    #
	    #
	    #
	    # * * * * * * * * * END UNUSED CODE BLOCK * * * * * * * * * *
	    #
	}

	# HERE'S SOME MORE UNUSED CODE
    #} else {
	#
	## Straightforward copy of our sequences into the BLAT file.
	## The only minor change is that we're going to change the names
	## to be gene|seqname
	#for (my $i=0; $i<$num_seqs; $i++) {
	    #print $BlatFile ">$gene\|$SeqNames[$i]\n";
	    #my @Seq = split(//,$Seqs[$i]);
	    #for (my $j=0; $j<scalar(@Seq); $j++) {
		#print $BlatFile "$Seq[$j]";
		#print $BlatFile "\n" if (($j+1) % 60 == 0);
	    #}
	    #print $BlatFile "\n" if (scalar(@Seq) % 60 != 0);
	    #print $BlatFile "\n";
	#}

    }

    # If we have any unmapped sequences, we'll extract them from the protein
    # sequence file (into individual files) and perform a Spaln search (informed
    # by GTF coordinates)
    #
    if ($num_unmapped) {

	# Write each unmapped sequence out to its own file
	my @UnmappedSeqNames;
	my @UnmappedSeqs;
	my @ProtFnames;
	my @ProtIndices;
	for (my $i=0; $i<$num_seqs; $i++) {

	    next if ($FullMaps[$i]);
	    
	    my $prot_fname = $gene_fname;
	    $prot_fname =~ s/\.fa$/\.$SeqNames[$i]\.prot\.in/;

	    my $ProtFile = OpenOutputFile($prot_fname);
	    print $ProtFile ">$SeqNames[$i]\n";
	    my @Seq = split(//,$Seqs[$i]);
	    for (my $j=0; $j<scalar(@Seq); $j++) {
		print $ProtFile "$Seq[$j]";
		print $ProtFile "\n" if (($j+1) % 60 == 0);
	    }
	    print $ProtFile "\n";
	    close($ProtFile);

	    push(@UnmappedSeqNames,$SeqNames[$i]);
	    push(@UnmappedSeqs,$Seqs[$i]);
	    push(@ProtFnames,$prot_fname);
	    push(@ProtIndices,$i);
	    
	}

	# Do that nasty Spaln searchin'!
	my ($spaln_hits_ref,$spaln_pcts_id_ref)
	    = SpalnSearch(\@UnmappedSeqNames,\@UnmappedSeqs,\@ProtFnames,
			  \@SpalnStarts,\@SpalnEnds,\@SpalnChrs,90.0);
	my @SpalnHitStrs = @{$spaln_hits_ref};
	my @SpalnPctsID  = @{$spaln_pcts_id_ref};
	    
	# Run back through our proteins, and either (1) record full maps or
	# (2) write the sequence out to our Blat file
	for (my $i=0; $i<scalar(@ProtIndices); $i++) {

	    # In any case, let's clean out the protein file
	    RunSystemCommand("rm \"$ProtFnames[$i]\"");

	    my $prot_index = $ProtIndices[$i];
	    if ($SpalnHitStrs[$i]) {

		# I knew you had it in you!
		$SpalnHitStrs[$i] =~ s/BLAT\+Spaln/Spaln/;
		$FullMaps[$prot_index] = $SpalnHitStrs[$i];
		$num_unmapped--;

	    } else {
		
		# TO BLAT WITH YOU, FOUL SEQUENCE!
		print $BlatFile ">$gene\|$SeqNames[$prot_index]\n";
		my @Seq = split(//,$Seqs[$prot_index]);
		for (my $j=0; $j<scalar(@Seq); $j++) {
		    print $BlatFile "$Seq[$j]";
		    print $BlatFile "\n" if (($j+1) % 60 == 0);
		}
		print $BlatFile "\n" if (scalar(@Seq) % 60 != 0);
		print $BlatFile "\n";

	    }
	    
	}

    }

    # We can at least say that we got one sequence to fully map, so let's pop
    # some champagne! Or, absent champagne, open up... [gene].quilter.out!
    #
    if ($num_unmapped < $num_seqs) {
	my $outfname = $gene_fname;
	$outfname =~ s/\.fa$/\.quilter\.out/;
	my $OutFile = OpenOutputFile($outfname);
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($FullMaps[$i]) {
		print $OutFile "$FullMaps[$i]\n";
	    }
	}
	close($OutFile);
    }

    # Clear out the files -- if you need to save them, return
    for (my $i=0; $i<$num_seqs; $i++) {
	foreach my $fname (split(/\;/,$SpliceGraphs[$i])) {
	    RunSystemCommand("rm \"$fname\"") if ($fname);
	}
	foreach my $fname (split(/\;/,$XWInputNames[$i])) {
	    RunSystemCommand("rm \"$fname\"") if ($fname);
	}
    }
    if (-e $nucl_fname) { RunSystemCommand("rm \"$nucl_fname\""); }

}





############################################################
#
#  Function: ParseFastMapOutput
#
sub ParseFastMapOutput
{
    my $cmd = shift;
    my $num_seqs = shift;

    # We'll organize our hits by sequence
    my @HitsBySeq;
    for (my $i=0; $i<$num_seqs; $i++) {
	$HitsBySeq[$i] = 0;
    }

    my $fm2output = OpenSystemCommand($cmd);

    my %DuplicateCheck;
    while (my $line = <$fm2output>) {

	# Orient around the start of the next hit
	if ($line =~ /Protein Num\s+\: (\d+)/) {

	    # 1. Protein ID number
	    my $prot_id = $1;

	    # 2. Exon number (in list supplied to FastMap2)
	    $line = <$fm2output>;
	    $line =~ /Exon Num\s+\: (\d+)/;
	    my $exon_num = $1;

	    # 3. Reading frame (0-2)
	    $line = <$fm2output>;
	    $line =~ /Reading Frame\s+\: (\d)/;
	    my $reading_frame = $1;

	    # 4. Score
	    $line = <$fm2output>;
	    $line =~ /Score\s+\: (\d+[\.\d+]?)/;
	    my $score = $1;
	    # Because I'm me, I need a '.0'
	    $score = $score.'0' if ($score =~ /\./);
	    
	    # 5. Amino range (inclusive)
	    $line = <$fm2output>;
	    $line =~ /Amino Range\s+\: (\d+)\.\.(\d+)/;
	    my $amino_start = $1;
	    my $amino_end   = $2;

	    # 6. Nucleotide range (inclusive)
	    $line = <$fm2output>;
	    $line =~ /Mapped Nucl\.s\s+\: (\d+)\.\.(\d+)/;
	    my $nucl_start = $1;
	    my $nucl_end   = $2;

	    # 7. Protein characters (we already have access...)
	    $line = <$fm2output>;
	    #$line =~ /Query Protein\s+\:\s+(.+)$/;
	    #my $prot_chars = $1;
	    #$prot_chars =~ s/\s//g;

	    # 8. Nucleotide sequence (17 lowercase buffer nucls on each side)
	    $line = <$fm2output>;
	    $line =~ /Nucleotide Seq\:\s+(\S+)/;
	    my $nucl_chars = $1;

	    # 9. Translated sequence (we could generate this...)
	    $line = <$fm2output>;
	    #$line =~ /Translated Seq\:\s+(.+)$/;
	    #my $trans_chars = $1;
	    #$trans_chars =~ s/\s//g;

	    # THAT'S A HIT! WE'RE TALKIN' PLATINUM, BABY!

	    # NOTE that we're going to have the amino coordinates first, for
	    #      sorting purposes.
	    my $hit = $amino_start.'|'.$amino_end.'|'.$nucl_start.'|'.$nucl_end;

	    # If we've already seen this hit (i.e., if there were overlapping
	    # but not perfectly duplicated entries in the GTF file) we don't
	    # need to learn about it a second time
	    next if ($DuplicateCheck{$prot_id.'&'.$hit});
	    $DuplicateCheck{$prot_id.'&'.$hit} = 1;

	    # Now we'll complete the hit with the rest of the data and
	    # toss it into the hit list!
	    $hit = $hit.'|'.$exon_num.'|'.$reading_frame.'|'.$score.'|'.$nucl_chars;
	    if ($HitsBySeq[$prot_id]) {
		$HitsBySeq[$prot_id] = $HitsBySeq[$prot_id].'&'.$hit;
	    } else {
		$HitsBySeq[$prot_id] = $hit;
	    }
	    
	}
	
    }
    close($fm2output);

    # For each sequence, sort its hits by start position (wrt amino acid sequence),
    # using an algorithm known as 'idiot sort'
    for (my $i=0; $i<$num_seqs; $i++) {
	if ($HitsBySeq[$i]) {

	    my %StartAminoToSeqs;
	    my @SeqHits = split(/\&/,$HitsBySeq[$i]);
	    for (my $j=0; $j<scalar(@SeqHits); $j++) {
		$SeqHits[$j] =~ /^(\d+)\|/;
		my $start_amino = $1;
		if ($StartAminoToSeqs{$start_amino}) {
		    $StartAminoToSeqs{$start_amino} = $StartAminoToSeqs{$start_amino}.'&'.$SeqHits[$j];
		} else {
		    $StartAminoToSeqs{$start_amino} = $SeqHits[$j];
		}
	    }

	    # L A Z Y   B O I
	    my $hit_str = '';
	    foreach my $start_amino (sort { $a <=> $b } keys %StartAminoToSeqs) {
		$hit_str = $hit_str.'&'.$StartAminoToSeqs{$start_amino};
	    }
	    $hit_str =~ s/^\&//;
	    $HitsBySeq[$i] = $hit_str;

	}
    }

    return(\@HitsBySeq);
    
}




############################################################
#
#  Function: GenSpliceGraph
#
sub GenSpliceGraph
{
    my $hitstr     = shift;
    my $seqstr     = shift;
    my $seq_id     = shift;
    my $chr        = shift;
    my $revcomp    = shift;
    my $gene_fname = shift;

    # We'll include the chromosome name in the output file name,
    # so we'll want to make sure strand direction is reflected
    $chr = $chr.'-revcomp' if ($revcomp);

    # Good to know how long the protein sequence is
    my $seq_len = length($seqstr);
    
    # Break the hit string into its individual hits
    my @Hits     = split(/\&/,$hitstr);
    my $num_hits = scalar(@Hits);

    # Now we can prep a lil' file for ExonWeaver
    my $weaver_in = $gene_fname;
    $weaver_in =~ s/\.fa$/\-$seq_id\.$chr\.weaver\.in/;
    my $WeaverFile = OpenOutputFile($weaver_in);

    # NOTE: I'm going to generally assume that between all my tools
    #  we have a standardized non-coding buffer of 7 characters
    # (rationale being that, when stitching two hits together, you
    #  can extend each hit out 1 additional character, and then find
    #  a split codon, without losing access to intron dinucleotides).
    
    # First, we'll provide some metadata, just so nobody's caught off-guard
    print $WeaverFile "Num Hits : $num_hits\n";
    print $WeaverFile "Seq Len  : $seq_len\n";
    print $WeaverFile "Sequence : $seqstr\n";

    # For each hit, we'll generate a (naive) representation of the
    # splice signal strength at each position along its associated
    # nucleotide sequence.  The score indicates the strength of
    # the position as the first or last nucleotide of the INTRON.
    foreach my $hit (@Hits) {
	print $WeaverFile "\n";
	PrintHitToWeaverInf($hit,$WeaverFile);
    }


    # Close up the file so we can put it to work!
    close($WeaverFile);

    # Also, in case you're like me, and you've hit this point and started wondering,
    # "Hey, why aren't we just piping FastMap2 output into ExonWeaver?" the main reason
    # is that we needed to divide FastMap2 output according to the protein sequences
    # and then sort them so that we have a guarantee of ascending start aminos.

    # We'll return the name of this file
    my $weaver_out = $weaver_in;
    $weaver_out =~ s/\.in$/\.out/;

    # PUT IT TO WORK!!!
    my $weaver_cmd = $srcdir."ExonWeaver --report-singles \"$weaver_in\" > \"$weaver_out\"";
    RunSystemCommand($weaver_cmd);

    # If there's any content to the output file, return its name -- otherwise,
    # delete the empty file and return 0.  Note that it'd be somewhat bizarre to
    # get no output, given that we've enabled --report-singles...
    if (-s $weaver_out) {
	return ($weaver_in,$weaver_out);
    } else {
	RunSystemCommand("rm \"$weaver_in\"");
	RunSystemCommand("rm \"$weaver_out\"") if (-e $weaver_out);
	return (0,0);
    }
    
}





############################################################
#
#
#
sub PrintHitToWeaverInf
{
    my $hit = shift;
    my $WeaverFile = shift;
    
    # Extract the nucleotide sequence for this hit
    my @HitData = split(/\|/,$hit);
    my $nucl_str = uc($HitData[7]);
    
    # Compute the strengths of the splice sites
    my ($three_prime_str,$five_prime_str) = GetSpliceStrengths($nucl_str);

    # One last thing -- make sure we have clear float formatting for the score
    if ($HitData[6] !~ /\.\d+$/) {
	if ($HitData[6] !~ /\./) { $HitData[6] = $HitData[6].'.'; }
	$HitData[6] = $HitData[6].'0';
    }
    
    # Now that we have that additional whiff of data in the mix,
    # let's get writing!
    print $WeaverFile "Hit Score   : $HitData[6]\n";
    print $WeaverFile "Amino Range : $HitData[0]\.\.$HitData[1]\n";
    print $WeaverFile "Nucl Range  : $HitData[2]\.\.$HitData[3]\n";
    print $WeaverFile "Nucleotides : $nucl_str\n";
    print $WeaverFile "3' SS Str.  : $three_prime_str\n";
    print $WeaverFile "5' SS Str.  : $five_prime_str\n";
	
}





############################################################
#
#  Function: GetSpliceStrengths
#
#  NOTE: This is currently just checking for the canonical dinucleotides.
#        Eventually we'll want it to be based on something more nuanced.
#
sub GetSpliceStrengths
{
    my $nucl_str = shift;
    my @Nucls = split(//,$nucl_str);

    # Because we can't evaluate dinucleotides at the first or last positions,
    # those default to 0 -- I'm pretending these are probabilities.
    my $three_prime_str = '0.00';
    my $five_prime_str  = '0.00';

    # Could you be more naive?
    my $canon_score    = 0.95;
    my $noncanon_score = 0.05;
    for (my $i=1; $i<scalar(@Nucls)-1; $i++) {
	if ($Nucls[$i-1].$Nucls[$i] eq 'AG') {
	    $three_prime_str = $three_prime_str.' '.$canon_score;
	    $five_prime_str  = $five_prime_str.' '.$noncanon_score;
	} elsif ($Nucls[$i-1].$Nucls[$i] eq 'GT') {
	    $three_prime_str = $three_prime_str.' '.$canon_score;
	    $five_prime_str  = $five_prime_str.' '.$noncanon_score;
	} else {
	    $three_prime_str = $three_prime_str.' '.$noncanon_score;
	    $five_prime_str  = $five_prime_str.' '.$noncanon_score;
	}
    }
    
    $three_prime_str = $three_prime_str.' 0.00';
    $five_prime_str  = $five_prime_str.' 0.00';

    return ($three_prime_str,$five_prime_str);
    
}





############################################################
#
#  Function: CheckForFullMaps
#
sub CheckForFullMaps
{
    my $splicegraphs = shift;
    my $seq_len = shift;

    # If we find any full mappings, we'll return those filenames
    my $full_maps = 0;
    foreach my $xw_outfname (split(/\;/,$splicegraphs)) {

	my $xwf = OpenInputFile($xw_outfname);
	while (my $line = <$xwf>) {
	    $line =~ s/\n|\r//g;
	    if ($line =~ /Amino Acid Range : 1\.\.$seq_len/) {
		# GOT ONE!
		if ($full_maps) { $full_maps = $full_maps.';'.$xw_outfname; }
		else            { $full_maps = $xw_outfname;                }
		last;
	    }
	}
	close($xwf);
	
    }

    # Couldn't have been easier!
    return $full_maps;

}





############################################################
#
#  Function: ExtractCanonXWMap
#
#  NOTE: This function assumes that we have a full mapping to a single chromosome.
#        Any full mappings using chimeric inputs to ExonWeaver will need to use
#        ExtractChimericXWMap.
#
sub ExtractCanonXWMap
{
    my $fname  = shift;
    my $seqlen = shift;
    my $seq_id = shift;
    my $chr    = shift;

    # Open the file and scan until we we're in the full mapping zone (FMZ)
    my $inf = OpenInputFile($fname);
    my $line = <$inf>;
    while ($line = <$inf>) {
	$line =~ s/\n|\r//g;
	last if ($line =~ /Amino Acid Range \: 1\.\.$seqlen$/);
    }

    # I'm not going to tolerate abuse of this function!
    if (eof($inf)) { close($inf); die "\n  ERROR:  '$fname' didn't fully map!\n\n"; }
    
    # Eat the nucleotide range and subsequent empty line
    $line = <$inf>;
    $line = <$inf>;

    # All lines until an empty line will tell us the coordinates of particular
    # specific ranges of exons.
    my @AminoRanges;
    my @NuclRanges;
    my $num_exons = 0;
    while ($line = <$inf>) {
	
	$line =~ s/\n|\r//g;
	last if (!$line);

	$line =~ /Exon \d+\: Aminos (\d+\.\.\d+)\, Nucls (\d+\.\.\d+)/;
	$AminoRanges[$num_exons] = $1;
	$NuclRanges[$num_exons]  = $2;
	$num_exons++;

    }

    # Now we're into the actual mapping, which we'll flatten into 3 linear arrays
    # (flatten in the sense that we're removing line formatting)
    my $trans_str = '';
    my $nucl_str  = '';
    my $prot_str  = '';
    $line = <$inf>;
    $line =~ s/\n|\r//g;
    while ($line) {
	
	# Each line begins with two spaces
	$line =~ s/^  //;
	$trans_str = $trans_str.$line;

	$line = <$inf>;
	$line =~ s/\n|\r//g;
	$line =~ s/^  //;
	$nucl_str = $nucl_str.$line;

	$line = <$inf>;
	$line =~ s/\n|\r//g;
	$line =~ s/^  //;
	$prot_str = $prot_str.$line;

	# We'll either grab 2 blank lines (end of mapping) or else cycle back
	# to a translated line
	$line = <$inf>;
	$line = <$inf>;
	$line =~ s/\n|\r//g;
	
    }

    # All finished with the fileage -- now it's time for the parseage!
    close($inf);

    # Break those baddies into arrays
    my @Trans = split(//,$trans_str);
    my @Nucls = split(//,$nucl_str);
    my @Prot  = split(//,$prot_str);

    # Might need to make a teensy change, for formatting
    my $strand = 1;
    if ($chr =~ /\-revcomp$/) {
	$chr =~ s/\-revcomp$/\[revcomp\]/;
	$strand = -1;
    }
    
    # Time to start constructing a hit string!
    my $hitstr = "Sequence ID: $seq_id\n";
    $hitstr = $hitstr."Map Method : FastMap2\n";
    $hitstr = $hitstr."Chromosome : $chr\n";
    $hitstr = $hitstr."Num Exons  : $num_exons\n";

    # Walk along, exon-by-exon, building up the official mapping!
    my $exon_num = 0;
    my $scan = 0;
    my $nucl_pos;
    while ($exon_num < $num_exons) {

	$hitstr = $hitstr."* Aminos $AminoRanges[$exon_num], $chr:$NuclRanges[$exon_num]\n";

	# Where does this exon start?
	$NuclRanges[$exon_num] =~ /^(\d+)\./;
	$nucl_pos = $1;

	# If we're at an internal exon, there may be flanking splice signal characters
	while ($Nucls[$scan] eq lc($Nucls[$scan])) {
	    $scan++;
	}

	# Now run until the end of the exon (either the end of the sequence or the
	# splice signal characters
	while ($scan < scalar(@Nucls) && $Nucls[$scan] eq uc($Nucls[$scan])) {
	    if ($Trans[$scan] =~ /[A-Z]/) {
		$hitstr = $hitstr.$nucl_pos.',';
	    }
	    $nucl_pos += $strand;
	    $scan++;
	}
	$hitstr =~ s/\,$//;
	$hitstr = $hitstr."\n";

	# That's that for that exon (that that that that)
	$exon_num++;

    }

    return $hitstr;
    
}






############################################################
#
#  Function: AttemptSpalnFill
#
sub AttemptSpalnFill
{
    my $xw_outfname       = shift;
    my $max_hits_list_str = shift;
    my $gap_starts_ref    = shift;
    my $gap_ends_ref      = shift;
    my $seq_str           = shift;
    my $seq_id            = shift;
    my $genome            = shift;
    my $gene_fname        = shift;

    my @GapStarts = @{$gap_starts_ref};
    my @GapEnds   = @{$gap_ends_ref};
    my $num_gaps  = scalar(@GapStarts);

    my $seq_len = length($seq_str);
    my @Seq = split(//,$seq_str);

    # Figure out what chromosome we're dealing with
    $xw_outfname =~ /\/([^\/]+)$/;
    my $hitfname = $1;
    $hitfname =~ /\.([^\.]+)\.weaver\.out/;
    my $chr = $1;
    my $revcomp = 0;
    if ($chr =~ /\-revcomp$/) {
	$chr =~ s/\-revcomp$//;
	$revcomp = 1;
    }

    # We'll kick things off by opening up the ExonWeaver output file and grabbing the
    # nucleotide ranges corresponding to each of our favorite hits.
    my $xwoutf = OpenInputFile($xw_outfname);
    my @MaxHits = split(/\,/,$max_hits_list_str); # These are [0..num_hits-1]
    my @MaxHitNuclRanges;
    my $hit_num = 0;
    my $exon_num = 0;
    while (my $line = <$xwoutf>) {

	if ($line =~ /\-\-\-\-\- Exons/) {
	    if ($exon_num == $MaxHits[$hit_num]) {

		# Scan to where the nucleotide range is hiding
		$line = <$xwoutf>; # blank line
		$line = <$xwoutf>; # exon id list
		$line = <$xwoutf>; # mapping score
		$line = <$xwoutf>; # amino range
		$line = <$xwoutf>; # Nucleotide Range!
		$line =~ s/\n|\r//g;
		$line =~ /Nucleotide Range \: (\d+\.\.\d+)$/;

		# Record that this hit had the given nucleotide range
		push(@MaxHitNuclRanges,$1);

		$hit_num++;
		last if ($hit_num == scalar(@MaxHits));

	    }
	    $exon_num++;
	}
	
    }
    close($xwoutf);

    # Now we can run through our list of gaps and determine where we'd
    # want to look in the genome to try to fill them in
    my @NuclSearchRanges;

    # First, if we couldn't map the first aminos of the sequence, we'll
    # pull in a fairly large search range.
    if ($GapStarts[0] == 1) {
	$MaxHitNuclRanges[0] =~ /^(\d+)\.\./;
	my $range_end   = $1;
	my $range_start;
	if ($revcomp) { $range_start = Min($range_end+1000000,$ChrSizes{$chr}); }
	else          { $range_start = Max($range_end-1000000,1);               }
	push(@NuclSearchRanges,$range_start.'..'.$range_end);
    }

    # We'll handle all of the guaranteed internal gaps
    for (my $i=0; $i<scalar(@MaxHits)-1; $i++) {
	$MaxHitNuclRanges[$i] =~ /\.\.(\d+)$/;
	my $range_start = $1;
	$MaxHitNuclRanges[$i+1] =~ /^(\d+)\.\./;
	my $range_end = $1;
	push(@NuclSearchRanges,$range_start.'..'.$range_end);
    }

    # If the final gap extends to the end of the protein sequence, we'll need to
    # handle it specially -- otherwise, we've got our ranges!
    if ($GapEnds[$num_gaps-1]==$seq_len) {
	$MaxHitNuclRanges[scalar(@MaxHits)-1] =~ /\.\.(\d+)$/;
	my $range_start = $1;
	my $range_end;
	if ($revcomp) { $range_end = Max($range_start-1000000,0);               }
	else          { $range_end = Min($range_start+1000000,$ChrSizes{$chr}); }
	push(@NuclSearchRanges,$range_start.'..'.$range_end);
    }

    
    # Heckin' dang! Let's do some SPALNery!
    #
    # Our goal is just to find the specific ranges of unindexed exons in the gap
    # areas that match up with the specific unmapped amino sequences.
    #
    # NOTE that we're fine with not fully filling in every gap -- any gap that we
    # can fill in is a success!

    
    # We're going to need to give Spaln some files, duh!
    my $temp_fname = $gene_fname;
    $temp_fname =~ s/\.fa$/\.spaln\.prot\.fa/;
    my $nucl_fname = $temp_fname;
    $nucl_fname =~ s/\.prot\.fa$/\.nucl\.fa/;
    my $spaln_fname = $nucl_fname;
    $spaln_fname =~ s/\.nucl\.fa$/\.spaln\.out/;

    my @GapHits;
    for (my $i=0; $i<$num_gaps; $i++) {

	open(my $tempf,'>',$temp_fname);
	print $tempf ">$i\n";
	for (my $j=$GapStarts[$i]-1; $j<$GapEnds[$i]; $j++) {
	    print $tempf "$Seq[$j]";
	}
	print $tempf "\n\n";
	close($tempf);
	
	my $nucl_range = $NuclSearchRanges[$i];
	my $sfetch_cmd = $sfetch." -range $nucl_range \"$genome\" \"$chr\" > \"$nucl_fname\"";
	RunSystemCommand($sfetch_cmd);

	my $spaln_cmd = $spaln."\"$nucl_fname\" \"$temp_fname\" 1>\"$spaln_fname\" 2>/dev/null";

	# Don't bail if Spaln does -- just move on to the next case
	next if (system($spaln_cmd));
	my $SpalnOut = OpenInputFile($spaln_fname);

	# Note that while the intuition is that a gap is a single wily
	# exon, it's possible for there to be multiple exons (enough to
	# require multiple 'join' lines), so we use an array to capture
	# each presumptive range
	my $coordlist_str = '';
	while (!eof($SpalnOut)) {

	    my $line = <$SpalnOut>;
	    while ($line && $line =~ /\;C (\S+)/) {

		my $line_coords = $1;

		# We might need to do some cleanup of this line, if it's
		# the first or last of the list.
		$line_coords =~ s/complement|join|\(|\)|\s//g;
		$coordlist_str = $coordlist_str.$line_coords;

		$line = <$SpalnOut>;

	    }

	    last if ($coordlist_str);
	    
	}
	close($SpalnOut);

	# If we didn't get any hints from Spaln, move onto the next range
	next if (!$coordlist_str);

	# OOOOhhOhoHohoh! Looks like we might have a coding region!
	# Time to see if FastMap can help us out with the specifics...
	my $mapcmd = $srcdir."FastMap2 \"$temp_fname\" 1 \"$nucl_fname\"";
	foreach my $coord_pair (split(/\,/,$coordlist_str)) {
	    $coord_pair =~ /^(\d+)\.\.(\d+)$/;
	    $mapcmd = $mapcmd.' '.$1.' '.$2;
	}
	my $hitref  = ParseFastMapOutput($mapcmd,1);
	my @SeqHits = @{$hitref};
	next if (!$SeqHits[0]);

	# As above, so now! We have hits stored in an &-separated list, with
	# individual hit data fields separated by |s, and ordered as such:
	#
	#   0. amino start position
	#   1. amino end   position
	#   2. nucl  start position
	#   3. nucl  end   position
	#   4. exon number
	#   5. reading frame (0-2)
	#   6. score (BLOSUM 62, half-bit)
	#   7. nucleotide sequence (w/ 17 lowercase buffer nucls on each side)
	#
	
	# We'll need to adjust the coordinates for both the aminos and nucleotides.
	$nucl_range =~ /(\d+)\.\./;
	my $gap_nucl_start = $1;
	    
	foreach my $hit (split(/\&/,$SeqHits[0])) {

	    my @HitData = split(/\|/,$hit);

	    $HitData[0] += $GapStarts[$i] - 1;
	    $HitData[1] += $GapStarts[$i] - 1;

	    if ($revcomp) {
		$HitData[2] = $gap_nucl_start - $HitData[2] + 1;
		$HitData[3] = $gap_nucl_start - $HitData[3] + 1;
	    } else {
		$HitData[2] += $gap_nucl_start - 1;
		$HitData[3] += $gap_nucl_start - 1;
	    }

	    # Record the updated hit!
	    $hit = $HitData[0].'|'.$HitData[1].'|'.$HitData[2].'|'.$HitData[3];
	    $hit = $hit.'|'.$HitData[4].'|'.$HitData[5].'|'.$HitData[6].'|'.$HitData[7];
	    push(@GapHits,$hit);
	    
	}
	
    }

    # A bit of cleanup and futility checking
    RunSystemCommand("rm \"$temp_fname\"")  if (-e $temp_fname);
    RunSystemCommand("rm \"$nucl_fname\"")  if (-e $nucl_fname);
    RunSystemCommand("rm \"$spaln_fname\"") if (-e $spaln_fname);
    return 0 if (!scalar(@GapHits));
    

    # Uber chill, bro!  Looks like we've got (hopefully) a list of
    # nucleotide coordinates corresponding to the presumptive exons
    # that we were able to get out of running SPALN.
    #
    # Next up, we'll pull in the original FastMap2 hits and add them
    # to GapHits.  This will give us a full input to push off to
    # ExonWeaver, which will hopefully give us a full-sequence mapping...


    # Open up the original ExonWeaver input file
    my $xw_infname = $xw_outfname;
    $xw_infname =~ s/out$/in/;
    my $xwinf = OpenInputFile($xw_infname);

    while (my $line = <$xwinf>) {

	# Is this the start of a new exon?
	if ($line =~ /Hit Score   \: (\S+)/) {

	    my $hit_score = $1 + 0.0;

	    $line = <$xwinf>; # Amino range
	    $line =~ /Amino Range \: (\d+)\.\.(\d+)/;
	    my $hit_amino_start = $1;
	    my $hit_amino_end   = $2;
	    
	    $line = <$xwinf>; # Nucl  range
	    $line =~ /Nucl Range  \: (\d+)\.\.(\d+)/;
	    my $hit_nucl_start = $1;
	    my $hit_nucl_end   = $2;
	    
	    $line = <$xwinf>; # Nucleotide Seq.
	    $line =~ /Nucleotides \: (\S+)/;
	    my $hit_nucl_str = $1;

	    # We don't store the 3' / 5' SS strength vals -- these will get
	    # recomputed momentarily.  Also, note that we use '-' for the exon
	    # and reading frame info.
	    my $hit = $hit_amino_start.'|'.$hit_amino_end.'|'.$hit_nucl_start;
	    $hit    = $hit.'|'.$hit_nucl_end.'|-|-|'.$hit_score.'|'.$hit_nucl_str;
	    push(@GapHits,$hit);
	    
	}
	
    }

    close($xwinf);

    # Sort those stinky hits!
    my %StartAminoToHits;
    foreach my $hit (@GapHits) {
	$hit =~ /^(\d+)\|/;
	my $start_amino = $1;
	if ($StartAminoToHits{$start_amino}) {
	    $StartAminoToHits{$start_amino} = $StartAminoToHits{$start_amino}.'&'.$hit;
	} else {
	    $StartAminoToHits{$start_amino} = $hit;
	}
    }
    my $final_exon_str = '';
    foreach my $hit_amino_start (sort { $a <=> $b } keys %StartAminoToHits) {
	foreach my $hit (split(/\&/,$StartAminoToHits{$hit_amino_start})) {
	    $final_exon_str = $final_exon_str.'&'.$hit;
	}
    }
    $final_exon_str =~ s/^\&//;

    # At this point, we can assume that we'll at least recover the existing output
    # data (if not grow it), so we clear the prior files (mainly so we can reuse
    # the file names...)
    RunSystemCommand("rm \"$xw_infname\"");
    RunSystemCommand("rm \"$xw_outfname\"");

    # Push off to ExonWeaver!
    ($xw_infname,$xw_outfname) =
	GenSpliceGraph($final_exon_str,$seq_str,$seq_id,$chr,$revcomp,$gene_fname);

    # Did we get a full mapping?
    my $hit_str = CheckForFullMaps($xw_outfname,$seq_len);

    # If we did get a full mapping, give credit to SPALN for filling in the gaps
    if ($hit_str) {
	$hit_str =~ s/Map Method \: FastMap2/Map Method \: FastMap2\+Spaln/;
    }

    # DONE
    return $hit_str;

}






############################################################
#
#  Function: AttemptChimericXWMap
#
sub AttemptChimericXWMap
{
    my $orig_xw_inf_str   = shift;
    my $orig_xw_outf_str  = shift;
    my $canon_chr_fname   = shift;
    my $max_exon_list_str = shift;
    my $amino_starts_ref  = shift;
    my $amino_ends_ref    = shift;
    my $seq_str           = shift;
    my $seq_id            = shift;
    my $gene_fname        = shift;

    my @XWOutNames     = split(/\;/,$orig_xw_outf_str);
    my @XWInNames      = split(/\;/,$orig_xw_inf_str);
    my @AminoGapStarts = @{$amino_starts_ref};
    my @AminoGapEnds   = @{$amino_ends_ref};

    # I'm going to build a list of all the exons we want to use as inputs to
    # ExonWeaver.  We'll need a method for sorting these by their starting aminos, too.
    my @FullExonList;
    my %ExonInfoToListIndex;
    my @ChrsByExon;

    # Next up, we'll run through each of the non-canonical chromosome's input
    # files and find all of the putative exons that fall into the gaps between
    # canonical hits
    my $num_exons = 0;
    for (my $i=0; $i<scalar(@XWInNames); $i++) {

	# Which chromosome is this?
	$XWInNames[$i] =~ /\.([^\.]+)\.weaver\.in$/;
	my $chr = $1;

	# We can treat the canonical chromosome somewhat specially, since it has
	# a list pre-generated
	my @ExonList;
	if ($XWOutNames[$i] eq $canon_chr_fname) {

	    my $exon_list_ref
		= ExtractXWInputExons($canon_chr_fname,$max_exon_list_str);

	    @ExonList = @{$exon_list_ref};
	    for (my $j=0; $j<scalar(@ExonList); $j++) {

		my $exon = $ExonList[$j];
		push(@FullExonList,$exon);
		push(@ChrsByExon,$chr);
		
		$exon =~ /Amino Range \: (\d+)\.\.(\d+)/;
		my $start_amino = $1;
		my $end_amino = $2;
		my $key = $start_amino.'|'.$end_amino.'|'.$i;
		
		# We need to have the lowest entry be '1' por supuesto
		$num_exons++;
		$ExonInfoToListIndex{$key} = $num_exons;

	    }
	

	} else {

	    # Noncanonical chromosome!  We'll need to search the input file for
	    # exons that fit in the canonical gaps (with a bit of wiggle room for
	    # overlap).
	    my $xwinf = OpenInputFile($XWInNames[$i]);
	    my $gap_num = 0;
	    while (my $line = <$xwinf>) {

		# Exon is kicking off!
		if ($line =~ /^Hit Score/) {

		    # Begin the exon tracking
		    my $exon = $line;

		    $line = <$xwinf>;
		    $line =~ /Amino Range \: (\d+)\.\.(\d+)/;
		    my $start_amino = $1;
		    my $end_amino   = $2;
		    $exon = $exon.$line;

		    # Scan to the next gap whose start position is no more than
		    # 2 aminos higher than this hit's start;
		    while ($gap_num < scalar(@AminoGapStarts) && $AminoGapStarts[$gap_num] > $start_amino + 2) {
			$gap_num++;
		    }

		    # If we've exceeded our array, we're done with this file
		    last if ($gap_num == scalar(@AminoGapStarts));

		    # Is this hit in the gap?
		    if ($start_amino >= $AminoGapStarts[$gap_num]-2 &&
			$end_amino   <= $AminoGapEnds[$gap_num]+2) {

			# Excellent!  Build it up and add it to the list
			$line = <$xwinf>; # Nucl. Range
			$exon = $exon.$line;
			$line = <$xwinf>; # Nucl. String
			$exon = $exon.$line;
			$line = <$xwinf>; # 3' SS Strengths
			$exon = $exon.$line;
			$line = <$xwinf>; # 5' SS Strengths
			$exon = $exon.$line."\n";

			push(@FullExonList,$exon);
			push(@ChrsByExon,$chr);

			my $key = $start_amino.'|'.$end_amino.'|'.$i;

			$num_exons++;
			$ExonInfoToListIndex{$key} = $num_exons;
			
		    }

		}

	    }

	    close($xwinf);

	}
	
    }

    # In case we *only* hit to a noncanonical chromosome, we might not have pulled
    # any exons in for this check, so we'll want to jump off.
    return 0 if ($num_exons == 0);

    # NOTE that we don't check to see whether there's really any chance of getting
    # a chimeric hit, since we'll want to hang onto the original hits this sequence
    # had to the canonical chromosome anyways.

    # Now, we should be able to sort our keys and put of all our
    # exons into a file in order of start amino
    my $chimera_in = $gene_fname;
    $chimera_in =~ s/\.fa$/\-$seq_id\.chimera\.weaver\.in/;
    my $chimera_out = $chimera_in;
    $chimera_out =~ s/\.in$/\.out/;
    my $ChimeraFile = OpenOutputFile($chimera_in);

    # Need to print out some metadata to the Chimeric file
    my $seq_len = length($seq_str);
    print $ChimeraFile "Num Hits : $num_exons\n";
    print $ChimeraFile "Seq Len  : $seq_len\n";
    print $ChimeraFile "Sequence : $seq_str\n\n";

    # Now we'll print out each exon, going in order of start amino
    foreach my $key (sort keys %ExonInfoToListIndex) {
	my $exon = $FullExonList[$ExonInfoToListIndex{$key}-1];
	print $ChimeraFile "$exon";
    }

    # Close up the file and give it to ExonWeaver!  Be sure to let them know
    # that you don't want to worry about chromosome consistency...
    #
    # NOTE: We aren't interested in single-exon hits this time -- if we don't
    #       get any decent splicing on the canonical exon, we'll pass this sequence
    #       entirely over to BLAT+SPALN
    #
    close($ChimeraFile);
    my $weaver_cmd = $srcdir."ExonWeaver --allow-inconsistency \"$chimera_in\" > \"$chimera_out\"";
    RunSystemCommand($weaver_cmd);

    # We don't need to hold onto this file anymore
    RunSystemCommand("rm \"$chimera_in\"");
    
    # If we don't have any content in the output file, bail on the chimeric search
    if (!(-s $chimera_out)) {
	RunSystemCommand("rm \"$chimera_out\"") if (-e $chimera_out);
	return 0;
    }

    # If we have a full mapping, we'll need to extract it!
    my $full_map = CheckForFullMaps($chimera_out,$seq_len);
    if ($full_map) {
	$full_map = ExtractChimericXWMap($chimera_out,\@ChrsByExon,$seq_len,$seq_id);
    }

    # No more need for this file -- burn it and punt the full map (if we have one)
    RunSystemCommand("rm \"$chimera_out\"") if (-e $chimera_out);
    return $full_map;

}




############################################################
#
#  Function: ExtractXWInputExons
#
sub ExtractXWInputExons
{
    my $fname = shift;
    my $exon_list_str = shift;

    my @ExonList = split(/\,/,$exon_list_str);
    my $inf = OpenInputFile($fname);

    my @ExonData;
    my $exon_index = 0;
    my $exon_count = 0;
    while (my $line = <$inf>) {

	# Are we at the start of another exon?
	if ($line =~ /^Hit Score/) {

	    # Because the indexing method is [1..num_seqs], we increment early
	    $exon_count++;

	    # Is this the next canon hit exon?
	    if ($exon_count == $ExonList[$exon_index]) {
	    
		# You bet it is! Take note of its start position in the amino sequence
		# and store it.
		my $canon_exon_str = $line;
		$line = <$inf>; # Amino Range
		$canon_exon_str = $canon_exon_str.$line;
		$line = <$inf>; # Nucl. Range
		$canon_exon_str = $canon_exon_str.$line;
		$line = <$inf>; # Nucl. String
		$canon_exon_str = $canon_exon_str.$line;
		$line = <$inf>; # 3' SS Strengths
		$canon_exon_str = $canon_exon_str.$line;
		$line = <$inf>; # 5' SS Strengths
		$canon_exon_str = $canon_exon_str.$line."\n";
		
		push(@ExonData,$canon_exon_str);
		$exon_index++;
		last if ($exon_index == scalar(@ExonList));
		
	    }

	}

    }

    close($inf);

    return(\@ExonData);
    
}




############################################################
#
#  Function: ExtractChimericXWMap
#
#  NOTE: To give the opposite note from ExtractCanonXWMap, this function doesn't
#        assume that all hits are to the same strand / chromosome.
#
sub ExtractChimericXWMap
{
    my $xw_outfname = shift;
    my $chr_list_ref = shift;
    my $seqlen = shift;
    my $seq_id = shift;

    my @ChrsByExon = @{$chr_list_ref};
    
    # Open the file and scan until we we're in the full mapping zone (FMZ)
    my $inf = OpenInputFile($xw_outfname);
    my $line = <$inf>;
    while ($line = <$inf>) {
	$line =~ s/\n|\r//g;
	last if ($line =~ /Amino Acid Range \: 1\.\.$seqlen$/);
    }

    # I'm not going to tolerate abuse of this function!
    if (eof($inf)) { close($inf); die "\n  ERROR:  '$xw_outfname' didn't fully map!\n\n"; }
    
    # Eat the nucleotide range and subsequent empty line
    $line = <$inf>;
    $line = <$inf>;

    # All lines until an empty line will tell us the coordinates of particular
    # specific ranges of exons.
    # We'll need to know the indices of the exons so we can trace them back to their
    # chromosomes.
    my @ExonIndices;
    my @AminoRanges;
    my @NuclRanges;
    my $num_exons = 0;
    while ($line = <$inf>) {
	
	$line =~ s/\n|\r//g;
	last if (!$line);

	$line =~ /Exon (\d+)\: Aminos (\d+\.\.\d+)\, Nucls (\d+\.\.\d+)/;
	$ExonIndices[$num_exons] = $1-1;
	$AminoRanges[$num_exons] = $2;
	$NuclRanges[$num_exons]  = $3;
	$num_exons++;

    }

    # Now we're into the actual mapping, which we'll flatten into 3 linear arrays
    # (flatten in the sense that we're removing line formatting)
    my $trans_str = '';
    my $nucl_str  = '';
    my $prot_str  = '';
    $line = <$inf>;
    $line =~ s/\n|\r//g;
    while ($line) {
	
	# Each line begins with two spaces
	$line =~ s/^  //;
	$trans_str = $trans_str.$line;

	$line = <$inf>;
	$line =~ s/\n|\r//g;
	$line =~ s/^  //;
	$nucl_str = $nucl_str.$line;

	$line = <$inf>;
	$line =~ s/\n|\r//g;
	$line =~ s/^  //;
	$prot_str = $prot_str.$line;

	# We'll either grab 2 blank lines (end of mapping) or else cycle back
	# to a translated line
	$line = <$inf>;
	$line = <$inf>;
	$line =~ s/\n|\r//g;
	
    }

    # All finished with the fileage -- now it's time for the parseage!
    close($inf);

    # Break those baddies into arrays
    my @Trans = split(//,$trans_str);
    my @Nucls = split(//,$nucl_str);
    my @Prot  = split(//,$prot_str);

    # Time to start constructing a hit string!
    my $hitstr = "Sequence ID: $seq_id\n";
    $hitstr = $hitstr."Map Method : FastMap2\n";
    $hitstr = $hitstr."Chromosome : Chimeric\n";
    $hitstr = $hitstr."Num Exons  : $num_exons\n";

    # Because we use this function when filling gaps with BLAT, it's possible
    # that we actually don't have a chimeric hit here, so we'll track whether
    # multiple chromosomes are actually used...
    my %ChrTracker;

    # Walk along, exon-by-exon, building up the official mapping!
    my $exon_num = 0;
    my $scan = 0;
    my $nucl_pos;
    while ($exon_num < $num_exons) {

	# Need to determine the chromosome for each exon
	my $chr = $ChrsByExon[$ExonIndices[$exon_num]];

	# Might need to make a teensy change, for formatting
	my $strand = 1;
	if ($chr =~ /\-revcomp$/) {
	    $chr =~ s/\-revcomp$/\[revcomp\]/;
	    $strand = -1;
	}
	$ChrTracker{$chr} = 1;
    
	$hitstr = $hitstr."* Aminos $AminoRanges[$exon_num], $chr:$NuclRanges[$exon_num]\n";

	# Where does this exon start?
	$NuclRanges[$exon_num] =~ /^(\d+)\./;
	$nucl_pos = $1;

	# If we're at an internal exon, there may be flanking splice signal characters
	while ($Nucls[$scan] eq lc($Nucls[$scan])) {
	    $scan++;
	}

	# Now run until the end of the exon (either the end of the sequence or the
	# splice signal characters
	while ($scan < scalar(@Nucls) && $Nucls[$scan] eq uc($Nucls[$scan])) {
	    if ($Trans[$scan] =~ /[A-Z]/) {
		$hitstr = $hitstr.$nucl_pos.',';
	    }
	    $nucl_pos += $strand;
	    $scan++;
	}
	$hitstr =~ s/\,$//;
	$hitstr = $hitstr."\n";

	# That's that for that exon (that)
	$exon_num++;

    }

    # If we only had one chromosome, this ain't no chimera!
    my @Chrs = keys %ChrTracker;
    if (scalar(@Chrs) == 1) {
	my $chr = $Chrs[0];
	$hitstr =~ s/Chromosome \: Chimeric/Chromosome \: $chr/;
    }

    print "$hitstr\n";
    
}






############################################################
#
#  Function: FindMaximalXWHitSet
#
sub FindMaximalXWHitSet
{
    my $outf_list_str    = shift;
    my $target_chr_fname = shift;

    my @OutFileNames = split(/\;/,$outf_list_str);

    # First, find which files correspond to our target chromosome
    my $target_index = 0;
    while ($OutFileNames[$target_index] ne $target_chr_fname) {
	$target_index++;
    }

    # Open up the target chromosome's XW output file and pull in
    # a list of the hits our sequence had to it, with exon indices.
    # Keep in mind that these exon indices are [1..num_exons],
    # with respect to the order of exons in the input file.
    my @TargetHitStarts;
    my @TargetHitEnds;
    my @TargetHitExons;
    my $xwoutf   = OpenInputFile($OutFileNames[$target_index]);
    my $num_hits = 0;
    while (my $line = <$xwoutf>) {

	next if ($line !~ /\-\-\-\-\- Exons/);

	$line = <$xwoutf>; # blank line
	$line = <$xwoutf>; # Exon list!
	$line =~ /Spliced Exon IDs \: (\S+)/;
	$TargetHitExons[$num_hits] = $1;

	$line = <$xwoutf>; # score line
	$line = <$xwoutf>; # Amino range!
	$line =~ /Amino Acid Range \: (\d+)\.\.(\d+)/;
	$TargetHitStarts[$num_hits] = $1;
	$TargetHitEnds[$num_hits]   = $2;

	$num_hits++;
	
    }
    close($xwoutf);

    # Where are there holes in our mapping?
    # NOTE: Because the output method in XW runs linearly through the list of exons
    #       to see which start hits, they should be pre-sorted according to amino
    #       start position, so we can take advantage of that when determining the
    #       set of non-overlapping hits that maximally cover the protein
    my @MaximalHitSet;
    push(@MaximalHitSet,0); # Prime it!
    my $num_top_hits = 1;
    my $top_recent_hit = 0;
    for (my $i=1; $i<$num_hits; $i++) {

	# If this overlaps with our most recent top hit, see which covers more of the
	# sequence.
	if ($TargetHitStarts[$i] <= $TargetHitEnds[$top_recent_hit]) {

	    my $recent_coverage = 1 + $TargetHitEnds[$top_recent_hit] - $TargetHitStarts[$top_recent_hit];
	    my $new_coverage = 1 + $TargetHitEnds[$i] - $TargetHitStarts[$i];
	    if ($new_coverage > $recent_coverage) {
		$MaximalHitSet[$num_top_hits-1] = $i;
		$top_recent_hit = $i;
	    }

	} else {

	    # No overlap means new top!
	    push(@MaximalHitSet,$i);
	    $top_recent_hit = $i;	
	    $num_top_hits++;
    
	}
    }

    # Swag!  Now that we have our maximal set, what does it *mean*?
    my $exon_list_str  = $TargetHitExons[$MaximalHitSet[0]];
    my $amino_list_str = $TargetHitStarts[$MaximalHitSet[0]].'..'.$TargetHitEnds[$MaximalHitSet[0]];
    my $hit_list_str   = $MaximalHitSet[0];
    for (my $i=1; $i<$num_top_hits; $i++) {
	$exon_list_str  = $exon_list_str.','.$TargetHitExons[$MaximalHitSet[$i]];
	$amino_list_str = $amino_list_str.','.$TargetHitStarts[$MaximalHitSet[$i]].'..'.$TargetHitEnds[$MaximalHitSet[$i]];
	$hit_list_str   = $hit_list_str.','.$MaximalHitSet[$i];
    }

    # Pass back up the list of exon indices and corresponding amino coverage
    return($exon_list_str,$amino_list_str,$hit_list_str);
    
}






############################################################
#
#  Function: GetAminoHitGaps
#
sub GetAminoHitGaps
{
    my $amino_ranges_str = shift;
    my $seq_len = shift;

    my @AminoRanges = split(/\,/,$amino_ranges_str);
    my $num_ranges  = scalar(@AminoRanges);

    my @GapStarts;
    my @GapEnds;
    my $num_gaps = 0;

    # If the initial range starts at a position higher than 1 we'll need to
    # kick off the gaps with a starter gap
     if ($AminoRanges[0] !~ /^1\.\./) {
	$GapStarts[0] = 1;
	$AminoRanges[0] =~ /^(\d+)\.\./;
	$GapEnds[0] = $1-1;
	$num_gaps++;
    }
    
    # Run along the maximal ranges and see how much of the sequence is left
    # uncovered.
    for (my $i=0; $i<$num_ranges-1; $i++) {
	$AminoRanges[$i] =~ /\.\.(\d+)$/;
	$GapStarts[$num_gaps] = $1+1;
	$AminoRanges[$i+1] =~ /^(\d+)\.\./;
	$GapEnds[$num_gaps] = $1+1;
	$num_gaps++;
    }
    
    # If the last range doesn't end with the length of the protein
    # sequence, there's one more gap to consider!
    if ($AminoRanges[$num_ranges-1] !~ /\.\.$seq_len$/) {
	$AminoRanges[$num_ranges-1] =~ /\.\.(\d+)$/;
	$GapStarts[$num_gaps] = $1+1;
	$GapEnds[$num_gaps] = $seq_len;
	$num_gaps++;
    }

    return(\@GapStarts,\@GapEnds);
    
}





############################################################
#
#  Function: SaveTopXWInputs
#
sub SaveTopXWInputs
{
    my $xw_infname = shift;

    # Get the name without any of the directory path stuff
    $xw_infname =~ /\/([^\/]+)$/;
    my $fname = $1;

    # Break the name up into its components
    $fname =~ /^(\S+)\-(\d+)\.([^\.]+)\.weaver\.in/;
    my $gene  = $1;
    my $seqid = $2;
    my $chr   = $3;

    # We'll change over to the representation I prefer
    if ($chr =~ /\-revcomp/) {
	$chr =~ s/\-revcomp//;
	$chr = $chr.'[revcomp]';
    }

    # What are we going to call the file we're storing these inputs in?
    my $save_name = $seq_dirname.$seqid.'.partial.tmp';

    # Open up the save file and write out the name of the canonical chromosome
    my $SaveFile = OpenOutputFile($save_name);
    print $SaveFile "Chromosome: $chr\n\n";

    # Next, open up the original file and copy over all of the important stuff
    my $inf = OpenInputFile($xw_infname);
    while (my $line = <$inf>) {
	if ($line =~ /^Hit Score/) {
	    print $SaveFile "$line"; # Hit Score
	    $line = <$inf>;
	    print $SaveFile "$line"; # Amino Range
	    $line = <$inf>;
	    print $SaveFile "$line"; # Nucl Range
	    $line = <$inf>;
	    print $SaveFile "$line"; # Nucleotides
	    print $SaveFile "\n";
	}
    }
    close($inf);
    
    close($SaveFile);

    # Too easy!

}





############################################################
#
#  Function: RecordMaximalHits
#
sub RecordMaximalHits
{
    my $infname = shift;
    my $outfname = shift;
    my $hit_index_list = shift;

    # Open up the original ExonWeaver output file and the file we'll be transferring
    # a subset of those hits into
    my $inf  = OpenInputFile($infname);
    my $outf = OpenOutputFile($outfname);

    # ExonWeaver output files always start with a blank line, and who am I to
    # question that?
    print $outf "\n";

    # Which hits are we copying?
    my @HitIndices = split(/\,/,$hit_index_list);

    # DO THAT TRANSFER, BABY!
    my $hit_num = 0;
    my $exon_num = 0;
    while (my $line = <$inf>) {

	# I think this is fine...
	while ($line =~ /\-\-\-\-\- Exons/) {

	    if ($exon_num == $HitIndices[$hit_num]) {

		# This is one of the hits we love!
		print $outf "$line";
		$line = <$inf>;
		while ($line !~ /\-\-\-\-\- Exons/) {
		    print $outf "$line";
		    last if (feof($inf));
		    $line = <$inf>;
		}

		# Hit done!
		$hit_num++;
		
	    } else {
		$line = <$inf>; # Need to prime to escape the loop
	    }

	    # Have we picked up all the hits we set out to?
	    last if ($hit_num > scalar(@HitIndices));

	    # Nope, but we can at least move on to the next exon!
	    $exon_num++;
	    
	}

    }

    close($outf);
    close($inf);
    
}





############################################################
#
#  Function: RunBlatOnFileSet
#
sub RunBlatOnFileSet
{

    my $seqfiles_ref  = shift;
    my $genome        = shift;
    my $gene_name_fmt = shift; # 'seqname' or 'filename'

    my @SeqFileNames = @{$seqfiles_ref};

    # Since there's going to be a little bit of file I/O at the start,
    # we'll give a quick update
    DispProgQuilter('blatprep');
    
    # Get your mind out of the gutter, the 'cum' stands for 'cumulative'
    #
    # ... wink
    #
    my $cum_blat_fname = $seq_dirname.'BLAT.fa';
    my $CumBlatFile    = OpenOutputFile($cum_blat_fname);
    my %BlatNameGuide;
    my %BlatGenes;

    # We'll keep track of the number of gene families being supplied to BLAT,
    # so we can divide the workload between threads after BLAT's done.
    my $num_blat_genes = 0;
    my $num_blat_seqs  = 0; # NOTE that we'll name sequences from 1..num_blat_seqs

    # Run through the sequences that we want to concatenate for our BLAT run
    foreach my $fname (@SeqFileNames) {
    
	# Open up the file
	my $BlatFile = OpenInputFile($fname);

	# As we parse this file, we'll want to track what gene family each sequence
	# belongs to, along with its name
	my ($gene,$seqname);
	
	# If we're getting the gene name from the filename, it's regex time!
	if ($gene_name_fmt eq 'filename') {
	    $fname =~ /\/([^\/]+)\.fa/;
	    $gene  = $1;
	}

	# Transfer its contents into the big file, renaming sequences as needed
	while (my $line = <$BlatFile>) {
	    
	    $line =~ s/\n|\r//g;

	    # Depending on whether we're getting the gene family from the filename
	    # or the sequence name we'll parse 'new sequence' lines differently
	    if ($line =~ /^\>(.+)$/) {
		
		$seqname = $1;
		if ($gene_name_fmt eq 'seqname') {
		    $seqname =~ /^([^\|]+)\|([^\|]+)$/;
		    $gene    = $1;
		    $seqname = $2;
		}

		$num_blat_seqs++;
		print $CumBlatFile ">$num_blat_seqs\n";

		# Who are you, again?
		$BlatNameGuide{$num_blat_seqs} = $seqname.'&'.$gene;
		$BlatGenes{$gene} = 1;

	    } else {

		# Just a regular ol' sequence line
		print $CumBlatFile "$line\n";
		
	    }
	}

	close($BlatFile);

	# If we've indicated that we're getting our information from the sequence
	# name, that means that these files have been generated by Quilter, and
	# so we'll want to clear them
	RunSystemCommand("rm \"$fname\"") if ($gene_name_fmt eq 'seqname');

    }

    close($CumBlatFile);

    # Right on! Now it's time to run BLAT on our big concatenated file!
    # In order to do that, we'd better know where to find BLAT...
    my $BLAT;
    my $UnameCmd = OpenSystemCommand('uname -a |');
    my $uname = <$UnameCmd>;
    close($UnameCmd);
    if    (uc($uname) =~ /^LINUX /)  { $BLAT = $srcdir.'../inc/blat/blat.linux.x86_64';  }
    elsif (uc($uname) =~ /^DARWIN /) { $BLAT = $srcdir.'../inc/blat/blat.macOSX.x86_64'; }
    else                             { $BLAT = $srcdir.'../inc/blat/blat.macOSX.i386';   }

    # I think that now is the time for honesty... We're running BLAT!
    DispProgQuilter('blatrunning');
    
    # What do you want your output file to be named? Just kidding, I'll
    # figure it out for you.
    my $blat_outfname = $cum_blat_fname;
    $blat_outfname =~ s/\.fa$/\.out/;

    # Assemble and run the BLAT command!
    my $blat_cmd  = $BLAT.' -tileSize=7 -minIdentity=90 -maxIntron=1';
    $blat_cmd     = $blat_cmd.' -t=dnax -q=prot -out=blast8 -minScore=40';
    $blat_cmd     = $blat_cmd.' 1>/dev/null 2>&1';
    $blat_cmd     = $blat_cmd.' '.$genome.' '.$cum_blat_fname.' '.$blat_outfname;
    RunSystemCommand($blat_cmd);

    # Once that's all over with, we can go ahead and clear the cumulative file
    RunSystemCommand("rm \"$cum_blat_fname\"");
    
    # We'll return the BLAT outfile, the name guide and the number of families BLAT'd
    return($blat_outfname,\%BlatNameGuide,\%BlatGenes);

}








############################################################
#
#  Function: GenBlatMaps
#
sub GenBlatMaps
{
    my $blat_outfname      = shift;
    my $blat_nameguide_ref = shift;
    my $blat_genes_ref     = shift;
    my $num_cpus           = shift;

    my %BlatNameGuide = %{$blat_nameguide_ref};
    my %BlatGenes = %{$blat_genes_ref};

    # If we have fewer genes than cpus we'll reduce our threading accordingly
    my @BlatGeneList = keys %BlatGenes;
    my $num_blat_genes = scalar(@BlatGeneList);
    $num_cpus = Min($num_cpus,$num_blat_genes);

    # Spin off all our happy friends!
    my $threadID = SpawnProcesses($num_cpus);

    # Which range of gene name indices is this thread tasked with?
    my $gene_index_start =  $threadID    * int($num_blat_genes / $num_cpus);
    my $gene_index_end   = ($threadID+1) * int($num_blat_genes / $num_cpus);
    if ($threadID == $num_cpus-1) { $gene_index_end = $num_blat_genes; }

    # Make a hash of genes that this thread is in charge of for quick reference
    my %ThreadGenes;
    for (my $i=$gene_index_start; $i<$gene_index_end; $i++) {
	$ThreadGenes{$BlatGeneList[$i]} = 1;
    }

    # Run through the BLAT output file and pull in all the lines pertaining to
    # your sequences, little thread!
    my $blatf = OpenInputFile($blat_outfname);
    my %BlatData;
    while (my $line = <$blatf>) {

	$line =~ s/\n|\r//g;
	if ($line =~ /^(\S+)/) {

	    # Grab the original naming info for this sequence
	    my $seq_num = $1;
	    $BlatNameGuide{$seq_num} =~ /^([^\&]+)\&([^\&]+)$/;
	    my $seqname = $1;
	    my $gene = $2;

	    # If this isn't one of ours, leave it for a different friend to play with
	    next if (!$ThreadGenes{$gene});
	    
	    # If we have partial information from ExonWeaver (derived from GTF)
	    # then our search might only pertain to certain parts of the actual
	    # protein sequence -- check if this is the case!
	    my $partial = '-';
	    if ($seqname =~ /^(\d+)(s\d+e\d+)$/) {
		$seqname = $1.'-partial'; # Just for quick recognition...
		$partial = $2;
	    }

	    # We no longer need the indirection of the numeric naming for BLAT
	    $line =~ s/^\S+//;

	    # Record this hit
	    my $data = $seqname.' '.$partial.' '.$line;
	    if ($BlatData{$gene}) { $BlatData{$gene} = $BlatData{$gene}.'&'.$data; }
	    else                  { $BlatData{$gene} = $data;                      }

	}
    }
    close($blatf);

    # We now have all of our BLAT output organized by gene family.
    #
    # Now we'll break each family into its individual sequences and see if we can get
    # some full mappings.
    #
    # NOTE that if there are any partials for a sequence we'll treat them as a wholly
    # different sequence
    my $genes_completed = 0;
    foreach my $gene (keys %BlatData) {

	my %BlatBySeq;
	foreach my $datum (split(/\&/,$BlatData{$gene})) {
	    $datum =~ /^(\S+) (.+)$/;
	    my $seqname = $1;
	    my $hitdata = $2;
	    if ($BlatBySeq{$seqname}) {
		$BlatBySeq{$seqname} = $BlatBySeq{$seqname}.'&'.$hitdata;
	    } else {
		$BlatBySeq{$seqname} = $hitdata;
	    }
	}

	# Read in this gene's sequences. We'll store to a hash because we don't expect
	# to go through the sequences in the same order they were written.
	my %Seqs;
	my $inf = OpenInputFile($seq_dirname.$gene.'.fa');
	my $seqname;
	while (my $line = <$inf>) {

	    $line =~ s/\n|\r//g;
	    next if (!$line);

	    if ($line =~ /\>(\S+)/) {
		$seqname = $1;
	    } else {
		$line = uc($line);
		if ($Seqs{$seqname}) { $Seqs{$seqname} = $Seqs{$seqname}.$line; }
		else                 { $Seqs{$seqname} = $line;                 }
	    }
	    
	}
	close($inf);

	# Now we can go sequence-by-sequence (again, treating partials specially)
	# looking for a full mapping.
	# We'll want to see partials before treating the full sequence, so we'll
	# flip the sort order of the sequence names.
	my @SeqNames = sort { $b cmp $a } keys %BlatBySeq;
	my @FullMaps;
	my $num_full_maps = 0;
	for (my $i=0; $i<scalar(@SeqNames); $i++) {

	    # Get the sequence name and pull its BLAT data
	    my $seqname  = $SeqNames[$i];
	    my @BlatHits = split(/\&/,$BlatBySeq{$seqname});
	    
	    if ($seqname =~ /\-partial$/) {

		# Revert to the original sequence name
		$seqname =~ s/\-partial$//;

		# Has BLAT given us the power to fill in the gaps in this sequence?
		my $seq = $Seqs{$seqname};
		$FullMaps[$i] = AttemptBlatFill($seqname,$seq,\@BlatHits);

		# NOTE
		# NOTE
		# NOTEY
		# NOTEY
	        # NOTE: For now, if we fill a partial, we'll just go ahead and skip
		#       the non-partial work... BUT THIS MIGHT BE SOMETHING TO CHANGE
		$FullMaps[++$i] = 0 if ($FullMaps[$i]);
		$num_full_maps++; # This won't increment otherwise
		
	    } else {

		# Dive right on in with Spaln!
		my $seq = $Seqs{$seqname};
		$FullMaps[$i] = BlatToSpalnSearch($seqname,$seq,\@BlatHits);

	    }

	    # Did we get one?
	    $num_full_maps++ if ($FullMaps[$i]);

	}

	# Regardless of whether or not we found anything, we're done with this gene!
	$genes_completed++;
	DispProgQuilter('blat2spaln|'.$threadID.'|'.$genes_completed);

	# Well well welly well...
	# Do we have any full maps to add to our '.quilter.out' file?
	next if (!$num_full_maps);
	
	# WELL WELL WELLY WELL
	# LET'S GET APPENDING

	# We don't know if there's already an output file for this gene family, so
	# we'll use the appending output option
	my $outfname = $seq_dirname.$gene.'.quilter.out';
	open(my $outf,'>>',$outfname)
	    || die "\n  ERROR:  Failed to open output file '$outfname'\n\n";
	for (my $i=0; $i<scalar(@SeqNames); $i++) {
	    if ($FullMaps[$i]) {
		print $outf "$FullMaps[$i]\n";
	    }
	}
	close($outf);

    }

    # Friend time has come to its conclusion -- everyone goes home except the master
    if ($threadID) { exit(0); }
    while (wait() != -1) {}

    # I think... we're done?!?!?!?
    
}






############################################################
#
#  Function: AttemptBlatFill
#
sub AttemptBlatFill
{
    my $seqname = shift;
    my $seq_str = shift;
    my $blathits_ref = shift;

    my @BlatHits = @{$blathits_ref};

    # Find the file with the stored ExonWeaver inputs to this gene's canon chromosome
    my $xwinfname = $seq_dirname.$seqname.'.partial.tmp';

    # Open it up and pull in the original hits, along with the chromosome info
    my $OrigXWInf = OpenInputFile($xwinfname);

    my $orig_chr = <$OrigXWInf>;
    $orig_chr =~ s/\n|\r//g;
    $orig_chr =~ s/^Chromosome\: //;
    $orig_chr =~ s/\[revcomp\]/\-revcomp/;

    my @OrigHits;
    while (my $line = <$OrigXWInf>) {

	if ($line =~ /Hit Score/) {

	    $line =~ /Hit Score   \: (\S+)/;
	    my $hit_score = $1;

	    $line = <$OrigXWInf>;
	    $line =~ /Amino Range \: (\d+)\.\.(\d+)/;
	    my $amino_start = $1;
	    my $amino_end   = $2;
	    
	    $line = <$OrigXWInf>;
	    $line =~ /Nucl Range  \: (\d+)\.\.(\d+)/;
	    my $nucl_start = $1;
	    my $nucl_end   = $2;
	    
	    $line = <$OrigXWInf>;
	    $line =~ /Nucleotides \: (\S+)/;
	    my $nucl_str = $1;

	    my $hitstr = $amino_start.'|'.$amino_end.'|'.$nucl_start.'|'.$nucl_end;
	    $hitstr    = $hitstr.'|-|-|'.$hit_score.'|'.$nucl_str;
	    push(@OrigHits,$hitstr);
	    
	}
	
    }
    close($OrigXWInf);

    # We'll clear this file out, since it's not going to be used for anything else
    RunSystemCommand("rm \"$xwinfname\"");

    # Next up, we'll parse each of the lines from our BLAT hit, and assume that it
    # nailed the ranges exactly right...
    my %ExtractedHits;
    foreach my $blathit (@BlatHits) {

	# Recall that we've added details about the nature of our partial hit to
	# the front of the hit.
	$blathit =~ /^\S+\s+(\S+)\s+(.*)$/;
	my $partial_info = $1;
	my $blat_outline = $2;
	my ($blat_chr,$blat_amino_start,$blat_amino_end,$blat_nucl_start,$blat_nucl_end,$blat_score)
	    = ParseBlatLine($blat_outline);
	
	# We'll check to make sure there don't appear to be any gaps
	if ((1+abs($blat_nucl_end-$blat_nucl_start))/3 != 1+$blat_amino_end-$blat_amino_start) {
	    next;
	}

	# We may have a mismatch here or there, but no gaps! Pull in the sequence
	# (with an additional 17 buffer nucleotides) and we'll be solid!
	my $sfetch_range;
	if ($blat_nucl_start > $blat_nucl_end) {
	    my $sfetch_range_start = $blat_nucl_start + 17;
	    my $sfetch_range_end   = $blat_nucl_end   - 17;
	    $sfetch_range = $sfetch_range_start.'..'.$sfetch_range_end;
	} else {
	    my $sfetch_range_start = $blat_nucl_start - 17;
	    my $sfetch_range_end   = $blat_nucl_end   + 17;
	    $sfetch_range = $sfetch_range_start.'..'.$sfetch_range_end;
	}

	my $sfetchcmd = $sfetch." -range $sfetch_range \"$genome\" \"$blat_chr\"";
	my $inf = OpenSystemCommand($sfetchcmd);
	my $line = <$inf>; # Eat the '>' line
	my $blat_nucl_str = '';
	while (my $line = <$inf>) {
	    $line =~ s/\n|\r//g;
	    $blat_nucl_str = $blat_nucl_str.uc($line);
	}
	close($inf);

	# Now we're past the point of using the chromosome for actual data extraction,
	# so we can change the name to be more human-informative.
	# NOTE that we'll be passing this off to ExtractChimericXWMap, which
	# anticipates '-revcomp' rather than '[revcomp]'
	$blat_chr = $blat_chr.'-revcomp' if ($blat_nucl_start > $blat_nucl_end);

	my $hitstr = $blat_amino_start.'|'.$blat_amino_end.'|'.$blat_nucl_start;
	$hitstr    = $hitstr.'|'.$blat_nucl_end.'|-|-|'.$blat_score.'|'.$blat_nucl_str;
	$hitstr    = $hitstr.'|'.$blat_chr;

	if ($ExtractedHits{$blat_amino_start}) {
	    $ExtractedHits{$blat_amino_start} = $ExtractedHits{$blat_amino_start}.'&'.$hitstr;
	} else {
	    $ExtractedHits{$blat_amino_start} = $hitstr;
	}
	
    }


    # We're nearing a position where we can weave some hits!  We'll just need to
    # merge our original hits and BLAT-based hits so that they're organized by
    # start amino acid position
    my @BlatStarts = sort { $a <=> $b } keys %ExtractedHits;
    return 0 if (scalar(@BlatStarts) == 0);

    my $orig_exon_pos = 0;
    my $blat_list_pos = 0;
    my @ChrsByExon; # Keeping with the nomenclature of 'ExtractChimericXWMap'
    my @FullHitList;
    while ($orig_exon_pos < scalar(@OrigHits) && $blat_list_pos < scalar(@BlatStarts)){
	if ($OrigHits[$orig_exon_pos] <= $BlatStarts[$blat_list_pos]) {
	    push(@FullHitList,$OrigHits[$orig_exon_pos]);
	    push(@ChrsByExon,$orig_chr);
	    $orig_exon_pos++;
	} else {
	    foreach my $blathit (split(/\&/,$ExtractedHits{$BlatStarts[$blat_list_pos]})) {
		# Pull this hit's chromosome
		$blathit =~ /^(\S+)\|([^\|]+)$/;
		my $blat_data = $1;
		my $blat_chr  = $2;
		push(@FullHitList,$blat_data);
		push(@ChrsByExon,$blat_chr);
	    }
	    $blat_list_pos++;
	}
    }
    while ($orig_exon_pos < scalar(@OrigHits)) {
	push(@FullHitList,$OrigHits[$orig_exon_pos]);
	push(@ChrsByExon,$orig_chr);
	$orig_exon_pos++;
    }
    while ($blat_list_pos < scalar(@BlatStarts)){
	foreach my $blathit (split(/\&/,$ExtractedHits{$BlatStarts[$blat_list_pos]})) {
	    # Pull this hit's chromosome
	    $blathit =~ /^(\S+)\|([^\|]+)$/;
	    my $blat_data = $1;
	    my $blat_chr  = $2;
	    push(@FullHitList,$blat_data);
	    push(@ChrsByExon,$blat_chr);
	}
	$blat_list_pos++;
    }

    # Now we can convert our hits into ExonWeaver input! We'll reuse the filename
    # that stored the original hits.
    my $num_exons  = scalar(@FullHitList);
    my $seq_len    = length($seq_str);
    my $WeaverFile = OpenOutputFile($xwinfname);
    print $WeaverFile "Num Hits : $num_exons\n";
    print $WeaverFile "Seq Len  : $seq_len\n";
    print $WeaverFile "Sequence : $seq_str\n";
    foreach my $hit (@FullHitList) {
	print $WeaverFile "\n";
	PrintHitToWeaverInf($hit,$WeaverFile);
    }
    close($WeaverFile);

    # Come up with a hip name for the output file
    my $xwoutfname = $xwinfname;
    $xwoutfname =~ s/\.tmp$/\.out/;

    # NOTE: Even though it's possible that our BLAT results will only have found
    #       hits to a single chromosome, we'll use the chimeric search to cover all
    #       possible bizzare splicing patterns.
    my $weaver_cmd = $srcdir."ExonWeaver --allow-inconsistency \"$xwinfname\" > \"$xwoutfname\"";
    RunSystemCommand($weaver_cmd);

    # Begone, input file!
    RunSystemCommand("rm \"$xwinfname\"");

    # Did we get anything at all?
    if (!(-s $xwoutfname)) {
	RunSystemCommand("rm \"$xwoutfname\"") if (-e $xwoutfname);
	return 0;
    }

    # You know the way this goes
    my $full_map = CheckForFullMaps($xwoutfname,$seq_len);
    if ($full_map) {
	$full_map = ExtractChimericXWMap($xwoutfname,\@ChrsByExon,$seq_len,$seqname);
	# Credit to BLAT!
	$full_map =~ s/FastMap2/BLAT\+FastMap2/;
    }
    
    # That's all!
    RunSystemCommand("rm \"$xwoutfname\"");
    return $full_map;
    
}






############################################################
#
#  Function: BlatToSpalnSearch
#
sub BlatToSpalnSearch
{
    my $seqname = shift;
    my $seq_str = shift;
    my $blathits_ref = shift;

    my @BlatHits = @{$blathits_ref};

    # We'll go ahead and kick things off by making a file with the protein sequence
    my $prot_fname = $seq_dirname.$seqname.'.blat2spaln.prot.in';
    my $ProtFile = OpenOutputFile($prot_fname);
    print $ProtFile ">$seqname\n";
    my @Seq = split(//,$seq_str);
    for (my $i=0; $i<scalar(@Seq); $i++) {
	print $ProtFile "$Seq[$i]";
	print $ProtFile "\n" if (($i+1) % 60 == 0);
    }
    print $ProtFile "\n";
    close($ProtFile);

    # Because our search function needs a list of protein sequences to search,
    # we have to do this dumb stuff
    my @ProtFnames;
    push(@ProtFnames,$prot_fname);
    my @SeqNames;
    push(@SeqNames,$seqname);
    my @SeqStrs;
    push(@SeqStrs,$seq_str);


    # Now we'll organize our hits into groups along the length of the protein.
    #
    # The idea is that each position of the protein can belong to one or more groups,
    # where a group is defined by overlapping hits to the same chromosome.
    #
    # Our groups will run [1..num_hit_groups]
    #
    my $seq_len = length($seq_str);
    my @Groupings;
    for (my $i=0; $i<$seq_len; $i++) {
	$Groupings[$i] = 0;
    }
    my @GroupHits; # The index in BlatHits
    my @GroupStarts;
    my @GroupEnds;
    my @GroupChrs;
    my $num_groups = 0;
    for (my $hit_id=0; $hit_id<scalar(@BlatHits); $hit_id++) {

	# We'll end up clearing out the first pieces of data in the blat hits,
	# since we know there aren't partial hits in this list.
	$BlatHits[$hit_id] =~ /^\-\s+(.*)$/;
	my $blat_outline = $1;
	my ($chr,$amino_start,$amino_end,$nucl_start,$nucl_end,$score)
	    = ParseBlatLine($blat_outline);

	# We'll want to know if this is the reverse complement
	if ($nucl_start > $nucl_end) { $chr = $chr.'-'; }
	else                         { $chr = $chr.'+'; }

	# Let's go ahead and convert the amino coordinates into easy mode
	$amino_start--;
	$amino_end--;

	# Revisionist hit-story!
	my $hit = $chr.'|'.$amino_start.'|'.$amino_end.'|'.$nucl_start.'|'.$nucl_end;
	$hit    = $hit.'|'.$score.'|'.$hit_id;
	$BlatHits[$hit_id] = $hit;

	# Does this hit overlap with any other hits to the same chromosome?
	my $overlap = 0;
	for (my $i=$amino_start; $i<=$amino_end; $i++) {
	    if ($Groupings[$i]) {
		foreach my $group (split(/\,/,$Groupings[$i])) {
		    if ($GroupChrs[$group] eq $chr) {
			$overlap = $group;
			last;
		    }
		}
	    }
	    last if ($overlap);
	}
	
	# Find your crew, lil' hit
	my $group;
	if ($overlap) {
	    
	    $group = $overlap;
	    $GroupHits[$group] = $GroupHits[$group].','.$hit_id;
	    $GroupStarts[$group] = Min($GroupStarts[$group],$amino_start);
	    $GroupEnds[$group] = Max($GroupEnds[$group],$amino_end);
	    
	} else {
	    
	    # First to lay claim (from this chromosome)!
	    $num_groups++;
	    $GroupHits[$num_groups] = $hit_id;
	    $GroupChrs[$num_groups] = $chr;
	    $GroupStarts[$num_groups] = $amino_start;
	    $GroupEnds[$num_groups] = $amino_end;
	    $group = $num_groups;
	    
	}

	# Make sure every position is appropriately marked for this group!
	for (my $i=$amino_start; $i<=$amino_end; $i++) {
	    if (!$Groupings[$i]) {
		$Groupings[$i] = $group;
	    } elsif ($Groupings[$i] !~ /^$group|\,$group/) {
		$Groupings[$i] = $Groupings[$i].','.$group;
	    }
	}
	
    }


    # Next up, because it's possible that groups could have merged, we'll
    # do a quick scan through the array to merge overlapping groups from the
    # same chromosome
    my %StartsToGroups;
    my $final_num_groups = 0;
    for (my $group=1; $group<=$num_groups; $group++) {

	# Which chromosome is this?
	my $chr = $GroupChrs[$group];
	next if (!$chr);

	for (my $check=$group+1; $check<=$num_groups; $check++) {

	    next if ($GroupChrs[$check] ne $chr);

	    # OOOOOH, same chromosome! Do they overlap at all?

	    # NOTE: While this 'gap' seems large, the way to think about it
	    #  (i.m.o.) is "What size of Blat hit to another chromosome would
	    #  be enough to consider interrupting this run?"
	    my $maxgap = 25;
	    if (($GroupStarts[$group] >= $GroupStarts[$check] - $maxgap &&
		 $GroupStarts[$group] <= $GroupEnds[$check] + $maxgap)
		||
		($GroupEnds[$group] >= $GroupStarts[$check] - $maxgap &&
		 $GroupEnds[$group] <= $GroupEnds[$check] + $maxgap)) {

		# OVERLAP!!!! Combine those groups!
		$GroupStarts[$group] = Min($GroupStarts[$group],$GroupStarts[$check]);
		$GroupEnds[$group] = Max($GroupEnds[$group],$GroupEnds[$check]);
		$GroupHits[$group] = $GroupHits[$group].','.$GroupHits[$check];
		$GroupChrs[$check] = 0;

		# So, this is annoying, but we'll have to go back and make
		# sure that we don't now transitively overlap with an earlier group...
		$check = $group;
		
	    }
	    
	}

	# Switch over to the 'final_num_groups' index
	$GroupStarts[$final_num_groups] = $GroupStarts[$group];
	$GroupEnds[$final_num_groups] = $GroupEnds[$group];
	$GroupChrs[$final_num_groups] = $GroupChrs[$group];
	$GroupHits[$final_num_groups] = $GroupHits[$group];

	# Good stuff! Now we'll record the start position of this group, so we can
	# go through the hits in sorted order along the protein.
	if ($StartsToGroups{$GroupStarts[$final_num_groups]}) {
	    $StartsToGroups{$GroupStarts[$final_num_groups]}
	    = $StartsToGroups{$GroupStarts[$final_num_groups]}.','.$final_num_groups;
	} else {
	    $StartsToGroups{$GroupStarts[$final_num_groups]} = $final_num_groups;
	}
	
	# That's a keeper!
	$final_num_groups++;
	
    }


    # Radical! Now, we'll go through all of our groups and determine their
    # search ranges on the chromosome.
    my @GroupNuclRanges;
    my @GroupSorting;
    my %ChrGroupSorting;
    my %CoverageByChr;
    foreach my $start_amino (sort {$a <=> $b} keys %StartsToGroups) {
	foreach my $group (split(/\,/,$StartsToGroups{$start_amino})) {

	    # And that's your place!
	    push(@GroupSorting,$group);

	    # But this is also your place!
	    my $chr = $GroupChrs[$group];
	    if ($ChrGroupSorting{$chr}) {
		$ChrGroupSorting{$chr} = $ChrGroupSorting{$chr}.','.$group;
	    } else {
		$ChrGroupSorting{$chr} = $group;
	    }

	    # While we're at it (wrt bookkeeping), how much of the protein do the
	    # hits to this chromosome cover?
	    my $group_coverage = $GroupEnds[$group] - $GroupStarts[$group] + 1;
	    if ($CoverageByChr{$chr}) { $CoverageByChr{$chr} += $group_coverage; }
	    else                      { $CoverageByChr{$chr}  = $group_coverage; }

	    # Make a list of all nucleotide start and end positions
	    my @NuclStarts;
	    my @NuclEnds;
	    foreach my $hit_id (split(/\,/,$GroupHits[$group])) {
		my $hit = $BlatHits[$hit_id];
		$hit =~ /^[^\|]+\|\d+\|\d+\|(\d+)\|(\d+)\|/;
		push(@NuclStarts,$1);
		push(@NuclEnds,$2);
	    }

	    # Perform a series of clusterings until the span of each set of clustered
	    # hits is less than a preset distance
	    my $max_span = 1000000; # 1 Mb	
	    my $ranges   = ClusterNuclRanges(\@NuclStarts,\@NuclEnds,$max_span);
	    $GroupNuclRanges[$group] = $ranges;

	}
    }


    

    ##########################
    #                        #
    #    S E A R C H    1    #
    #                        #
    ##########################
    
    
    # Swell! Now we'll go through each of our chromosomes, and any that have more
    # than 75% coverage will be Spalned
    my $top_pct_id  = 0;
    my $top_hit_str = 0;
    foreach my $chr (keys %CoverageByChr) {

	# Do you have solid coverage?
	next if ($CoverageByChr{$chr} < 3 * $seq_len / 4);

	# Woot! So, if your coverage is so great, where's it coming from?
	my @RangeStarts;
	my @RangeEnds;
	foreach my $group (split(/\,/,$ChrGroupSorting{$chr})) {
	    foreach my $range (split/\,/,$GroupNuclRanges[$group]) {
		$range =~ /^(\d+)\.\.(\d+)$/;
		push(@RangeStarts,$1);
		push(@RangeEnds,$2);
	    }
	}

	# We'll need to know the chromosome's true identity
	my $revcomp = 0;
	$revcomp = 1 if ($chr =~ /\-$/);
	my $true_chr = $chr;
	$true_chr =~ s/\S$//;
	my $chr_len = $ChrSizes{$true_chr};

	# First off, we'll try a spaln search on a range that assumes no funny business
	# vis-a-vis splicing.
	my ($starts_ref,$ends_ref,$num_ranges,$sum_len)
	    = GetSpalnSfetchRanges(\@RangeStarts,\@RangeEnds,$chr_len,1);
	my @SpalnStarts = @{$starts_ref};
	my @SpalnEnds   = @{$ends_ref};
	
	# We won't do a search that would require pulling in >15Mb
	next if ($sum_len > 15000000);

	# Because we're using a function that plays friendly with a multiple-chromosome
	# version of this index building, we'll need to note the chromosomes for each
	# range
	my @SpalnChrs;
	for (my $i=0; $i<$num_ranges; $i++) {
	    $SpalnChrs[$i] = $chr;
	}

	# Let's see what we see!
	my ($hit_strs_ref,$hit_pcts_id_ref)
	    = SpalnSearch(\@SeqNames,\@SeqStrs,\@ProtFnames,\@SpalnStarts,\@SpalnEnds,\@SpalnChrs,90.0);
	my $hit_str = @{$hit_strs_ref}[0];
	my $hit_pct_id = @{$hit_pcts_id_ref}[0];

	# If we got a hit, then we'll be happy about that!
	if ($hit_pct_id > $top_pct_id) {
	    $top_pct_id  = $hit_pct_id;
	    $top_hit_str = $hit_str;
	}
    }

    # Any chance we got a hit on easy mode?
    if ($top_pct_id) {
	RunSystemCommand("rm \"$prot_fname\"");
	return $top_hit_str;
    }


    ##########################
    #                        #
    #    S E A R C H    2    #
    #                        #
    ##########################
    
    # Alright, we've tried the 'standard' approach to splicing -- now let's see if
    # there's anything out-of-order, but just on a chromosomal level...
    foreach my $chr (keys %CoverageByChr) {

	# Do you have solid coverage?
	next if ($CoverageByChr{$chr} < 3 * $seq_len / 4);

	# Woot! So, if your coverage is so great, where's it coming from?
	my @RangeStarts;
	my @RangeEnds;
	foreach my $group (split(/\,/,$ChrGroupSorting{$chr})) {
	    foreach my $range (split/\,/,$GroupNuclRanges[$group]) {
		$range =~ /^(\d+)\.\.(\d+)$/;
		push(@RangeStarts,$1);
		push(@RangeEnds,$2);
	    }
	}

	# We'll need to know the chromosome's true identity
	my $revcomp = 0;
	$revcomp = 1 if ($chr =~ /\-$/);
	my $true_chr = $chr;
	$true_chr =~ s/\S$//;
	my $chr_len = $ChrSizes{$true_chr};

	# First off, we'll try a spaln search on a range that assumes no funny business
	# vis-a-vis splicing.

	# Hmm, doesn't look like we got a solid hit the first time.
	# Let's see if there's a funkiness to the order (weird splicing alert!)
	my ($starts_ref,$ends_ref,$num_ranges,$sum_len)
	    = GetSpalnSfetchRanges(\@RangeStarts,\@RangeEnds,$chr_len,0);
	my @SpalnStarts = @{$starts_ref};
	my @SpalnEnds   = @{$ends_ref};

	# We still won't do a search that would require pulling in >15Mb
	next if ($sum_len > 15000000);

	my @SpalnChrs;
	for (my $i=0; $i<$num_ranges; $i++) {
	    $SpalnChrs[$i] = $chr;
	}

	# Let's try this again!
	my ($hit_strs_ref,$hit_pcts_id_ref)
	    = SpalnSearch(\@SeqNames,\@SeqStrs,\@ProtFnames,\@SpalnStarts,\@SpalnEnds,\@SpalnChrs,90.0);
	my $hit_str = @{$hit_strs_ref}[0];
	my $hit_pct_id = @{$hit_pcts_id_ref}[0];

	# Last call for spaln-ohol!
	if ($hit_pct_id > $top_pct_id) {
	    $top_pct_id  = $hit_pct_id;
	    $top_hit_str = $hit_str;
	}

    }

    
    # If we had a good hit, that's a good thing, don't you see?
    if ($top_pct_id) {
	RunSystemCommand("rm \"$prot_fname\"");
	return $top_hit_str;
    }


    ##########################
    #                        #
    #    S E A R C H    3    #
    #                        #
    ##########################
    
    # Hmmmmm.....
    # Alright, let's go wild 'n' crazy and see if we can jumble our chromosomes
    # together in some sort of pleasing fashion...
    my @RangeStarts;
    my @RangeEnds;
    my @RangeChrs;
    my $sum_len = 0;
    foreach my $group (@GroupSorting) {

	my @ChrStarts;
	my @ChrEnds;
	foreach my $range (split(/\,/,$GroupNuclRanges[$group])) {
	    $range =~ /^(\d+)\.\.(\d+)$/;
	    my $start = $1;
	    my $end = $2;
	    push(@ChrStarts,$start);
	    push(@ChrEnds,$end);
	}

	# We'll organize groups on a per-chromosome basis
	my $group_chr = $GroupChrs[$group];
	my $true_chr = $group_chr;
	$true_chr =~ s/\S$//;
	my $chr_len = $ChrSizes{$true_chr};
	
	my ($starts_ref,$ends_ref,$num_ranges) =
	    GetSpalnSfetchRanges(\@ChrStarts,\@ChrEnds,$chr_len,1);
	@ChrStarts = @{$starts_ref};
	@ChrEnds   = @{$ends_ref};

	for (my $i=0; $i<$num_ranges; $i++) {
	    $sum_len += abs($ChrEnds[$i]-$ChrStarts[$i]);
	    push(@RangeStarts,$ChrStarts[$i]);
	    push(@RangeEnds,$ChrEnds[$i]);
	    push(@RangeChrs,$group_chr);
	}
    }

    # Make sure we aren't pulling in an inordinate amount of sequence (15 Mb)...
    if ($sum_len > 15000000) {
	RunSystemCommand("rm \"$prot_fname\"");
	return 0;
    }

    # Last call for Spaln-ohol!
    my ($hit_strs_ref,$hit_pcts_id_ref)
	= SpalnSearch(\@SeqNames,\@SeqStrs,\@ProtFnames,\@RangeStarts,\@RangeEnds,\@RangeChrs,75.0);
    my $hit_str = @{$hit_strs_ref}[0];
    my $hit_pct_id = @{$hit_pcts_id_ref}[0];

    # However this breaks, it breaks with you in the bin!
    RunSystemCommand("rm \"$prot_fname\"");

    # And that's all there is!
    return 0 if (!$hit_pct_id);
    return $hit_str;
    
}







############################################################
#
#  Function: GetSpalnSfetchRanges
#
sub GetSpalnSfetchRanges
{
    my $starts_ref = shift;
    my $ends_ref = shift;
    my $chr_len = shift;
    my $condense = shift;

    my @RangeStarts = @{$starts_ref};
    my @RangeEnds   = @{$ends_ref};
    my $num_ranges  = scalar(@RangeStarts);

    my $strand = 1;
    $strand = -1 if ($RangeStarts[0] > $RangeEnds[0]);

    if ($condense) {
	# We'll run everything through our cluster-er first
	my $max_span = 1000000; # 1Mb
	my $ranges_str = ClusterNuclRanges(\@RangeStarts,\@RangeEnds,$max_span);
	my @Ranges = split(/\,/,$ranges_str);
	for (my $i=0; $i<scalar(@Ranges); $i++) {
	    $Ranges[$i] =~ /^(\d+)\.\.(\d+)$/;
	    $RangeStarts[$i] = $1;
	    $RangeEnds[$i] = $2;
	}
	$num_ranges = scalar(@Ranges);
    }
    
    # Give the boundaries some extra room to play
    $RangeStarts[0]           -= 50000 * $strand;
    $RangeEnds[$num_ranges-1] += 50000 * $strand;
    
    # We'll also just pull in 10Kb to each range...
    for (my $i=0; $i<$num_ranges; $i++) {
	$RangeStarts[$i] -= 10000 * $strand;
	$RangeEnds[$i]   += 10000 * $strand;
    }
    
    # Check for any overlapping ranges.  Note that it's possible for adjacent
    # ranges to not be in the expected order, so we can't do the easy implicit
    # check...
    if ($strand == 1) { # FWD (usually revcomp comes first)
	
	my @FinalStarts;
	my @FinalEnds;
	
	push(@FinalStarts,$RangeStarts[0]);
	push(@FinalEnds,$RangeEnds[0]);
	my $num_finals = 0; # We'll keep this 1 too low for now, for indexing
	
	for (my $i=1; $i<$num_ranges; $i++) {
	    if ($RangeStarts[$i] > $FinalStarts[$num_finals]
		&& $RangeStarts[$i] < $FinalEnds[$num_finals]) {
		$FinalEnds[$num_finals] = $RangeEnds[$i];
	    } else {
		push(@FinalStarts,$RangeStarts[$i]);
		push(@FinalEnds,$RangeEnds[$i]);
		$num_finals++;
	    }
	}
	
	# Copy it over!
	for (my $i=0; $i<=$num_finals; $i++) {
	    $RangeStarts[$i] = $FinalStarts[$i];
	    $RangeEnds[$i] = $FinalEnds[$i];
	}
	$num_ranges = $num_finals+1; # No longer 1 too low!
	
	# Make sure that you aren't overstepping any boundaries
	for (my $i=0; $i<$num_ranges; $i++) {
	    $RangeStarts[$i] = Max(1,$RangeStarts[$i]);
	    $RangeEnds[$i] = Min($RangeEnds[$i],$chr_len);
	}
	
    } else { # REV (usually revcomp comes first, like I already told you)
	
	my @FinalStarts;
	my @FinalEnds;
	
	push(@FinalStarts,$RangeStarts[0]);
	push(@FinalEnds,$RangeEnds[0]);
	my $num_finals = 0; # We'll keep this 1 too low for now, for indexing
	
	for (my $i=1; $i<$num_ranges; $i++) {
	    if ($RangeStarts[$i] < $FinalStarts[$num_finals]
		&& $RangeStarts[$i] > $FinalEnds[$num_finals]) {
		$FinalEnds[$num_finals] = $RangeEnds[$i];
	    } else {
		push(@FinalStarts,$RangeStarts[$i]);
		push(@FinalEnds,$RangeEnds[$i]);
		$num_finals++;
	    }
	}
	
	# Copy it over!
	for (my $i=0; $i<=$num_finals; $i++) {
	    $RangeStarts[$i] = $FinalStarts[$i];
	    $RangeEnds[$i] = $FinalEnds[$i];
	}
	$num_ranges = $num_finals+1; # No longer 1 too low!
	
	# Make sure that you aren't overstepping any boundaries
	for (my $i=0; $i<$num_ranges; $i++) {
	    $RangeStarts[$i] = Min($RangeStarts[$i],$chr_len);
	    $RangeEnds[$i] = Max(1,$RangeEnds[$i]);
	}
	
    }

    # Compute the sum length of the ranges
    my $sum_len = 0;
    for (my $i=0; $i<$num_ranges; $i++) {
	$sum_len += abs($RangeEnds[$i]-$RangeStarts[$i]);
    }

    # Return the adjusted ranges
    return (\@RangeStarts,\@RangeEnds,$num_ranges,$sum_len);

}







############################################################
#
#  Function: SpalnSearch
#
sub SpalnSearch
{
    my $seqnames_ref = shift;
    my $seq_strs_ref = shift;
    my $prot_fnames_ref = shift;
    my $range_starts_ref = shift;
    my $range_ends_ref = shift;
    my $range_chrs_ref = shift;
    my $min_pct_id = shift;

    my @SeqNames = @{$seqnames_ref};
    my @SeqStrings = @{$seq_strs_ref};
    my @ProtFnames = @{$prot_fnames_ref};
    
    my @RangeStarts = @{$range_starts_ref};
    my @RangeEnds = @{$range_ends_ref};
    my @RangeChrs = @{$range_chrs_ref};

    my $num_ranges = scalar(@RangeChrs);

    # We'll need to track where jump positions are in our Frankensteinian seq.
    # The format will be 'reported-pos':'actual-nucl-pos' so we scan until
    # 'reported-pos' is higher than the index we have, and then add the difference
    # to 'actual-nucl-pos'
    #
    # Also, note that if there's some out-of-order splicery indicated, this will
    # still take us to the right place!
    #
    my @JumpList;
    my $reported_pos = 1;
    for (my $i=0; $i<$num_ranges; $i++) {
	my $jump_info = $reported_pos.':'.$RangeStarts[$i].'|'.$RangeChrs[$i];
	push(@JumpList,$jump_info);
	$reported_pos += abs($RangeEnds[$i]-$RangeStarts[$i])+1;
    }
    
	
    
    # Wowee!  Time to pull in a bunch of nucleotides, concatenate them into one
    # great big sequence, and do a Spaln search!
    my $nucl_fname = $ProtFnames[0];
    $nucl_fname =~ s/\.prot\.in/\.nucl\.in/;
    my $NuclFile = OpenOutputFile($nucl_fname);
    print $NuclFile ">nucls\n";
    
    for (my $i=0; $i<$num_ranges; $i++) {
	
	my $range = $RangeStarts[$i].'..'.$RangeEnds[$i];
	my $chr = $RangeChrs[$i];
	$chr =~ s/\S$//;
	my $sfetchcmd = $sfetch." -range $range \"$genome\" \"$chr\"";
	my $inf = OpenSystemCommand($sfetchcmd);
	my $line = <$inf>; # eat header
	while ($line = <$inf>) {
	    $line =~ s/\n|\r//g;
	    next if (!$line);
	    print $NuclFile "$line\n";
	}
	close($inf);
	
    }
    
    # Ready, cap'n!
    close($NuclFile);

    my @SpalnHitStrs;
    my @SpalnPctIDs;
    for (my $i=0; $i<scalar(@ProtFnames); $i++) {

	my $seqname = $SeqNames[$i];
	my @Seq = split(//,$SeqStrings[$i]);
	my $prot_fname = $ProtFnames[$i];

	# Assemble the Spaln command!
	my $spaln_fname = $prot_fname;
	$spaln_fname =~ s/\.prot\.in/\.spaln\.out/;
	my $spaln_cmd = $spaln."\"$nucl_fname\" \"$prot_fname\" 1>\"$spaln_fname\" 2>/dev/null";

	# Sometimes Spaln doesn't like an input, so we don't bail if the
	# system call goes badly
	if (system($spaln_cmd)) {
	    push(@SpalnHitStrs,0);
	    push(@SpalnPctIDs,0);
	    next;
	}
    
	# Don't tell me we did all that for NOTHING?!
	if (!(-e $spaln_fname)) {
	    push(@SpalnHitStrs,0);
	    push(@SpalnPctIDs,0);
	    next;
	}
    
	# UGH, I wish I could quit you, Spaln parsing
	# Note that we'll open the file here, just so it's easier to return
	# if we run into parsing errors.
	my $SpalnFile = OpenInputFile($spaln_fname);
	my ($nucl_ranges_ref,$amino_ranges_ref,$centers_ref,$spaln_pct_id)
	    = ParseSpalnOutput($SpalnFile,\@Seq);
	    # = OPSO_ParseSpalnOutput($SpalnFile,\@Seq);
	close($SpalnFile);
	
	# Get that spaln output file THE HECK OUTTA HERE!
	RunSystemCommand("rm \"$spaln_fname\"") if (-e $spaln_fname);
	
	# If we didn't have a good return on investment, we're done-zo!
	if ($spaln_pct_id < $min_pct_id) {
	    push(@SpalnHitStrs,0);
	    push(@SpalnPctIDs,0);
	    next;
	}
	
	# All of our nucleotide coordinates need to be adjusted...
	my @UnadjNuclRanges = @{$nucl_ranges_ref};
	my @UnadjCenters = @{$centers_ref};

	my @HitChrs;
	my @ChrsUsed;
	my $num_chrs_used = 0;
	my @NuclRanges;
	my @AminoRanges = @{$amino_ranges_ref};
	my @CodonCenters;
	my $num_exons = 0;
	while ($num_exons < scalar(@UnadjNuclRanges)) {
	    
	    $UnadjNuclRanges[$num_exons] =~ /^(\d+)\.\.(\d+)$/;
	    my $range_start = $1;
	    my $range_end   = $2;
	    my $start_chr;
	    my $end_chr;
	    my $hit_chr;
	    
	    ($range_start,$start_chr) = AdjustNuclCoord($range_start,\@JumpList);
	    ($range_end,$end_chr) = AdjustNuclCoord($range_end,\@JumpList);
	    $NuclRanges[$num_exons] = $range_start.'..'.$range_end;
	    
	    # Clean up the chromosome -- NOW WITH CHIMERIC EXON CAPABILITIES!
	    if ($start_chr ne $end_chr) {
		
		if ($end_chr =~ /\-$/) { $end_chr =~ s/\-$/\[revcomp\]/; }
		else                   { $end_chr =~ s/\+$//;            }
		if ($start_chr =~ /\-$/) { $start_chr =~ s/\-$/\[revcomp\]/; }
		else                     { $start_chr =~ s/\+$//;            }
		$hit_chr = $start_chr.'/'.$end_chr;
		
	    } else {
		
		if ($end_chr =~ /\-$/) { $end_chr =~ s/\-$/\[revcomp\]/; }
		else                   { $end_chr =~ s/\+$//;            }
		$hit_chr = $end_chr;
		
	    }
	    push(@HitChrs,$hit_chr);

	    # Check if the series of chromosomes used in this hit needs updating
	    if (!$num_chrs_used || $ChrsUsed[$num_chrs_used-1] ne $hit_chr) {
		push(@ChrsUsed,$hit_chr);
		$num_chrs_used++;
	    }

	    # Acquire the actual codon centers
	    my $codon_center_str = '';
	    foreach my $coord (split(/\,/,$UnadjCenters[$num_exons])) {
		my ($codon_center,$codon_chr) = AdjustNuclCoord($coord,\@JumpList);
		$codon_center_str = $codon_center_str.','.$codon_center;
	    }
	    $codon_center_str =~ s/^\,//;
	    $CodonCenters[$num_exons] = $codon_center_str;
	
	    $num_exons++;
	    
	}
	
	my $chr;
	if ($num_chrs_used > 1) {
	    
	    $chr = 'Chimeric:';
	    foreach my $used_chr (@ChrsUsed) {
		$chr = $chr.$used_chr.'/';
	    }
	    $chr =~ s/\/$//;

	    # If we've set a threshold lower than 90%, we'll attach a note that this
	    # is a low quality mapping (if it's less than 90%).	    
	    $chr = $chr.' # low quality mapping: '.$spaln_pct_id.'% identity'
		if ($spaln_pct_id < 90.0);

	} else {
	    $chr = $ChrsUsed[0];
	}
	
	# Build up the hit string for this chromosome
	my $hitstr = "Sequence ID: $seqname\n";
	$hitstr    = $hitstr."Map Method : BLAT+Spaln\n";
	$hitstr    = $hitstr."Chromosome : $chr\n";
	$hitstr    = $hitstr."Num Exons  : $num_exons\n";
	for (my $i=0; $i<$num_exons; $i++) {
	    $hitstr= $hitstr."* Aminos $AminoRanges[$i], $HitChrs[$i]:$NuclRanges[$i]\n";
	    $hitstr= $hitstr."$CodonCenters[$i]\n";
	}
	
	# WOOF
	push(@SpalnHitStrs,$hitstr);
	push(@SpalnPctIDs,$spaln_pct_id);

    }

    # Get that nucleotide file OUTTA HERE!
    RunSystemCommand("rm \"$nucl_fname\"");
    
    # DOUBLE WOOF
    return (\@SpalnHitStrs,\@SpalnPctIDs);
    
}








############################################################
#
#  Function: ClusterNuclRanges
#
sub ClusterNuclRanges
{
    my $starts_ref = shift;
    my $ends_ref   = shift;
    my $max_span   = shift;

    my @Starts     = @{$starts_ref};
    my @Ends       = @{$ends_ref};
    my $num_ranges = scalar(@Starts);

    my $revcomp = 0;
    $revcomp = 1 if ($Starts[0] > $Ends[0]);

    # What are the minimum and maximum values?
    my $min;
    my $max;
    if ($revcomp) {
	$min = $Ends[0];
	$max = $Starts[0];
	for (my $i=1; $i<$num_ranges; $i++) {
	    $min = Min($min,$Ends[$i]);
	    $max = Max($max,$Starts[$i]);
	}
    } else {
	$min = $Starts[0];
	$max = $Ends[0];
	for (my $i=1; $i<$num_ranges; $i++) {
	    $min = Min($min,$Starts[$i]);
	    $max = Max($max,$Ends[$i]);
	}
    }

    # If we're already under the span, we're done!
    if ($max-$min < $max_span) {
	return $max.'..'.$min if ($revcomp);
	return $min.'..'.$max;
    }

    # Nope, recursion time!  Group sequences according to who they're closest to.
    # Note that this is a little fast 'n' loose, but given the span is as large as
    # it is, we aren't worried about (what I assume would be) minor roughness
    my @LeftStarts;
    my @LeftEnds;
    my @RightStarts;
    my @RightEnds;
    for (my $i=0; $i<$num_ranges; $i++) {
	my $center = ($Starts[$i] + $Ends[$i]) / 2;
	if ($max - $center > $center - $min) {
	    # Closer to min than max (left)
	    push(@LeftStarts,$Starts[$i]);
	    push(@LeftEnds,$Ends[$i]);
	} else {
	    # Closer to max than min (right)
	    push(@RightStarts,$Starts[$i]);
	    push(@RightEnds,$Ends[$i]);
	}
    }

    # RECURSE
    my $left_ranges_str  = ClusterNuclRanges(\@LeftStarts,\@LeftEnds,$max_span);
    my $right_ranges_str = ClusterNuclRanges(\@RightStarts,\@RightEnds,$max_span);

    # And away you go!
    return $right_ranges_str.','.$left_ranges_str if ($revcomp);
    return $left_ranges_str.','.$right_ranges_str;
    
}







############################################################
#
#  Function: AdjustNuclCoord
#
sub AdjustNuclCoord
{
    my $coord    = shift;
    my $jump_ref = shift;

    # We'll need to be a little sneaky to figure out which position is the last
    # *before* our guide's value is too high...
    my @Jumps = @{$jump_ref};
    my $index = scalar(@Jumps)-1;
    for (my $i=1; $i<scalar(@Jumps); $i++) {
	$Jumps[$i] =~ /^(\d+)\:/;
	my $relative_pos = $1;
	if ($relative_pos > $coord) {
	    $index = $i-1;
	    last;
	}
    }

    $Jumps[$index] =~ /^(\d+)\:(\d+)\|(\S+)$/;
    my $relative_pos = $1;
    my $absolute_pos = $2;
    my $chr = $3;

    my $revcomp = 0;
    $revcomp = 1 if ($chr =~ /\-$/);
    
    # We've found our closest jump point -- but we need to get particular about things!
    if ($revcomp) { $absolute_pos -= $coord - $relative_pos; }
    else          { $absolute_pos += $coord - $relative_pos; }

    return ($absolute_pos,$chr);
    
}







############################################################
#
#  Function: ParseBlatLine
#
sub ParseBlatLine
{
    my $line = shift;
    my @Data = split(/\s+/,$line);

    # NOTE: This assumes that we've pulled the sequence name off the front!
    my $chr         = $Data[0];
    # 1,2,3,4 aren't useful for me
    my $amino_start = $Data[5];
    my $amino_end   = $Data[6];
    my $nucl_start  = $Data[7];
    my $nucl_end    = $Data[8];
    my $score       = $Data[10];

    return($chr,$amino_start,$amino_end,$nucl_start,$nucl_end,$score);
    
}






############################################################
#
#  Function: ParseSpalnOutput
#
sub ParseSpalnOutput
{

    my $SpalnFile = shift;
    my $seq_ref   = shift;

    my @Seq = @{$seq_ref};

    # Alright, let's work our way through this inevitable disaster.

    # First off, let's get to the start of the alignment
    my $line;
    while ($line = <$SpalnFile>) {
	last if ($line =~ /^ALIGNMENT/);
    }
    return (0,0,0,0) if (eof($SpalnFile));

    # We'll read each region into a triple of sequences.
    # NOTE that these could contain multiple "exons," so I'm just calling them groups
    my @GroupTransSeqs;
    my @GroupNuclSeqs;
    my @GroupProtSeqs;
    my $group_trans_str;
    my $group_nucl_str;
    my $group_prot_str;
    my @GroupNuclStarts;
    my @GroupProtStarts;
    my $num_groups = 0;

    # We'll do some checking for inconsisntent skip distances, in
    my $nucl_pos    =  0; # Only advances by counting
    my $recent_skip = -1; # Did we just 'skip'?
    while ($line = <$SpalnFile>) {

	# This could be a skip... stay frosty!
	# The final line also replaces where the translated seq line would be,
	# so this is a good time to check if we're at the end of the road.
	$line = <$SpalnFile>;
	if ($line =~ /^\;\; skip/ || eof($SpalnFile)) {
	    push(@GroupTransSeqs,$group_trans_str);
	    push(@GroupNuclSeqs,$group_nucl_str);
	    push(@GroupProtSeqs,$group_prot_str);
	    $num_groups++;
	    last if (eof($SpalnFile));
	    $line =~ /^\;\; skip (\d+)/;
	    $recent_skip = $1;
	    next;
	}

	# Phew, no skip -- get ready to record some data!
	my $trans_line = $line;
	my $nucl_line  = <$SpalnFile>;
	last if (eof($SpalnFile));
	my $prot_line  = <$SpalnFile>;
	last if (eof($SpalnFile));

	# If we just skipped, we'll want to make sure that our position
	# makes sense.
	if ($recent_skip) {

	    # Does the nucleotide index agree with the skip distance?
	    # The only time we'll take its word over our skip counting
	    # during a disagreement is if the skip is 60 nt's, in which
	    # case we'll assume that the skip was unnecessary.
	    $nucl_line =~ /^\s*(\d+)/;
	    my $line_nucl_pos = $1;
	    $GroupNuclStarts[$num_groups] = $line_nucl_pos;

	    # A catch for the possibility of the line's reported position
	    # not matching our bookkeeping
	    if ($recent_skip != -1 && $nucl_pos + $recent_skip != $line_nucl_pos) {
		if ($recent_skip == 60 && $line_nucl_pos == $nucl_pos) {
		    $GroupNuclStarts[$num_groups] = $nucl_pos;
		} else {
		    $GroupNuclStarts[$num_groups] = $nucl_pos + $recent_skip;
		}
	    }
	    $nucl_pos = $GroupNuclStarts[$num_groups];

	    # We won't contest the protein position, just in case there's a legit
	    # skip
	    $prot_line =~ /^\s+(\d+)/;
	    $GroupProtStarts[$num_groups] = $1;

	    # Reset the group sequence strings
	    $group_trans_str = '';
	    $group_nucl_str  = '';
	    $group_prot_str  = '';

	    # Do you want to do all this special work all over again?!
	    # I sure don't!
	    $recent_skip = 0;

	}

	my @Trans = split(//,$trans_line);
	my @Nucl  = split(//,$nucl_line);
	my @Prot  = split(//,$prot_line);

	# Scan until we're at the nucleotide-y stuff
	my $scan = 0;
	while ($Nucl[$scan] =~ /\s/) { $scan++; }
	while ($Nucl[$scan] =~ /\d/) { $scan++; }
	while ($Nucl[$scan] =~ /\s/) { $scan++; }
	
	while ($Nucl[$scan] ne '|') {
	    if ($Nucl[$scan] !~ /\s/) {
		$group_trans_str = $group_trans_str.$Trans[$scan];
		$group_nucl_str  = $group_nucl_str.$Nucl[$scan];
		$group_prot_str  = $group_prot_str.$Prot[$scan];
		$nucl_pos++ if ($Nucl[$scan] =~ /[A-Za-z]/);
	    }
	    $scan++;
	}

    }

    # Now that we've pulled in all of the Spaln output, we'll run through
    # each of the segments and extract exons (in the style of our actual
    # output format)
    my @ExonNuclRanges;
    my @ExonAminoRanges;
    my @CodonCenters;
    my $num_matches    = 0;
    my $num_mismatches = 0;
    for (my $i=0; $i<$num_groups; $i++) {

	my @Trans = split(//,$GroupTransSeqs[$i]);
	my @Nucl  = split(//,$GroupNuclSeqs[$i]);
	my @Prot  = split(//,$GroupProtSeqs[$i]);

	my $nucl_pos = $GroupNuclStarts[$i];
	my $prot_pos = $GroupProtStarts[$i];

	# Scan through the nucleotide sequence, looking for points where
	# we switch between upper and lower case characters (exons / introns)
	my $scan = 0;	
	while ($scan < scalar(@Nucl)) {

	    while ($scan < scalar(@Nucl) && $Nucl[$scan] eq lc($Nucl[$scan])) {
		$nucl_pos++ if ($Nucl[$scan] =~ /[a-z]/);
		$scan++;
	    }

	    last if ($scan == scalar(@Nucl));

	    my $exon_nucl_start  = $nucl_pos;
	    my $exon_prot_start  = $prot_pos;
	    my $codon_center_str = '';

	    my $start_prot_pos = $prot_pos;
	    while ($scan < scalar(@Nucl) && $Nucl[$scan] eq uc($Nucl[$scan])) {
		if ($Prot[$scan] =~ /[A-Z]/) {
		    # Do we have a match?
		    if ($Prot[$scan] eq $Trans[$scan] || ($Prot[$scan] eq 'S' && $Trans[$scan] eq 'J')) {
			$num_matches++;
		    } else {
			$num_mismatches++;
		    }
		    $codon_center_str = $codon_center_str.','.$nucl_pos;
		    $prot_pos++;
		}
		$nucl_pos++ if ($Nucl[$scan] =~ /[A-Z]/);
		$scan++;
	    }

	    next if ($start_prot_pos == $prot_pos);

	    my $exon_nucl_end = $nucl_pos-1;
	    my $exon_prot_end = $prot_pos-1;

	    # We'll end up with a highly undesirable leading comma, bringing shame upon
	    # our family name if left uncontested
	    $codon_center_str =~ s/^\,//;

	    # Turn our positions into ranges and push 'em onto the list
	    my $exon_nucl_range = $exon_nucl_start.'..'.$exon_nucl_end;
	    my $exon_prot_range = $exon_prot_start.'..'.$exon_prot_end;
	    push(@ExonNuclRanges,$exon_nucl_range);
	    push(@ExonAminoRanges,$exon_prot_range);
	    push(@CodonCenters,$codon_center_str);

	}

    }

    # If we don't have *any* matches, something's gone very wrong
    return(0,0,0,0) if ($num_matches == 0);

    # Something's also gone wrong if we don't have full coverage
    return(0,0,0,0) if ($num_matches+$num_mismatches != scalar(@Seq));

    # For curious minds, what was our percent identity?
    my $pct_id = int(1000.0 * $num_matches / ($num_matches + $num_mismatches)) / 10.0;

    # I think that's all we need to do for now... Hot dog!
    return(\@ExonNuclRanges,\@ExonAminoRanges,\@CodonCenters,$pct_id);

}






############################################################
#
#  Function: FinalFileCheck
#
sub FinalFileCheck
{
    my $dirname  = shift;
    my $num_cpus = shift;

    # This shouldn't be necessary, but we'll just make sure everything looks reasonable
    $dirname = ConfirmDirectory($dirname);
    
    # Make a list of all of the '.quilter.out' files
    my @OutFileList;
    my $Dir = OpenDirectory($dirname);
    while (my $fname = readdir($Dir)) {
	push(@OutFileList,$dirname.$fname) if ($fname =~ /\.quilter\.out$/);
    }
    closedir($Dir);

    # If we don't have any result files, that's some real rough stuff
    return if (scalar(@OutFileList) == 0);

    # We'll generate a bunch of processes to make things go down sooooo smooth
    $num_cpus = Min($num_cpus,scalar(@OutFileList));
    my $threadID = SpawnProcesses($num_cpus);

    # Off to the file mines with you!
    my $start_file_index =  $threadID    * int(scalar(@OutFileList)/$num_cpus);
    my $end_file_index   = ($threadID+1) * int(scalar(@OutFileList)/$num_cpus);
    $end_file_index = scalar(@OutFileList) if ($threadID == $num_cpus-1);

    my $num_cleaned_files = 0;
    for (my $i=$start_file_index; $i<$end_file_index; $i++) {

	my $fname = $OutFileList[$i];

	# We'll hash each of the sequences according to the seq name (which is
	# assumed to still be a numeric code at this point)
	my %HitsBySeq;
	my %Chrs;
	my $inf = OpenInputFile($fname);
	while (my $line = <$inf>) {

	    if ($line =~ /Sequence ID\: (\S+)/) {

		my $seqname = $1;
		my $fullhit = $line;

		$line = <$inf>; # Map Method
		$fullhit = $fullhit.$line;

		$line = <$inf>; # Chromosome
		$fullhit = $fullhit.$line;
		$line =~ /Chromosome \: (\S+)/;
		my $chr = $1;
		if ($Chrs{$chr}) { $Chrs{$chr}++; }
		else             { $Chrs{$chr}=1; }

		$line = <$inf>; # Number of Exons
		$fullhit = $fullhit.$line;
		$line =~ /Num Exons  \: (\d+)/;
		my $num_exons = $1;

		while ($num_exons) {
		    $line = <$inf>; # Exon overview
		    $fullhit = $fullhit.$line;
		    $line = <$inf>; # Codon centers
		    $fullhit = $fullhit.$line;
		    $num_exons--;
		}

		# Store it!
		$HitsBySeq{$seqname} = $fullhit;
	    }
	    
	}

	close($inf);

	# We now have all of the hits that we were able to pull -- let's see
	# which chromosome appears to be canonical.
	# We'll ignore 'Incomplete' and 'Chimeric' mappings, for obvious reasons
	my $top_chr = 'Undetermined';
	my $top_chr_count = 0;
	foreach my $chr (keys %Chrs) {
	    if ($chr !~ /Incomplete|Chimeric/ && $Chrs{$chr} > $top_chr_count) {
		$top_chr = $chr;
		$top_chr_count = $Chrs{$chr};
	    }
	}

	# We'll also pull in the full list of sequence names for this family
	# so we can record which sequences weren't mapped.
	#
	# We'll note sequences that mapped to noncanonical chromosomes or
	# mapped chimerically as unmapped, though with a distinction.
	$fname =~ s/\.quilter\.out/\.fa/;
	$inf = OpenInputFile($fname);
	while (my $line = <$inf>) {
	    if ($line =~ /\>(\S+)/) {
		my $seqname = $1;
		if (!$HitsBySeq{$seqname}) {
		    my $miss_str = "Sequence ID: $seqname\nMap Method : Unmapped\n";
		    $HitsBySeq{$seqname} = $miss_str;
		}
	    }
	}
	close($inf);

	# Now we'll write out our full output along with all of the 'miss' information
	# to the final files for these families.
	# We'll use the ugly opening code for the hitfile because we're overwriting.
	my $hitfname = $fname;
	$hitfname =~ s/\.fa$/\.quilter\.out/;
	open(my $hitf,'>',$hitfname) || die "\n  ERROR: Failed to open final output file '$hitfname'\n\n";

	# What's the canonical chromosome for this family?
	print $hitf "Canonical Chromosome: $top_chr\n\n";

	my $missfname = $hitfname;
	$missfname =~ s/\.quilter\./\.misses\./;
	my $missf = OpenOutputFile($missfname);

	foreach my $seqname (sort { $a <=> $b } keys %HitsBySeq) {

	    # Hit acquired! (or not... but why ruin the moment like that, dude?)
	    print $hitf "$HitsBySeq{$seqname}\n";

	    # If there's any sense by which this sequence is unmapped, record it
	    if ($HitsBySeq{$seqname} =~ /Map Method \: Unmapped/) {
		print $missf "Unmapped   : $seqname\n";
	    } else {
		$HitsBySeq{$seqname} =~ /Chromosome \: (\S+)/;
		my $chr = $1;
		if ($chr =~ /Incomplete|Chimeric/) {
		    print $missf "$chr: $seqname\n";
		}
	    }

	}
	
	close($hitf);
	close($missf);

	# If we're hashtag blessed enough to not have had any misses for this gene,
	# no need for a miss file!
	if (!(-s $missfname)) {
	    RunSystemCommand("rm \"$missfname\"") if (-e $missfname);
	}

	$num_cleaned_files++;
	DispProgMirage('cleanup|'.$threadID.'|'.$num_cleaned_files);
	
    }
    
    # Take a breather -- you've earned it!
    if ($threadID) { exit(0); }
    while (wait() != -1) {}
    
}





#
#  NOTE: The following functions are taken from the original Quilter,
#        which underwent a long period of incremental fixes to account
#        for errors in Spaln output.
#
#        The 'ParseSpalnOutput' function above is more along the lines of
#        what I'd like (and will probably be transitioned into a parsing
#        function for an eventual Spaln replacement), but for the sake of
#        not having to translate all of the error catching stuff (at least
#        not simultaneously) I'm going to stick with this uglier code.
#
#        OPSO_ := 'Original Parse Spaln Output'
#




#########################################################################
#
#  Function Name: ParseSPALNOutput
#
#  About: The name says it all.  Run a SPALN command and make sense of
#         its wisdom.
#
#    1.  Don't trust percent identities.
#    2.  J = S.
#
#  >> FOR ALEX:
#
#     New expected input is:
#
#          1. A file with the Spaln output
#          2. The original sequence, as an array
#
#     New expected output is:
#
#          1. A list of nucleotide coordinate ranges formatted as <start>..<end>
#          2. A list of amino acid coordinate ranges formatted as <start>..<end>
#          3. A list of strings with comma-seperated amino acid center coordinates
#          4. The percent identity of the mapping
#
sub OPSO_ParseSpalnOutput
{
    my ($i,$j,$k);

    my $SpalnFile = shift;
    my $seq_ref   = shift;

    my @Seq = @{$seq_ref};
    my $prot_len = scalar(@Seq);

    # Are we timing?
    my $timing     = shift;
    my $timingdata = shift;
    my $timeA;

    my ($line,$hitstring);

    # NOW we can go to the real business (the alignment lines)
    $line = readline($SpalnFile); # 'ALIGNMENT'
    while (!eof($SpalnFile) && $line !~ /ALIGNMENT/) {
	$line = readline($SpalnFile);
    }

    # Still not totally stoked on seeing an eof
    return(0,0,0,0) if (eof($SpalnFile));
    
    #
    # ALRIGHT, GENTS AND LADIES, it looks like we're actually
    # doing this thing.
    #

    # Scan through the SPALN output constructing 4 arrays:
    #
    #     - An array with the DNA sequence characters
    #     - An array with the protein sequence characters
    #     - An array mapping protein sequence indices to
    #       DNA sequence indices
    #     - An array mapping DNA sequence characters to
    #       positions in the genome
    #
    # The combination of these arrays will allow us to do
    # all of the work we need in terms of generating our
    # final mapping of the protein to the genome, correcting
    # for micro-exons.
    #
    # UPDATE (Nov. 2019): It looks like Spaln will, somewhat
    # regularly, start/end exons with extraneous gaps, so I'm
    # going to try to pay attention to the starts and ends of
    # exons to make sure these don't drag good hits into the
    # mud.
    #
    my $gap_run_len = 0;
    my $starting_exon = 1;
    my @FullNuclSeq;
    my @FullProtSeq;
    my @AAPositions;
    my @NuclPositions;
    my $trans_line;
    my $full_length = 0;
    my $num_aas     = 0; 
    my $current_pos = 0;
    my $last_end_pos; # To catch a possible SPALN error (false 'skip's)
    my $first_jump  = 0;
    while (!eof($SpalnFile)) {

	# Grab the next line, clean it up, make sure it's
	# meaningful.
	#
	$line = readline($SpalnFile);
	$line =~ s/\n|\r//g;
	next if (!$line);


	# If we're skipping around in the genome, adjust
	# the position variable accordingly.
	#
	if ($line =~ /skip\s+(\d+)\s+nt/) {
	    $current_pos += int($1);
	    $starting_exon = 1; # Whatever we see next will start an exon
	    next;
	}


	# If we've made it this far but don't match this
	# terminal format then we're looking at SPALN's
	# translation of the nucleotide sequence
	#
	if ($line !~ /\| \S+$/) {
	    $trans_line = $line;
	    next;
	}

	my $nn_line = $line;

	# EXCELLENT!  This is the next line of DNA characters
	#
	my @NextNucls = split(//,$line);

	# We also split up what would be the translation, so that
	# we can evaluate percent ID for ourselves.
	#
	my @TransLine = split(//,$trans_line);

	# The final check we need to do is to add/subtract
	# the initial 'jump' into the hit (but only if this
	# is the very beginning)
	#
	if ($first_jump == 0) {
	    if ($line =~ /^\s*(\d+)/) {
		$current_pos += $1;
		$first_jump   = 1;
	    } else {
		return(0,0,0,0);
	    }
	}


	# Does it look like SPALN lied to us about skipping?
	$nn_line =~ /^\s*(\d+)\s/;
	my $reported_pos = $1;
	if ($last_end_pos && $last_end_pos == $reported_pos) {
	    $current_pos = $last_end_pos;
	}

	
	# ... which means that the following line is the next
	# line of protein characters.
	#
	$line = readline($SpalnFile);
	$line =~ s/\n|\r//g;
	my @NextAAs = split(//,$line);


	# Advance to the actual DNA sequence.
	# NOTE:  it's possible for SPALN to have a line of
	#        '-' characters, which we'd want to throw
	#        away.  BUT we still need to count up the
	#        position indices when we encounter these. << Do we?
	#
	$i = 0;
	#while ($i < @NextNucls && $NextNucls[$i] !~ /[a-zA-Z]|\|/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] =~ /\s/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] !~ /\s/) { $i++; }
	while ($i < @NextNucls && $NextNucls[$i] =~ /\s/) { $i++; }
	


	# Some sort of SPALN weirdness...
	#
	next if ($i == @NextNucls || $NextNucls[$i] eq '|');


	# Run along the two lines storing all of the relevant
	# information that we can get out of them.
	#
	while ($NextNucls[$i] !~ /\|/) {
	    
	    push(@FullNuclSeq,$NextNucls[$i]);
	    push(@NuclPositions,$current_pos);

	    # If we're centered on a codon, take note.
	    #
	    if ($NextAAs[$i] =~ /\w/) {		
		push(@FullProtSeq,$NextAAs[$i]);
		push(@AAPositions,$full_length);
		$num_aas++;
	    }
	    
	    # Also, note a *SUBSTITUTION* mismatch (sorry for yelling)
	    if ($NextAAs[$i] =~ /\S/
		&& $NextAAs[$i] ne $TransLine[$i] 
		&& !(uc($TransLine[$i]) eq 'J' && uc($NextAAs[$i]) eq 'S')) { # Sometimes SPALN calls an 'S' a 'J'

		if ($NextAAs[$i] eq '-' && !$starting_exon) {
		    $gap_run_len++ if (!$starting_exon);
		}
		
	    } else {
		$gap_run_len=0;
	    }
	    # I'm using that format because we might see more delightful tricks


	    # Increment the position in the genome
	    if ($NextNucls[$i] =~ /\w/) {
		$current_pos++;
	    }

	    $full_length++;
	    $i++;
	    
	}

	$last_end_pos = $current_pos;
	
    }
    

    # Selenocysteine does some crazy stuff, man.  I translated some once
    # and I'm still coming down.
    if ($num_aas != $prot_len) {
	
	my ($fps_ref,$aap_ref) =
	    OPSO_SelenocysteineCheck(\@FullNuclSeq,\@FullProtSeq,\@AAPositions,\@Seq);
	@FullProtSeq = @{$fps_ref};
	@AAPositions = @{$aap_ref};
	$num_aas = scalar(@FullProtSeq);

	return (0,0,0,0) if ($num_aas != $prot_len);	

    }

    # Alrighty, now we oughta have the most complete story of the translated sequence
    # to determine the percent identity.
    my $num_matches = 0;
    for (my $i=0; $i<$prot_len; $i++) {
	if (uc($FullProtSeq[$i]) eq uc($Seq[$i])) {
	    $num_matches++;
	}
    }
    my $spaln_pct_id = int(1000.0 * $num_matches / $prot_len) / 10.0;


    # ********************************************************************
    # *** 'IF ZERO' WILL TOGGLE MICRO EXON SEARCHING OFF *****************
    # ********************************************************************
    #
    #
    #if (0) {

    # Run along SPALN's output looking for 'exons' with
    # improbably low amino acid counts.  Currently, we
    # say that an exon should have at least 4 amino acids.
    #
    my @MicroExonRuns;
    my $micro_start = -1;
    my $micro_end   = 0;
    
    $i = 0;
    while ($i < $num_aas) {


	# Figure out the genome position of the current amino
	# acid, so that we can infer exon boundaries (large
	# changes in nucleotide position).
	#
	my $aa_loc   = $NuclPositions[$AAPositions[$i]];
	my $aa_start = $i;
	my $exon_len = 1;

	
	# Run through this exon, counting how many amino acids
	# it has.
	#
	$i++;
	while ($i < $num_aas && abs($NuclPositions[$AAPositions[$i]] - $NuclPositions[$AAPositions[$i-1]]) < 5) {
	    $exon_len++;
	    $i++;
	}


	# If we have too few amino acids, then we're working
	# with a micro exon, baby!  Observe that we chain together
	# runs of contiguous micro-exons.
	#
	if ($exon_len < 4) {
	    
	    if ($micro_start == -1) { $micro_start = $aa_start; }
	    $micro_end = $i;
	    
	} elsif ($micro_start != -1) {
	    
	    push(@MicroExonRuns,$micro_start.'-'.$micro_end);
	    $micro_start = -1;
	    
	}
	
    }
    
    
    # Next we'll run through each of the indicated micro-exon
    # regions ('regions' because they can be chained together)
    # and try to re-align the amino acids to the ends of the
    # 'major' exons that flank them.
    #
    # NOTE: Be careful that we don't step over a breakpoint by
    #       checking whether NuclPositions[i] and NP[i+1] are
    #       1 off from one another.
    #
    foreach my $micro_run_str (@MicroExonRuns) {

	
	# Split the indices out of the string
	#
	my @MicroRun = split(/\-/,$micro_run_str);
	$micro_start = $MicroRun[0];
	$micro_end   = $MicroRun[1]; # One past actual end, so we can do strict '<'


	# Identify the amino acid sequence that constitutes
	# this micro exon (or chain of micro exons)
	#
	my @MicroAAs;
	for ($j = $micro_start; $j < $micro_end; $j++) {
	    push(@MicroAAs,$FullProtSeq[$j]);
	}

	# Figure out how long the micro exon is (keeping in mind that '_end' is
	# one position past the actual end index)
	#
	my $micro_len = $micro_end-$micro_start;

	
	# Look at the preceding exon
	#
	my @RearExt;
	my $rear_ext_len = 0;
	if ($micro_start > 0) {

	    
	    # Get situated so that we can start trying to
	    # transfer amino acids from our micro exons to
	    # the back end of the preceding exon.
	    #
	    # NOTE (for making sense of this):
	    # --------------------------------
	    # NuclPositions are the chromosome indices, AAPositions are the
	    # array indices that SPALN has identified as codon centers -- thus
	    # looking up the NuclPosition of an AAPosition gets you the genome
	    # index for the AAPosition's codon center.
	    #
	    my $rear_codon_center = $NuclPositions[$AAPositions[$micro_start-1]];
	    my $nuclseq_pos = $AAPositions[$micro_start-1]+3; # Which array index are we centered around?
	    while ($rear_ext_len < $micro_len) {
		
		# Make sure we don't overstep (by checking if we've passed over
		# an stretch of genome indices indicating an intron)
		my $step_check = $NuclPositions[$nuclseq_pos];
		$step_check -= 3 * ($rear_ext_len+1);
		
		last if ($step_check != $rear_codon_center);


		# Figure out what the next codon is
		my @NextCodon = ();
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos-1]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos+1]);

		
		# If we have a match we are very happy.  If
		# we do not have a match we are very sad.
		#
		# Why would we do 'micro_start+(rear_ext_len+1)'?
		#
		if (uc($FullProtSeq[$micro_start+$rear_ext_len]) eq uc(TranslateCodon(\@NextCodon))) {
		    push(@RearExt,$nuclseq_pos);
		    $nuclseq_pos += 3;
		    $rear_ext_len++;
		} else {
		    last;
		}
		
	    }
	    
	}
	
	# If we found a complete run, we're sooooo stoked!
	if ($rear_ext_len == $micro_len) {
	    for ($i=0; $i<$micro_len; $i++) {
		$AAPositions[$micro_start+$i] = $RearExt[$i];
	    }
	    next;
	}

	
	# Look at the exon ahead
	my @FwdExt;
	my $fwd_ext_len = 0;
	if ($micro_end < $num_aas) {
	    
	    # Try to find a good extension to the front of
	    # the upcoming exon.
	    my $fwd_codon_center = $NuclPositions[$AAPositions[$micro_end]];
	    my $nuclseq_pos = $AAPositions[$micro_end]-3;
	    while ($fwd_ext_len < $micro_len) {

		
		# Make sure we don't overstep
		my $step_check = $NuclPositions[$nuclseq_pos];
		$step_check += 3 * ($fwd_ext_len+1);
		
		last if ($step_check != $fwd_codon_center);


		# What does the next codon encode?
		my @NextCodon = ();
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos-1]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos]);
		push(@NextCodon,$FullNuclSeq[$nuclseq_pos+1]);
		

		# If we have a match, in terms of what the codon
		# encodes, then we can continue our considerations.
		# Otherwise we cannot.
		#
		# Note that the first nucleotide position in the
		# 'FwdExt' array corresponds to the last amino acid
		# in the micro-exon region.
		#
		if (uc($FullProtSeq[$micro_end-($fwd_ext_len+1)]) eq uc(TranslateCodon(\@NextCodon))) {
		    push(@FwdExt,$nuclseq_pos);
		    $nuclseq_pos -= 3;
		    $fwd_ext_len++;
		} else {
		    last;
		}
		
	    }
	    
	}

	
	# Did we find a complete run?
	if ($fwd_ext_len == $micro_len) {
	    for ($i=0; $i<$micro_len; $i++) {
		$AAPositions[$micro_end-($i+1)] = $FwdExt[$i];
	    }
	    next;
	}
	

	# Append as much as you can to the preceding
	# major exon first.  We don't put any work into
	# evaluating whether amino acids that could go
	# either way should be fitted to one or the other
	# side.
	#
	$i = 0;
	while ($i < $rear_ext_len) {
	    $AAPositions[$micro_start+$i] = $RearExt[$i];
	    $i++;
	}

	
	# Skip ahead to the stuff we're appending to the
	# upcoming major exon, if you need to.
	if ($i < $micro_len - $fwd_ext_len) {
	    $i = $micro_len - $fwd_ext_len;
	}


	# Append this stuff onto the forward friend.
	#
	$j = $fwd_ext_len - 1;
	while ($i < $micro_len) {
	    $AAPositions[$micro_end-($micro_len-$i)] = $FwdExt[$j];
	    $i++;
	    $j--;
	}
	
    }
    

    #}
    #
    #
    # ********************************************************************
    # *** 'IF ZERO' TOGGLE END *******************************************
    # ********************************************************************

        
    # Radical, dude!  Now all we need to do is just convert
    # our list of amino acid positions into a list of
    # their corresponding middle amino acids (with '*'s to
    # denote splice site boundaries).
    #
    # We'll try to be careful about how we place codons that are
    # called as insertions into the genome.
    #
    my $prev_codon_c = $NuclPositions[$AAPositions[0]];
    my $codon_center_str = $prev_codon_c;
    my $amino_range_str  = '1';
    my @Centers;
    my @AminoRanges;
    my @NuclRanges;
    for ($i=1; $i<$num_aas; $i++) {

	# A gap of more than 4 is taken to deserve a splice site marker
	#
	if ($NuclPositions[$AAPositions[$i]] - $NuclPositions[$AAPositions[$i-1]] > 4) {
	    
	    $codon_center_str =~ s/\,$//;
	    push(@Centers,$codon_center_str);
	    $codon_center_str = '';

	    $amino_range_str = $amino_range_str.'..'.$i;
	    push(@AminoRanges,$amino_range_str);

	    # While we have the two amino end points, let's get some nucleotide ranges
	    # like we're gall-darn maniacs.
	    #
	    # NOTE: I'm not sure how this works in cases of micro-exons, since there
	    #       may not be perfect capitalization... At least it'll give a sense...
	    #
	    $amino_range_str =~ /(\d+)\.\.(\d+)/;
	    my $range_start = $AAPositions[$1-1];
	    my $range_end   = $AAPositions[$2-1];

	    while ($range_start && $FullNuclSeq[$range_start-1] eq uc($FullNuclSeq[$range_start-1])) {
		$range_start--;
	    }
	    $range_start = $NuclPositions[$range_start];

	    while ($range_end < scalar(@FullNuclSeq)-1 && $FullNuclSeq[$range_end+1] eq uc($FullNuclSeq[$range_end+1])) {
		$range_end++;
	    }
	    $range_end = $NuclPositions[$range_end];

	    push(@NuclRanges,$range_start.'..'.$range_end);

	    $amino_range_str = $i+1;
	    
	}

	# Are we in an insertion (wrt genome)?
	#
	my $next_center;
	if ($FullNuclSeq[$AAPositions[$i]] !~ /\w/) {
	    $next_center = $prev_codon_c;
	} else {
	    $prev_codon_c = $NuclPositions[$AAPositions[$i]];
	    $next_center  = $prev_codon_c;
	}

	# Are we starting a new run?
	if ($codon_center_str) {
	    $codon_center_str = $codon_center_str.','.$next_center;
	} else {
	    $codon_center_str = $next_center;
	}
	
    }

    push(@Centers,$codon_center_str);

    $amino_range_str = $amino_range_str.'..'.$num_aas;
    push(@AminoRanges,$amino_range_str);

    # UGH, sorry about the copy-paste, but daddy lazy today :p
    $amino_range_str =~ /(\d+)\.\.(\d+)/;
    my $range_start = $AAPositions[$1-1];
    my $range_end   = $AAPositions[$2-1];
    
    while ($range_start && $FullNuclSeq[$range_start-1] eq uc($FullNuclSeq[$range_start-1])) {
	$range_start--;
    }
    $range_start = $NuclPositions[$range_start];
    
    while ($range_end < scalar(@FullNuclSeq)-1 && $FullNuclSeq[$range_end+1] eq uc($FullNuclSeq[$range_end+1])) {
	$range_end++;
    }
    $range_end = $NuclPositions[$range_end];
    
    push(@NuclRanges,$range_start.'..'.$range_end);
    
    # OLD PARSE SPALN OUTPUT COMPLETE!
    return(\@NuclRanges,\@AminoRanges,\@Centers,$spaln_pct_id);

}





#########################################################################
#
#  Function Name: SelenocysteineCheck
#
#  About: This function checks if a selonocyteine (U, encoded by TGA)
#         explains an observed absence of aminos in SPALN output.
#
sub OPSO_SelenocysteineCheck
{
    my $nucl_seq_ref = shift;
    my $prot_seq_ref = shift;
    my $aa_pos_ref   = shift;
    my $orig_seq_ref = shift;

    my @NuclSeq = @{$nucl_seq_ref};
    my @ProtSeq = @{$prot_seq_ref};
    my @AAPos   = @{$aa_pos_ref};
    my @OrigSeq = @{$orig_seq_ref};

    my $aa_pos = 0;
    while ($aa_pos < scalar(@OrigSeq)) {

	# Looks like we've discovered a selenocysteine in our midst!
	if ($OrigSeq[$aa_pos] eq 'U' && $aa_pos < scalar(@ProtSeq) && $ProtSeq[$aa_pos] ne 'U') {

	    # Is this consistent with stepping forward 3 from the previous
	    # amino's position?  3 back from forthcoming amino?
	    my $nailed_it = 0;
	    if ($aa_pos) {
		
		# Should be to the right of the last position...
		my $nucl_pos = $AAPos[$aa_pos-1] + 2;
		my $codon = uc($NuclSeq[$nucl_pos].$NuclSeq[$nucl_pos+1].$NuclSeq[$nucl_pos+2]);
		if ($codon eq 'TGA') {

		    # Note:  Splice moves everything at the indicated position to
		    #        the right and sticks in what you want there.
		    splice(@ProtSeq,$aa_pos,0,'U');
		    splice(@AAPos,$aa_pos,0,$nucl_pos+1);
		    $nailed_it = 1;
		} 

	    } 

	    if (!$nailed_it && $aa_pos < scalar(@ProtSeq)-1) {
		
		# ... Or left of the next position
		my $nucl_pos = $AAPos[$aa_pos] - 4;
		my $codon = uc($NuclSeq[$nucl_pos].$NuclSeq[$nucl_pos+1].$NuclSeq[$nucl_pos+2]);
		if ($codon eq 'TGA') {
		    splice(@ProtSeq,$aa_pos-1,0,'U');
		    splice(@AAPos,$aa_pos-1,0,$nucl_pos+1);
		} 

	    }

	}

	$aa_pos++;

    }

    return (\@ProtSeq,\@AAPos);

}











# EOF











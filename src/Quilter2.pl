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
sub GetBM { my $lib = $0; $lib =~ s/Quilter2.pl$//; return $lib; }
use lib GetBM();
use BureaucracyMirage;

sub PrintUsage;
sub ParseArgs;
sub ParseChromSizes;
sub UseGTF;
sub ParseGTF;
sub UseFastMap;
sub ParseFastMapOutput;
sub GenSpliceGraph;
sub CheckForFullMaps;
sub ExtractCanonHWMap;
sub AttemptChimericHWMap;
sub ExtractHWInputExons;
sub ExtractChimericHWMap;
sub FindMaximalHWHitSet;
sub RecordMaximalHits;
sub RunBlatOnFileSet;

if (@ARGV < 3) { PrintUsage(0); }

# As you'd expect, we'll parse the commandline arguments right up top
my $opts_ref = ParseArgs();
my %Opts     = %{$opts_ref};

# As you'd expect, we want to know what there is to know about our proteins
my $species_dirname = ConfirmDirectory($ARGV[0]);
my $seq_dirname = ConfirmDirectory($species_dirname.'seqs');
my $ali_dirname = CreateDirectory($species_dirname.'alignments');

# As you'd expect, next up is confirming the genome
my $genome = ConfirmFile($ARGV[1]);

# As you'd expect, we'll want to parse the gtf, if one's been provided
my $gtfname = $ARGV[2];

# As long as everything's going as expected, figure out where we are
my $srcdir = $0;
$srcdir =~ s/Quilter2.pl$//;
my $sfetch = $srcdir.'../inc/easel/miniapps/esl-sfetch';

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

# If we're skipping straight to BLAT (no gtf) then... let's!
# Otherwise, do the first round of real tough-guy work!
if ($gtfname ne '-') {

    UseGTF($seq_dirname,$ali_dirname,$genome,$gtfname,$num_cpus,\%Opts);

} else {

    # We'll need to put all of the files that we want to use in our
    # BLAT search into a list, because that's what I've decided we
    # need to do.
    my @BlatFileNames;
    my $SeqDir = OpenDirectory($seq_dirname);
    while (my $fname = readdir($SeqDir)) {
	next if ($fname !~ /\.fa$/);
	push(@BlatFileNames,$seq_dirname.$fname);
    }
    closedir($SeqDir);

    # We're running our BLAT search using file names as the way to derive
    # the gene families that sequences belong to, and we're doing it now!
    my ($blat_outfname,$blat_nameguide_ref,$num_blat_genes)
	= RunBlatOnFileSet('filename',@BlatFileNames,$genome);
    my %BlatNameGuide = %{$blat_nameguide_ref};

}

1;





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
    my $ali_dirname = shift;
    my $genome      = shift;
    my $gtfname     = shift;
    my $num_cpus    = shift;
    my $opts_ref    = shift;
    my %Opts = %{$opts_ref};

    # What was the name of our species, again?
    $seq_dirname =~ /\/([^\/]+)\/seqs\/$/;
    my $species = $1;

    # UHHHHHH, how exactly did you plan on using a GTF without PARSING IT FIRST?!
    my $gtf_ref = ParseGTF($gtfname,$species);
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

    # TIMING TEST!
    my $runtime = 0.0;
    
    # Work it, thready!
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
	
    }

    # Close up the blat file, and (if we're so blessed) clear it if it's empty
    close($BlatFile);
    if (!(-s $blatfname)) {
	RunSystemCommand("rm \"$blatfname\"") if (-e $blatfname);
    }
    
    # Wait for everyone to finish
    if ($threadID) { exit(0); }
    while (wait() != -1) {}

    die;

    # Great stuff, team!  Now, how's about we pull all the blat files together
    # and get real blatty up in this?
    my @BlatFileNames;
    for (my $i=0; $i<$num_cpus; $i++) {

	# It isn't guaranteed that this file will exist, in which case we
	# have nothing to do
	$blatfname = $seq_dirname.$threadID.'.blat.fa';
	next if (!(-e $blatfname));

	# Add it to the pile!
	push(@BlatFileNames,$blatfname);

    }

    # It's the dlark knlight! BLATMAN!
    my ($blat_outfname,$blat_nameguide_ref,$num_blat_genes)
	= RunBlatOnFileSet('seqname',\@BlatFileNames,$genome);
    my %BlatNameGuide = %{$blat_nameguide_ref};
    
}





############################################################
#
#  Function: ParseGTF
#
sub ParseGTF
{
    my $gtfname = shift;
    my $species = shift;

    my $GTFile = OpenInputFile($gtfname);

    # Essentially, what we want to do is pull all of the GTF entries for our
    # species into a hash
    my %GTF;
    my %DuplicateCheck;
    while (my $line = <$GTFile>) {

        next if ($line =~ /^\#/ || lc($line) !~ /\S+\s+\S+\s+exon/);

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

    # If there isn't an entry for this gene, fuhget about it!
    return if (!$GTF{$gene});

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
    # HitWeaver outputs.
    my @SpliceGraphs;
    my @HWInputNames;
    for (my $i=0; $i<$num_seqs; $i++) {
	$HWInputNames[$i] = 0;
	$SpliceGraphs[$i] = 0;
    }

    # Now we can go chromosome-by-chromosome and try to build the
    # best splice graph for our sequences we can on a single chromosome
    foreach my $chr (keys %RangesByChr) {

	my @Ranges = split(/\//,$RangesByChr{$chr});

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
	    }

	    # The cutest buffer you've ever seen
	    $search_start += 20;
	    $search_end   -= 20;
	    
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
	    }

	    # A teensy tiny little buffer
	    $search_start -= 20;
	    $search_end   += 20;
	    
	}

	# Extract the appropriate genomic region for our search
	RunSystemCommand($sfetch." -c $search_start\.\.$search_end \"$genome\" \"$chr\" > \"$nucl_fname\"");

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
	#   6. nucleotide sequence (w/ 17 lowercase buffer nucls on each side)
	#
	# Moreover, they've been run through Perl's sort, so we should be in
	# order of starting position in the amino acid sequence.
	#
	# Our next task is to ID the best splice graph that we can for each of
	# the sequences, individually.

	my @HitsBySeq = @{$hitref};
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($HitsBySeq[$i]) {

		# Get the name of the files provided to and produced by HitWeaver
		my ($hw_infname,$hw_outfname) =
		    GenSpliceGraph($HitsBySeq[$i],$Seqs[$i],$SeqNames[$i],$chr,$revcomp,$gene_fname);

		# If we didn't get any spliced alignments out of ol' HW, move along
		next if (!$hw_outfname);

		# That's lookin' like a graph to me!  Toss it on the pile!
		if ($SpliceGraphs[$i]) {
		    $HWInputNames[$i] = $HWInputNames[$i].';'.$hw_infname;
		    $SpliceGraphs[$i] = $SpliceGraphs[$i].';'.$hw_outfname;
		} else {
		    $HWInputNames[$i] = $hw_infname;
		    $SpliceGraphs[$i] = $hw_outfname;
		}
	    }
	}
    }

    # We'll save removing the nucleotide file until the end, in case we end up
    # reusing the file name for SPALN-assisted search
    
    
    # SpliceGraphs is now an array of file names (or zeros) for HitWeaver output,
    # which we can read in to check for full-protein splicings (or else to see if
    # there's at least agreement on some canonical exons...)

    
    # First off, we'll see which sequences have full-protein mappings and
    # record the mapped chromosomes, so we can determine a canonical
    # chromosome
    my @FullMaps;
    my %FullMapsByChr;
    for (my $i=0; $i<$num_seqs; $i++) {
	if ($SpliceGraphs[$i]) {

	    $FullMaps[$i] = CheckForFullMaps($SpliceGraphs[$i],length($Seqs[$i]));

	    # If there weren't any full mappings, skip to the next sequence.
	    # Otherwise, get the names of the chromosomes (from the file names)
	    # that had full mappings
	    next if (!$FullMaps[$i]);

	    foreach my $hw_fname (split(/\;/,$FullMaps[$i])) {

		$hw_fname =~ /\.([^\.]+)\.weaver\.out/;
		my $chr = $1;

		if ($FullMapsByChr{$chr}) { $FullMapsByChr{$chr}++; }
		else                      { $FullMapsByChr{$chr}=1; }

	    }
	    
	} else {
	    $FullMaps[$i] = 0;
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
	my @TopChrMaps;
	my @SeqHitStrs;
	for (my $i=0; $i<$num_seqs; $i++) {

	    # You gotta earn that hit string!
	    $SeqHitStrs[$i] = 0;

	    # Determine which camp we fall in
	    my $top_chr_fname = 0;
	    if ($SpliceGraphs[$i] =~ /([^\;]+\.$top_chr\.weaver\.out)/) {
		$top_chr_fname = $1;
		$TopChrMaps[$i] = $top_chr_fname;
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
		$SeqHitStrs[$i] =
		    ExtractCanonHWMap($TopChrMaps[$i],$seq_len,$SeqNames[$i],$top_chr);

	    } elsif ($top_chr_fname) {

		# We know this sequence is at least partially hitting to the
		# expected chromosome.  We'll see if we can use SPALN to find
		# putative exons to fill in the gaps in the 
		my ($max_exon_str,$max_amino_str,$max_hits_str) =
		    FindMaximalHWHitSet($SpliceGraphs[$i],$top_chr_fname);

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
		# SPALN putative exon ranges, but then plug them into a new HitWeaver
		# run, since we want clean splice boundary delineation between our
		# exons.
		#
		$SeqHitStrs[$i] =
		    AttemptSpalnFill($HWInputNames[$i],$top_chr_fname,$max_exon_str,
				     $max_hits_str,\@GapStarts,\@GapEnds,$Seqs[$i],
				     $SeqNames[$i],$genome,$gene_fname);

		# Alriiiight, let's see if we can make a monster!
		# We'll either get 0 if we didn't get a full mapping, or else
		# we'll get the hit string.
		$SeqHitStrs[$i] =
		    AttemptChimericHWMap($HWInputNames[$i],$SpliceGraphs[$i],
					 $top_chr_fname,$max_exon_str,\@GapStarts,
					 \@GapEnds,$Seqs[$i],$SeqNames[$i],$gene_fname);

		# At this point, we either have a full hit to multiple chromosomes
		# (TRANS SPLICING!) or we don't, as is the way of logic
		#
		# If we don't have a full mapping, we'll do some prep work for BLAT.
		# In addition to writing out the full sequence, we'll also put out
		# the parts of the sequence that correspond to the gaps left by our
		# partially mapping.
		#
		if (!$SeqHitStrs[$i]) {

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
			print $BlatFile ">$gene\|$SeqNames[$i]\n";
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
		    my $max_hits_fname = $top_chr_fname;
		    $max_hits_fname =~ s/\.weaver\.out$/\-partial\.weaver\.out/;
		    RecordMaximalHits($top_chr_fname,$max_hits_fname,$max_hits_str);
		    
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
	    
	}

	# We can at least say that we got one sequence to fully map, so let's pop
	# some champagne! Or, absent champagne, open up... [gene].quilter.out!
	my $outfname = $gene_fname;
	$outfname =~ s/\.fa$/\.quilter\.out/;
	my $OutFile = OpenOutputFile($outfname);
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($SeqHitStrs[$i]) {
		print $OutFile "$SeqHitStrs[$i]\n";
	    }
	}
	close($OutFile);

    } else {

	# Straightforward copy of our sequences into the BLAT file.
	# The only minor change is that we're going to change the names
	# to be gene|seqname
	for (my $i=0; $i<$num_seqs; $i++) {
	    print $BlatFile ">$gene\|$SeqNames[$i]\n";
	    my @Seq = split(//,$Seqs[$i]);
	    for (my $j=0; $j<scalar(@Seq); $j++) {
		print $BlatFile "$Seq[$j]";
		print $BlatFile "\n" if (($j+1) % 60 == 0);
	    }
	    print $BlatFile "\n" if (scalar(@Seq) % 60 != 0);
	    print $BlatFile "\n";
	}

    }

    # Clear out the files -- if you need to save them, return
    #return; # DEBUGGING
    for (my $i=0; $i<$num_seqs; $i++) {
	foreach my $fname (split(/\;/,$SpliceGraphs[$i])) {
	    RunSystemCommand("rm \"$fname\"") if ($fname);
	}
	foreach my $fname (split(/\;/,$HWInputNames[$i])) {
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
	    foreach my $start_amino (sort {$a <=> $b} keys %StartAminoToSeqs) {
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

    # Now we can prep a lil' file for Weaver (working title)!
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

	# Extract the nucleotide sequence for this hit
	my @HitData = split(/\|/,$hit);
	my $nucl_str = uc($HitData[7]);

	# Compute the strengths of the splice sites
	my ($three_prime_str,$five_prime_str) = GetSpliceStrengths($nucl_str);

	# Now that we have that additional whiff of data in the mix,
	# let's get writing!
	print $WeaverFile "\n";
	print $WeaverFile "Hit Score   : $HitData[6]\n";
	print $WeaverFile "Amino Range : $HitData[0]\.\.$HitData[1]\n";
	print $WeaverFile "Nucl Range  : $HitData[2]\.\.$HitData[3]\n";
	print $WeaverFile "Nucleotides : $nucl_str\n";
	print $WeaverFile "3' SS Str.  : $three_prime_str\n";
	print $WeaverFile "5' SS Str.  : $five_prime_str\n";
	
    }


    # Close up the file so we can put it to work!
    close($WeaverFile);

    # Also, in case you're like me, and you've hit this point and started wondering,
    # "Hey, why aren't we just piping FastMap2 output into HitWeaver?" the main reason
    # is that we needed to divide FastMap2 output according to the protein sequences
    # and then sort them so that we have a guarantee of ascending start aminos.

    # We'll return the name of this file
    my $weaver_out = $weaver_in;
    $weaver_out =~ s/\.in$/\.out/;

    # PUT IT TO WORK!!!
    my $weaver_cmd = $srcdir."HitWeaver --report-singles \"$weaver_in\" > \"$weaver_out\"";
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
    foreach my $hw_outfname (split(/\;/,$splicegraphs)) {

	my $hwf = OpenInputFile($hw_outfname);
	while (my $line = <$hwf>) {
	    $line =~ s/\n|\r//g;
	    if ($line =~ /Amino Acid Range : 1\.\.$seq_len/) {
		# GOT ONE!
		if ($full_maps) { $full_maps = $full_maps.';'.$hw_outfname; }
		else            { $full_maps = $hw_outfname;                }
		last;
	    }
	}
	close($hwf);
	
    }

    # Couldn't have been easier!
    return $full_maps;

}





############################################################
#
#  Function: ExtractCanonHWMap
#
#  NOTE: This function assumes that we have a full mapping to a single chromosome.
#        Any full mappings using chimeric inputs to HitWeaver will need to use
#        ExtractChimericHWMap.
#
sub ExtractCanonHWMap
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
    my $hw_inf_str        = shift;
    my $hw_outfname       = shift;
    my $max_exon_list_str = shift;
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
    $hw_outfname =~ /\/([^\/]+)$/;
    my $hitfname = $1;
    $hitfname =~ /^[^\.]+\.([^\.]+)\./;
    my $chr = $1;
    my $revcomp = 0;
    if ($chr =~ /\-revcomp$/) {
	$chr =~ s/\-revcomp$//;
	$revcomp = 1;
    }

    # We'll kick things off by opening up the HitWeaver output file and grabbing the
    # nucleotide ranges corresponding to each of our favorite hits.
    my $hwoutf = OpenInputFile($hw_outfname);
    my @MaxHits = split(/\,/,$max_hits_list_str); # These are [0..num_hits-1]
    my @MaxHitNuclRanges;
    my $hit_num = 0;
    my $exon_num = 0;
    while (my $line = <$hwoutf>) {

	if ($line =~ /\-\-\-\-\- Exons/) {
	    if ($exon_num == $MaxHits[$hit_num]) {

		# Scan to where the nucleotide range is hiding
		$line = <$hwoutf>; # blank line
		$line = <$hwoutf>; # exon id list
		$line = <$hwoutf>; # amino range
		$line = <$hwoutf>; # Nucleotide Range!
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
    close($hwoutf);

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

    my @GapExonRanges;
    for (my $i=0; $i<$num_gaps; $i++) {

	open(my $tempf,'>',$temp_fname);
	print $tempf ">$i\n";
	for (my $j=$GapStarts[$i]-1; $j<$GapEnds[$i]; $j++) {
	    print $tempf "$Seq[$j]";
	}
	print $tempf "\n\n";
	close($tempf);
	
	my $nucl_range = $NuclSearchRanges[$i];
	my $sfetch_cmd = $sfetch." -c $nucl_range \"$genome\" \"$chr\" > $nucl_fname";
	RunSystemCommand($sfetch_cmd);

	my $spaln_cmd = $spaln." -Q3 -O1 -S3 -ya3 \"$nucl_fname\" \"$temp_fname\"";
	
    }
    
    
    
}






############################################################
#
#  Function: AttemptChimericHWMap
#
sub AttemptChimericHWMap
{
    my $orig_hw_inf_str   = shift;
    my $orig_hw_outf_str  = shift;
    my $canon_chr_fname   = shift;
    my $max_exon_list_str = shift;
    my $amino_starts_ref  = shift;
    my $amino_ends_ref    = shift;
    my $seq_str           = shift;
    my $seq_id            = shift;
    my $gene_fname        = shift;

    my @HWOutNames     = split(/\;/,$orig_hw_outf_str);
    my @HWInNames      = split(/\;/,$orig_hw_inf_str);
    my @AminoGapStarts = split(/\,/,$amino_starts_ref);
    my @AminoGapEnds   = split(/\,/,$amino_ends_ref);

    # I'm going to build a list of all the exons we want to use as inputs to
    # HitWeaver.  We'll need a method for sorting these by their starting aminos, too.
    my @FullExonList;
    my %ExonInfoToListIndex;
    my @ChrsByExon;

    # Next up, we'll run through each of the non-canonical chromosome's input
    # files and find all of the putative exons that fall into the gaps between
    # canonical hits
    my $num_exons = 0;
    for (my $i=0; $i<scalar(@HWInNames); $i++) {

	# Which chromosome is this?
	$HWInNames[$i] = /\.([^\.]+)\.weaver\.in$/;
	my $chr = $1;

	# We can treat the canonical chromosome somewhat specially, since it has
	# a list pre-generated
	my @ExonList;
	if ($HWOutNames[$i] eq $canon_chr_fname) {

	    my $exon_list_ref
		= ExtractHWInputExons($canon_chr_fname,$max_exon_list_str);

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
	    my $hwinf = OpenInputFile($HWInNames[$i]);
	    my $gap_num = 0;
	    while (my $line = <$hwinf>) {

		# Exon is kicking off!
		if ($line =~ /^Hit Score/) {

		    # Begin the exon tracking
		    my $exon = $line;

		    $line = <$hwinf>;
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
			$line = <$hwinf>; # Nucl. Range
			$exon = $exon.$line;
			$line = <$hwinf>; # Nucl. String
			$exon = $exon.$line;
			$line = <$hwinf>; # 3' SS Strengths
			$exon = $exon.$line;
			$line = <$hwinf>; # 5' SS Strengths
			$exon = $exon.$line."\n";

			push(@FullExonList,$exon);
			push(@ChrsByExon,$chr);

			my $key = $start_amino.'|'.$end_amino.'|'.$i;

			$num_exons++;
			$ExonInfoToListIndex{$key} = $num_exons;
			
		    }

		}

	    }

	    close($hwinf);

	}
	
    }

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

    # Close up the file and give it to HitWeaver!  Be sure to let them know
    # that you don't want to worry about chromosome consistency...
    #
    # NOTE: We aren't interested in single-exon hits this time -- if we don't
    #       get any decent splicing on the canonical exon, we'll pass this sequence
    #       entirely over to BLAT+SPALN
    #
    close($ChimeraFile);
    my $weaver_cmd = $srcdir."HitWeaver --allow-inconsistency \"$chimera_in\" > \"$chimera_out\"";
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
	$full_map = ExtractChimericHWMap($chimera_out,\@ChrsByExon,$seq_len,$seq_id);
    }

    # No more need for this file -- burn it and punt the full map (if we have one)
    RunSystemCommand("rm \"$chimera_out\"") if (-e $chimera_out);
    return $full_map;

}




############################################################
#
#  Function: ExtractHWInputExons
#
sub ExtractHWInputExons
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
#  Function: ExtractChimericHWMap
#
#  NOTE: To give the opposite note from ExtractCanonHWMap, this function doesn't
#        assume that all hits are to the same strand / chromosome.
#
sub ExtractChimericHWMap
{
    my $hw_outfname = shift;
    my $chr_list_ref = shift;
    my $seqlen = shift;
    my $seq_id = shift;

    my @ChrsByExon = @{$chr_list_ref};
    
    # Open the file and scan until we we're in the full mapping zone (FMZ)
    my $inf = OpenInputFile($hw_outfname);
    my $line = <$inf>;
    while ($line = <$inf>) {
	$line =~ s/\n|\r//g;
	last if ($line =~ /Amino Acid Range \: 1\.\.$seqlen$/);
    }

    # I'm not going to tolerate abuse of this function!
    if (eof($inf)) { close($inf); die "\n  ERROR:  '$hw_outfname' didn't fully map!\n\n"; }
    
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
    $hitstr = $hitstr."Chromosome : CHIMERIC!\n";
    $hitstr = $hitstr."Num Exons  : $num_exons\n";

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

    print "$hitstr\n";
    
}






############################################################
#
#  Function: FindMaximalHWHitSet
#
sub FindMaximalHWHitSet
{
    my $outf_list_str    = shift;
    my $target_chr_fname = shift;

    my @OutFileNames = split(/\;/,$outf_list_str);

    # First, find which files correspond to our target chromosome
    my $target_index = 0;
    while ($OutFileNames[$target_index] ne $target_chr_fname) {
	$target_index++;
    }

    # Open up the target chromosome's HW output file and pull in
    # a list of the hits our sequence had to it, with exon indices.
    # Keep in mind that these exon indices are [1..num_exons],
    # with respect to the order of exons in the input file.
    my @TargetHitStarts;
    my @TargetHitEnds;
    my @TargetHitExons;
    my $hwoutf   = OpenInputFile($OutFileNames[$target_index]);
    my $num_hits = 0;
    while (my $line = <$hwoutf>) {

	next if ($line !~ /\-\-\-\-\- Exons/);

	$line = <$hwoutf>; # blank line
	$line = <$hwoutf>; # Exon list!
	$line =~ /Spliced Exon IDs \: (\S+)/;
	$TargetHitExons[$num_hits] = $1;

	$line = <$hwoutf>; # score line
	$line = <$hwoutf>; # Amino range!
	$line =~ /Amino Acid Range \: (\d+)\.\.(\d+)/;
	$TargetHitStarts[$num_hits] = $1;
	$TargetHitEnds[$num_hits]   = $2;

	$num_hits++;
	
    }
    close($hwoutf);

    # Where are there holes in our mapping?
    # NOTE: Because the output method in HW runs linearly through the list of exons
    #       to see which start hits, they should be pre-sorted according to amino
    #       start position, so we can take advantage of that when determining the
    #       set of non-overlapping hits that maximally cover the protein
    my @MaximalHitSet;
    push(@MaximalHitSet,0); # Prime it!
    my $num_top_hits = 1;
    my $top_recent_hit = 0;
    for (my $i=0; $i<$num_hits; $i++) {

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
    for (my $i=0; $i<$num_hits; $i++) {
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
#  Function: RecordMaximalHits
#
sub RecordMaximalHits
{
    my $infname = shift;
    my $outfname = shift;
    my $hit_index_list = shift;

    # Open up the original HitWeaver output file and the file we'll be transferring
    # a subset of those hits into
    my $inf  = OpenInputFile($infname);
    my $outf = OpenOutputFile($outfname);

    # HitWeaver output files always start with a blank line, and who am I to
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

    my $gene_name_loc = shift; # 'seqname' or 'filename'
    my $seqfiles_ref  = shift;
    my $genome        = shift;

    my @SeqFileNames = @{$seqfiles_ref};
    
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
	my $BlatFile = OpenOutputFile($fname);

	# As we parse this file, we'll want to track what gene family each sequence
	# belongs to, along with its name
	my ($gene,$seqname);
	
	# If we're getting the gene name from the filename, it's regex time!
	if ($gene_name_loc eq 'filename') {
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
		if ($gene_name_loc eq 'seqname') {
		    $seqname =~ /^([^\|]+)\|([^\|]+)$/;
		    $gene    = $1;
		    $seqname = $2;
		}
		$num_blat_seqs++;
		print $CumBlatFile ">$num_blat_seqs\n";

		# Who are you, again?
		$BlatNameGuide{$num_blat_seqs} = $seqname;
		
		# If this is the first time we've seen a member of this gene family,
		# take note
		if (!$BlatGenes{$gene}) {
		    $BlatGenes{$gene} = 1;
		    $num_blat_genes++;
		}

	    } else {

		# Just a regular ol' sequence line
		print $CumBlatFile "$line\n";
		
	    }
	}

	close($BlatFile);

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
    return($blat_outfname,\%BlatNameGuide,$num_blat_genes);

}






# EOF











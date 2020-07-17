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
sub PrintHitToWeaverInf;
sub GetSpliceStrengths;
sub CheckForFullMaps;
sub ExtractCanonHWMap;
sub AttemptSpalnFill;
sub AttemptChimericHWMap;
sub ExtractHWInputExons;
sub ExtractChimericHWMap;
sub FindMaximalHWHitSet;
sub SaveTopHWInputs;
sub GetAminoHitGaps;
sub RecordMaximalHits;
sub RunBlatOnFileSet;
sub GenBlatMaps;
sub AttemptBlatFill;
sub BlatToSpalnSearch;
sub SpalnSearchChr;
sub GroupRanges;
sub AdjustNuclCoord;
sub ParseBlatLine;
sub ParseSpalnOutput;
sub FinalFileCheck;



# Help?
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
my $spaln  = $srcdir.'../inc/spaln2.3.3/src/spaln';

# Spaln requires certain environment variables to be set, so why don'tcha set 'em?
$spaln =~ /^(.*)spaln$/;
$ENV{'ALN_TAB'} = $1.'../table';
$ENV{'ALN_DBS'} = $1.'../seqdb';

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


# If we're using a GTF for our first round of searching, hop to it!
# Otherwise, jump straight to building up one big file for BLAT search.
my $blat_naming;
my @BlatFileNames;
if ($gtfname ne '-') {

    UseGTF($seq_dirname,$ali_dirname,$genome,$gtfname,$num_cpus,\%Opts);

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

    # Cleanup on line 135! << CRITICAL: KEEP THIS LINE NUMBER CORRECT, FOR JOKE
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

    # GTF used!
    
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

	    # Buff it up! Until you can feel it!
	    $search_start += 50;
	    $search_end   -= 50;
	    
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

	    # Buff it up! And you don't even need it!
	    $search_start -= 50;
	    $search_end   += 50;
	    
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
		    AttemptSpalnFill($top_chr_fname,$max_hits_str,\@GapStarts,\@GapEnds,
				     $Seqs[$i],$SeqNames[$i],$genome,$gene_fname);

		# At this point, if we still don't have a full mapping, it seems like
		# something nonstandard could be going on, so we'll quickly check if
		# the full set of hits (across chromosomes) could be useful...
		if (!$SeqHitStrs[$i]) {
		    $SeqHitStrs[$i] =
			AttemptChimericHWMap($HWInputNames[$i],$SpliceGraphs[$i],
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
		    # HitWeaver -- this could be optimized, but I'm pretty sure we'd
		    # see diminishing returns...
		    my $top_hwinfname = $top_chr_fname;
		    $top_hwinfname =~ s/\.out$/\.in/;
		    SaveTopHWInputs($top_hwinfname);
		    
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
	print $WeaverFile "\n";
	PrintHitToWeaverInf($hit,$WeaverFile);
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
	    die "\n  Scan error: $fname\n  ($scan in '@Nucls')\n\n" if ($scan >= scalar(@Nucls)); # DEBUGGING
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
    my $hw_outfname       = shift;
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
		$line = <$hwoutf>; # mapping score
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
	my $sfetch_cmd = $sfetch." -c $nucl_range \"$genome\" \"$chr\" > $nucl_fname";
	RunSystemCommand($sfetch_cmd);

	my $spaln_cmd = $spaln." -Q3 -O1 -S3 -ya3 \"$nucl_fname\" \"$temp_fname\" 2>/dev/null";
	my $SpalnOut  = OpenSystemCommand($spaln_cmd);

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
    RunSystemCommand("rm \"$temp_fname\"") if (-e $temp_fname);
    RunSystemCommand("rm \"$nucl_fname\"") if (-e $nucl_fname);
    return 0 if (!scalar(@GapHits));
    

    # Uber chill, bro!  Looks like we've got (hopefully) a list of
    # nucleotide coordinates corresponding to the presumptive exons
    # that we were able to get out of running SPALN.
    #
    # Next up, we'll pull in the original FastMap2 hits and add them
    # to GapHits.  This will give us a full input to push off to
    # HitWeaver, which will hopefully give us a full-sequence mapping...


    # Open up the original HitWeaver input file
    my $hw_infname = $hw_outfname;
    $hw_infname =~ s/out$/in/;
    my $hwinf = OpenInputFile($hw_infname);

    while (my $line = <$hwinf>) {

	# Is this the start of a new exon?
	if ($line =~ /Hit Score   \: (\S+)/) {

	    my $hit_score = $1 + 0.0;

	    $line = <$hwinf>; # Amino range
	    $line =~ /Amino Range \: (\d+)\.\.(\d+)/;
	    my $hit_amino_start = $1;
	    my $hit_amino_end   = $2;
	    
	    $line = <$hwinf>; # Nucl  range
	    $line =~ /Nucl Range  \: (\d+)\.\.(\d+)/;
	    my $hit_nucl_start = $1;
	    my $hit_nucl_end   = $2;
	    
	    $line = <$hwinf>; # Nucleotide Seq.
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

    close($hwinf);

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
    RunSystemCommand("rm \"$hw_infname\"");
    RunSystemCommand("rm \"$hw_outfname\"");

    # Push off to HitWeaver!
    ($hw_infname,$hw_outfname) =
	GenSpliceGraph($final_exon_str,$seq_str,$seq_id,$chr,$revcomp,$gene_fname);

    # Did we get a full mapping?
    my $hit_str = CheckForFullMaps($hw_outfname,$seq_len);

    # If we did get a full mapping, give credit to SPALN for filling in the gaps
    if ($hit_str) {
	$hit_str =~ s/Map Method \: FastMap2/Map Method \: FastMap2\+Spaln/;
    }

    # DONE
    return $hit_str;

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
	$HWInNames[$i] =~ /\.([^\.]+)\.weaver\.in$/;
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
#  Function: SaveTopHWInputs
#
sub SaveTopHWInputs
{
    my $hw_infname = shift;

    # Get the name without any of the directory path stuff
    $hw_infname =~ /\/([^\/]+)$/;
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
    my $inf = OpenInputFile($hw_infname);
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

    my $seqfiles_ref  = shift;
    my $genome        = shift;
    my $gene_name_fmt = shift; # 'seqname' or 'filename'

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
    my $gene_index_start =  $threadID    * ($num_blat_genes / $num_cpus);
    my $gene_index_end   = ($threadID+1) * ($num_blat_genes / $num_cpus);
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
	    
	    # If we have partial information from HitWeaver (derived from GTF)
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

    # Find the file with the stored HitWeaver inputs to this gene's canon chromosome
    my $hwinfname = $seq_dirname.$seqname.'.partial.tmp';

    # Open it up and pull in the original hits, along with the chromosome info
    my $OrigHWInf = OpenInputFile($hwinfname);

    my $orig_chr = <$OrigHWInf>;
    $orig_chr =~ s/\n|\r//g;
    $orig_chr =~ s/^Chromosome\: //;
    $orig_chr =~ s/\[revcomp\]/\-revcomp/;

    my @OrigHits;
    while (my $line = <$OrigHWInf>) {

	if ($line =~ /Hit Score/) {

	    $line =~ /Hit Score   \: (\S+)/;
	    my $hit_score = $1;

	    $line = <$OrigHWInf>;
	    $line =~ /Amino Range \: (\d+)\.\.(\d+)/;
	    my $amino_start = $1;
	    my $amino_end   = $2;
	    
	    $line = <$OrigHWInf>;
	    $line =~ /Nucl Range  \: (\d+)\.\.(\d+)/;
	    my $nucl_start = $1;
	    my $nucl_end   = $2;
	    
	    $line = <$OrigHWInf>;
	    $line =~ /Nucleotides \: (\S+)/;
	    my $nucl_str = $1;

	    my $hitstr = $amino_start.'|'.$amino_end.'|'.$nucl_start.'|'.$nucl_end;
	    $hitstr    = $hitstr.'|-|-|'.$hit_score.'|'.$nucl_str;
	    push(@OrigHits,$hitstr);
	    
	}
	
    }
    close($OrigHWInf);

    # We'll clear this file out, since it's not going to be used for anything else
    RunSystemCommand("rm \"$hwinfname\"");

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

	my $sfetchcmd = $sfetch." -c $sfetch_range \"$genome\" \"$blat_chr\"";
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
	# NOTE that we'll be passing this off to ExtractChimericHWMap, which
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
    my @ChrsByExon; # Keeping with the nomenclature of 'ExtractChimericHWMap'
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

    # Now we can convert our hits into HitWeaver input! We'll reuse the filename
    # that stored the original hits.
    my $num_exons  = scalar(@FullHitList);
    my $seq_len    = length($seq_str);
    my $WeaverFile = OpenOutputFile($hwinfname);
    print $WeaverFile "Num Hits : $num_exons\n";
    print $WeaverFile "Seq Len  : $seq_len\n";
    print $WeaverFile "Sequence : $seq_str\n";
    foreach my $hit (@FullHitList) {
	print $WeaverFile "\n";
	PrintHitToWeaverInf($hit,$WeaverFile);
    }
    close($WeaverFile);

    # Come up with a hip name for the output file
    my $hwoutfname = $hwinfname;
    $hwoutfname =~ s/\.tmp$/\.out/;

    # NOTE: Even though it's possible that our BLAT results will only have found
    #       hits to a single chromosome, we'll use the chimeric search to cover all
    #       possible bizzare splicing patterns.
    my $weaver_cmd = $srcdir."HitWeaver --allow-inconsistency \"$hwinfname\" > \"$hwoutfname\"";
    RunSystemCommand($weaver_cmd);

    # Begone, input file!
    RunSystemCommand("rm \"$hwinfname\"");

    # Did we get anything at all?
    if (!(-s $hwoutfname)) {
	RunSystemCommand("rm \"$hwoutfname\"") if (-e $hwoutfname);
	return 0;
    }

    # You know the way this goes
    my $full_map = CheckForFullMaps($hwoutfname,$seq_len);
    if ($full_map) {
	$full_map = ExtractChimericHWMap($hwoutfname,\@ChrsByExon,$seq_len,$seqname);
	# Credit to BLAT!
	$full_map =~ s/FastMap2/BLAT\+FastMap2/;
    }
    
    # That's all!
    RunSystemCommand("rm \"$hwoutfname\"");
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

    # We'll start off by building up a hash of hits, broken apart by chromosome
    my @BlatHits = @{$blathits_ref};
    my %BlatHitsByChr;
    for (my $i=0; $i<scalar(@BlatHits); $i++) {

	# We'll end up clearing out the first pieces of data in the blat hits,
	# since we know there aren't partial hits in this list.
	$BlatHits[$i] =~ /^\-\s+(.*)$/;
	my $blat_outline = $1;
	my ($blat_chr,$blat_amino_start,$blat_amino_end,$blat_nucl_start,$blat_nucl_end,$blat_score)
	    = ParseBlatLine($blat_outline);

	# We'll want to know if this is the reverse complement
	my $revcomp = 0;
	$revcomp = 1 if ($blat_nucl_start > $blat_nucl_end);

	# Hash it!
	my $chr = $blat_chr.'-'.$revcomp;
	my $hit = $blat_amino_start.'|'.$blat_amino_end.'|'.$blat_nucl_start;
	$hit    = $hit.'|'.$blat_nucl_end.'|'.$blat_score;
	if ($BlatHitsByChr{$chr}) {
	    $BlatHitsByChr{$chr} = $BlatHitsByChr{$chr}.'&'.$hit;
	} else {
	    $BlatHitsByChr{$chr} = $hit;
	}
	
    }

    # We'll go ahead and make a file with the protein sequence
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

    # Now we can go through each of the chromosomes and see how high-quality of a
    # mapping we can get out of Spaln
    my @HitsByChr;
    my @HitPctsID;
    my @ChrList;
    my @FullHitIndices;
    my $tot_exons = 0;
    my $num_hits = 0;
    foreach my $blat_chr (keys %BlatHitsByChr) {

	$blat_chr =~ /^(\S+)\-(\d)$/;
	my $chr = $1;
	my $revcomp = $2;

	my @HitList = split(/\&/,$BlatHitsByChr{$blat_chr});
	my ($hitstr,$pct_id)
	    = SpalnSearchChr($seqname,\@Seq,$chr,$revcomp,\@HitList,$prot_fname);

	# Did we get a little ol' snack?
	next if (!$hitstr);

	# Note that these hits are fully formatted (so they all start with the whole
	# "sequence id," "map method," "chromosome," and "num exons" dump).
	$num_hits++;
	push(@HitsByChr,$hitstr);
	push(@HitPctsID,$pct_id);
	push(@ChrList,$blat_chr);

	# Quick question for ya: Does this hit give us full coverage?

	# Pull out the data and gather all of the amino ranges
	my @HitLines = split(/\n/,$hitstr);
	$HitLines[3] =~ /Num Exons  \: (\d+)/;
	my $num_exons = $1;
	my @ExonAminos;
	for (my $i=4; $i<4+($num_exons*2); $i+=2) {
	    $HitLines[$i] =~ /Aminos (\d+\.\.\d+)\,/;
	    push(@ExonAminos,$1);
	}

	# Increment the total number of exons across all hits
	$tot_exons += $num_exons;

	# We can initially check if the first exon starts the protein and the final
	# exon ends the protein
	my $seq_len = scalar(@Seq);
	if ($ExonAminos[0] =~ /^1\.\./ && $ExonAminos[$num_exons-1] =~ /\.\.$seq_len$/) {
	    # Run through and see if we have continuity...
	    my $continuity = 1;
	    for (my $i=1; $i<$num_exons; $i++) {
		$ExonAminos[$i-1] =~ /\.\.(\d+)$/;
		my $last_end = $1;
		$ExonAminos[$i] =~ /^(\d+)\.\./;
		my $this_start = $1;
		if ($last_end != $this_start-1) {
		    $continuity = 0;
		    last;
		}
	    }

	    # If we have continuity, this is a full hit!
	    push(@FullHitIndices,$num_hits) if ($continuity);
	    
	}
	

    }

    # No more searching for you, mister/miss!
    RunSystemCommand("rm \"$prot_fname\"");

    # If we didn't get anything worth reporting, we'll just have to cry all night
    return 0 if ($num_hits == 0);

    
    # TOTALLY TUBULAR!
    #
    # If we have one full-length mapping, we can just pass it up the chain.
    #
    # If we have multiple full-length mappings, pass up the one with the
    # highest percent-identity.
    #
    # If we haven't identified a solid full-sequence mapping, we'll mark the 
    # chromosome as 'Incomplete' and dump all of our exons as a giant mess.
    #
    if (scalar(@FullHitIndices)) {

	# Find the highest percent-identity hit
	my $top_hit_id = $FullHitIndices[0];
	my $top_pct_id = $HitPctsID[$top_hit_id];
	for (my $i=1; $i<scalar(@FullHitIndices); $i++) {
	    my $hit_id = $FullHitIndices[$i];
	    my $pct_id = $HitPctsID[$i];
	    if ($pct_id > $top_pct_id) {
		$top_hit_id = $hit_id;
		$top_pct_id = $pct_id;
	    }
	}

	# WINNER!
	return $HitsByChr[$top_hit_id];
	
    } else {

	# FRANKENSTEIN!
	my $final_hit = "Sequence ID: $seqname\n";
	$final_hit = $final_hit."Map Method : BLAT\+Spaln\n";
	$final_hit = $final_hit."Chromosome : Incomplete\n";
	$final_hit = $final_hit."Num Exons  : $tot_exons\n";
	for (my $i=0; $i<$num_hits; $i++) {
	    my @HitLines = split(/\n/,$HitsByChr[$i]);
	    for (my $j=4; $j<scalar(@HitLines); $j++) {
		$final_hit = $final_hit.$HitLines[$j]."\n";
	    }
	}
	return $final_hit;

    }
    
}








############################################################
#
#  Function: SpalnSearchChr
#
sub SpalnSearchChr
{
    
    my $seqname = shift;
    my $seq_ref = shift;
    my $chr     = shift;
    my $revcomp = shift;
    my $blathits_ref = shift;
    my $prot_fname = shift;

    my @Seq = @{$seq_ref};
    my @BlatHits = @{$blathits_ref};

    # NOTE: This is somewhat overcomplicated-looking, but I think it's the right
    #       approach to pulling in nucleotides for our Spaln search.  Basically,
    #       we try to get 

    # Start off by computing the per-position BLAT coverage
    my $seq_len = scalar(@Seq);
    my @BlatCoverage;
    my @CoverageID; # Who's claiming this sequence segment?
    for (my $i=0; $i<$seq_len; $i++) {
	$BlatCoverage[$i] = 0;
	$CoverageID[$i]   = -1;
    }

    for (my $i=0; $i<scalar(@BlatHits); $i++) {
	my @HitData = split(/\|/,$BlatHits[$i]);
	for (my $j=$HitData[0]-1; $j<$HitData[1]; $j++) {
	    $BlatCoverage[$j]++;
	    $CoverageID[$j]=$i;
	}
    }

    # We'll mark gap-adjacent positions and overlap positions in a way that provides
    # a clean ontology.
    my @GapAdjacency;
    my $num_gaps = 0;
    for (my $i=0; $i<$seq_len; $i++) {

	# Are we preceding a gap?
	if ($BlatCoverage[$i] == 0) {

	    # We're seeing a new gap!
	    $num_gaps++;

	    # Hey! You were adjacent the whole time!
	    $GapAdjacency[$i-1] = $num_gaps if ($i); 

	    # Who else is adjacent to this gap run? (or a member)
	    $GapAdjacency[$i++] = $num_gaps;
	    while ($i<$seq_len && $BlatCoverage[$i]==0) {
		$GapAdjacency[$i++] = $num_gaps;
	    }
	    
	    # I see you being adjacent over there!
	    $GapAdjacency[$i] = $num_gaps if ($i<$seq_len);
	    
	} else {
	    $GapAdjacency[$i] = 0;
	}
	
    }

    my @Overlaps;
    my $num_overlaps = 0;
    for (my $i=0; $i<$seq_len; $i++) {

	# Are we in an overlapping region?
	if ($BlatCoverage[$i] > 1) {
	    # New overlap!
	    $num_overlaps++;
	    while ($i<$seq_len && $BlatCoverage[$i] > 1) {
		$Overlaps[$i++] = $num_overlaps;
	    }
	}

	# Note that we'll have overstepped by one if we were just in an overlap
	# region, so we don't use an 'else'
	$Overlaps[$i] = 0;
	    
    }

    # Because we could have gaps and overlaps right next to one another,
    # we'll do a bit of final ontologizing by simply making groups where
    # there are runs of unclear coverage.
    my @HitGrouping;
    my @HitGroupStarts; # inclusive
    my @HitGroupEnds;   # exclusive
    my $num_hit_groups = 0;
    for (my $i=0; $i<$seq_len; $i++) {
	if ($GapAdjacency[$i] || $Overlaps[$i]) {
	    $num_hit_groups++;
	    push(@HitGroupStarts,$i);
	    while ($i<$seq_len && ($GapAdjacency[$i] || $Overlaps[$i])) {
		$HitGrouping[$i++] = $num_hit_groups;
	    }
	    push(@HitGroupEnds,$i);
	}
	$HitGrouping[$i] = 0;
    }

    # Next up, we'll figure out the full span of chromosomal sequence that
    # is implicated in groups of hits.  If the range is large (>5Mb) then
    # do some work to split them into more manageable sequence chunks.
    my @GroupNuclRanges;
    for (my $i=0; $i<$num_hit_groups; $i++) {

	my $group_start = $HitGroupStarts[$i];
	my $group_end   = $HitGroupEnds[$i];

	# Which blat hits are implicated in this group?
	# NOTE that even if there's a gap starting or ending the sequence, we're
	# guaranteed to touch at least one hit, so we don't need to worry about
	# getting too funky.
	my @GroupedHits;
	my @GroupNuclStarts;
	my @GroupNuclEnds;
	for (my $j=0; $j<scalar(@BlatHits); $j++) {
	    my @HitData = split(/\|/,$BlatHits[$j]);
	    # Check the first and final positions to see if we're in this group
	    if ($HitGrouping[$HitData[0]-1]==$i || $HitGrouping[$HitData[1]-1]==$i) {
		push(@GroupedHits,$j);
		push(@GroupNuclStarts,$HitData[2]);
		push(@GroupNuclEnds,$HitData[3]);
	    }
	}
	
	# What we're going to end up doing here is a series of clusterings until
	# the span of each set of grouped hits is less than a preset distance
	my $max_span = 5000000; # 5 Mb	
	my $group_ranges_str
	    = GroupRanges(\@GroupNuclStarts,\@GroupNuclEnds,$max_span,$revcomp);
	push(@GroupNuclRanges,$group_ranges_str); # These come out pre-sorted!

    }


    # WOOF, this is some madness, eh?
    #
    # Next up, we'll run through the whole sequence and build a string representing
    # the collection of ranges that provide coverage suggested-ish by BLAT
    #
    # NOTE that these are still going to be tight-ish ranges -- we'll worry about
    # extending out a bit in the next step...
    
    my $range_set_str = '';
    for (my $i=0; $i<$seq_len; $i++) { # This might be better as a 'while' loop... 

	# Are we pulling from a group?
	if ($HitGrouping[$i]) {

	    # Since we already did all that tough work, this is pretty light lifting!
	    $range_set_str = $range_set_str.','.$GroupNuclRanges[$HitGrouping[$i]-1];
	    $i = $HitGroupEnds[$HitGrouping[$i]-1]-1;

	} else {

	    # Find the hit who covers this position.  It should be an amino start,
	    # because of logic.
	    #
	    # Sadly, these aren't sorted, so for now I'm just going to risk a full
	    # scan (it's SCANdalous).
	    #
	    # (SCANdalous!)
	    #
	    my $j=0;
	    while ($j<scalar(@BlatHits)) {
		my @HitData = split(/\|/,$BlatHits[$j]);
		if ($HitData[0] == $i) {
		    $range_set_str = $range_set_str.','.$HitData[2].'..'.$HitData[3];
		    $i = $HitData[1]-1;
		    last;
		}
		$j++;
	    }
	    
	}
	
    }

    # Don't need a leading comma
    $range_set_str =~ s/^\,//;

    # Next up, we'll convert our set of tight-ish ranges into slightly wider ranges,
    # especially on the ends if they're gapped...
    my @RangeStarts;
    my @RangeEnds;
    foreach my $range (split(/\,/,$range_set_str)) {
	$range =~ /^(\d+)\.\.(\d+)$/;
	push(@RangeStarts,$1);
	push(@RangeEnds,$2);
    }
    my $num_ranges = scalar(@RangeStarts);

    # First off, if either end is gapped, pull in 1Mb
    # We'll worry about exceeding the chromosome bounds LATER
    if ($BlatCoverage[0] == 0) {
	if ($revcomp) { $RangeStarts[0] += 1000000; }
	else          { $RangeStarts[0] -= 1000000; }
    }
    if ($BlatCoverage[$seq_len-1] == 0) {
	if ($revcomp) { $RangeEnds[$num_ranges-1] -= 1000000; }
	else          { $RangeEnds[$num_ranges-1] += 1000000; }
    }

    # Nice!  Now we'll just give each bound an extra 10Kb
    if ($revcomp) {
	for (my $i=0; $i<$num_ranges; $i++) {
	    $RangeStarts[$i] += 10000;
	    $RangeEnds[$i]   -= 10000;
	}
    } else {
	for (my $i=0; $i<$num_ranges; $i++) {
	    $RangeStarts[$i] -= 10000;
	    $RangeEnds[$i]   += 10000;
	}
    }

    # Our last little bit of prep work is to merge any ranges that overlap
    # (and then make sure we're in the safety zone)
    if ($revcomp) {

	my @FinalRangeStarts;
	my @FinalRangeEnds;

	push(@FinalRangeStarts,$RangeStarts[0]);
	push(@FinalRangeEnds,$RangeEnds[0]);
	my $final_num_ranges = 1;
	for (my $i=1; $i<$num_ranges; $i++) {
	    if ($RangeStarts[$i] > $FinalRangeEnds[$final_num_ranges-1]) {
		$FinalRangeEnds[$final_num_ranges-1] = $RangeEnds[$i];
	    } else {
		push(@FinalRangeStarts,$RangeStarts[$i]);
		push(@FinalRangeEnds,$RangeEnds[$i]);
		$final_num_ranges++;
	    }
	}

	for (my $i=0; $i<$final_num_ranges; $i++) {
	    $RangeStarts[$i] = $FinalRangeStarts[$i];
	    $RangeEnds[$i]   = $FinalRangeEnds[$i];
	}
	$num_ranges = $final_num_ranges;

	$RangeStarts[0] = Min($ChrSizes{$chr},$RangeStarts[0]);
	$RangeEnds[1]   = Max(1,$RangeEnds[0]);
	
    } else {

	my @FinalRangeStarts;
	my @FinalRangeEnds;

	push(@FinalRangeStarts,$RangeStarts[0]);
	push(@FinalRangeEnds,$RangeEnds[0]);

	my $final_num_ranges = 1;
	for (my $i=1; $i<$num_ranges; $i++) {
	    if ($RangeStarts[$i] < $FinalRangeEnds[$final_num_ranges-1]) {
		$FinalRangeEnds[$final_num_ranges-1] = $RangeEnds[$i];
	    } else {
		push(@FinalRangeStarts,$RangeStarts[$i]);
		push(@FinalRangeEnds,$RangeEnds[$i]);
		$final_num_ranges++;
	    }
	}

	for (my $i=0; $i<$final_num_ranges; $i++) {
	    $RangeStarts[$i] = $FinalRangeStarts[$i];
	    $RangeEnds[$i]   = $FinalRangeEnds[$i];
	}
	$num_ranges = $final_num_ranges;

	$RangeStarts[0] = Max(1,$RangeStarts[0]);
	$RangeEnds[$num_ranges] = Min($ChrSizes{$chr},$RangeEnds[$num_ranges]);
	
    }

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
	my $jump_info = $reported_pos.':'.$RangeStarts[$i];
	push(@JumpList,$jump_info);
	$reported_pos += abs($RangeEnds[$i]-$RangeStarts[$i])+1;
    }

    

    # Wowee!  Time to pull in a bunch of nucleotides, concatenate them into one
    # great big sequence, and do a Spaln search!
    my $nucl_fname = $prot_fname;
    $nucl_fname =~ s/\.prot\.in/\.nucl\.in/;
    my $NuclFile = OpenOutputFile($nucl_fname);

    # Just for fun, let's know how wide an area we're (patchily) covering
    my $span = $RangeStarts[0].'..'.$RangeEnds[$num_ranges-1];
    print $NuclFile ">$seqname\:$chr\:$span\n";

    for (my $i=0; $i<$num_ranges; $i++) {

	my $range = $RangeStarts[$i].'..'.$RangeEnds[$i];
	my $sfetchcmd = $sfetch." -c $range \"$genome\" \"$chr\"";
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

    # Assemble the Spaln command!
    my $spaln_fname = $prot_fname;
    $spaln_fname =~ s/\.prot\.in/\.spaln\.out/;
    my $spaln_cmd = $spaln." -Q3 -O1 -S3 -ya3 \"$nucl_fname\" \"$prot_fname\" 1>\"$spaln_fname\" 2>/dev/null";
    RunSystemCommand($spaln_cmd);

    # Get that nucleotide file OUTTA HERE!
    RunSystemCommand("rm \"$nucl_fname\"");

    # Don't tell me we did all that for NOTHING?!
    return (0,0) if (!(-e $spaln_fname));
    
    # UGH, I wish I could quit you, Spaln parsing
    # Note that we'll open the file here, just so it's easier to return
    # if we run into parsing errors.
    my $SpalnFile = OpenInputFile($spaln_fname);
    my ($spaln_nucl_ranges_ref,$spaln_amino_ranges_ref,$spaln_centers_ref,$spaln_pct_id)
	= ParseSpalnOutput($SpalnFile,\@Seq);
    close($SpalnFile);

    # Get that spaln output file THE HECK OUTTA HERE!
    RunSystemCommand("rm \"$spaln_fname\"") if (-e $spaln_fname);

    # If we didn't have a good return on investment, we're done-zo!
    return (0,0) if ($spaln_pct_id < 90.0);

    # All of our nucleotide coordinates need to be adjusted...
    my @UnadjNuclRanges = @{$spaln_nucl_ranges_ref};
    my @UnadjCenters = @{$spaln_centers_ref};

    my @NuclRanges;
    my @AminoRanges = @{$spaln_amino_ranges_ref};
    my @CodonCenters;
    my $num_exons = 0;
    while ($num_exons < scalar(@UnadjNuclRanges)) {

	$UnadjNuclRanges[$num_exons] =~ /^(\d+)\.\.(\d+)$/;
	my $range_start = $1;
	my $range_end   = $2;

	$range_start = AdjustNuclCoord($range_start,\@JumpList,$revcomp);
	$range_end   = AdjustNuclCoord($range_end,\@JumpList,$revcomp);
	$NuclRanges[$num_exons] = $range_start.'..'.$range_end;

	my $codon_center_str = '';
	foreach my $coord (split(/\,/,$UnadjCenters[$num_exons])) {
	    my $codon_center = AdjustNuclCoord($coord,\@JumpList,$revcomp);
	    $codon_center_str = $codon_center_str.','.$codon_center;
	}
	$codon_center_str =~ s/^\,//;
	$CodonCenters[$num_exons] = $codon_center_str;

	$num_exons++;
    }

    # Finally, let your hair down!
    $chr = $chr.'[revcomp]' if ($revcomp);

    # Build up the hit string for this chromosome
    my $hitstr = "Sequence ID: $seqname\n";
    $hitstr    = $hitstr."Map Method : BLAT+Spaln\n";
    $hitstr    = $hitstr."Chromosome : $chr\n";
    $hitstr    = $hitstr."Num Exons  : $num_exons\n";
    for (my $i=0; $i<$num_exons; $i++) {
	$hitstr = $hitstr."* Aminos $AminoRanges[$i], $chr:$NuclRanges[$i]\n";
	$hitstr = $hitstr."$CodonCenters[$i]\n";
    }

    # That's it, dude!
    return ($hitstr,$spaln_pct_id);
    
}








############################################################
#
#  Function: GroupRanges
#
sub GroupRanges
{
    my $starts_ref = shift;
    my $ends_ref   = shift;
    my $max_span   = shift;
    my $revcomp    = shift;

    my @Starts     = @{$starts_ref};
    my @Ends       = @{$ends_ref};
    my $num_ranges = scalar(@Starts);

    # Special case: single range gets returned
    if ($num_ranges == 1) {
	return $Starts[0].'..'.$Ends[0];
    }

    # Find (1.) the total span of the range, and (2.) the pair of sequences
    # with the greatest distance between them.
    my $span = 0;
    my $left_leader  = -1;
    my $right_leader = -1;
    my $top_gap_size = 0;
    for (my $i=0; $i<$num_ranges; $i++) {
	for (my $j=0; $j<$num_ranges; $j++) {
	    if (( $revcomp && $Ends[$i] > $Starts[$j]) ||
		(!$revcomp && $Ends[$i] < $Starts[$j])) {

		$span = Max($span,abs($Starts[$i]-$Ends[$j]));
		my $gap_size = abs($Ends[$i]-$Starts[$j]);
		if ($top_gap_size < $gap_size) {
		    $top_gap_size = $gap_size;
		    $left_leader  = $i;
		    $right_leader = $j;
		}

	    }
	}
    }

    # There's a *very* special case where overlap would lead us to not pick up
    # any hits, in which case we just figure out what the min and max are, and
    # return those.
    if ($top_gap_size == 0) {
	my $min = 2 ** 35;
	my $max = 0;
	for (my $i=0; $i<$num_ranges; $i++) {
	    $min = Min($min,$Starts[$i]);
	    $max = Max($max,$Ends[$i]);
	}
	if ($revcomp) { return $max.'..'.$min; }
	else          { return $min.'..'.$max; }
    }

    # If the span is below our maximum, we're good!
    if ($span < $max_span) {
	return $Starts[$left_leader].'..'.$Ends[$right_leader];
    }

    # Nope, recursion time!  Group sequences according to who they're closest to.
    # Note that this is a little fast 'n' loose, but given the span is as large as
    # it is, we aren't worried about (what I assume would be) minor roughness
    my @LeftStarts;
    my @LeftEnds;
    my @RightStarts;
    my @RightEnds;
    for (my $i=0; $i<$num_ranges; $i++) {
	my $left_dist = abs($Starts[$i] - $Starts[$left_leader]);
	my $right_dist = abs($Starts[$i] - $Starts[$right_leader]);
	if ($left_dist < $right_dist) {
	    push(@LeftStarts,$Starts[$i]);
	    push(@LeftEnds,$Ends[$i]);
	} else {
	    push(@RightStarts,$Starts[$i]);
	    push(@RightEnds,$Ends[$i]);
	}
    }

    # RECURSE
    my $left_ranges_str  = GroupRanges(\@LeftStarts,\@LeftEnds,$max_span,$revcomp);
    my $right_ranges_str = GroupRanges(\@RightStarts,\@RightEnds,$max_span,$revcomp);

    # The very last thing we'll do is give a bit of an extended buffer to each side,
    # by employing the *worst* variable names you've ever seen!
    $left_ranges_str =~ /^(\S+)\.\.(\d+)$/;
    my $left_left  = $1;
    my $left_right = $2;
    $right_ranges_str =~ /^(\d+)\.\.(\S+)$/;
    my $right_left  = $1;
    my $right_right = $2;

    # Given that the current max span is 5Mb, I'll allow 500Kb -- and pretend this is
    # more generally a reasonable computation.
    if ($revcomp) {
	$left_right -= int($max_span/10);
	$right_left += int($max_span/10);
    } else {
	$left_right += int($max_span/10);
	$right_left -= int($max_span/10);
    }
    $left_ranges_str  = $left_left.'..'.$left_right;
    $right_ranges_str = $right_left.'..'.$right_right;

    # And away you go!
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
    my $revcomp  = shift;

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

    $Jumps[$index] =~ /^(\d+)\:(\d+)$/;
    my $relative_pos = $1;
    my $absolute_pos = $2;

    # We've found our closest jump point -- but we need to get particular about things!
    if ($revcomp) { $absolute_pos -= $coord - $relative_pos; }
    else          { $absolute_pos += $coord - $relative_pos; }

    return $absolute_pos;
    
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
	    $nucl_line =~ /^\s+(\d+)/;
	    my $line_nucl_pos = $1;

	    # A catch for the first "skip"
	    if ($recent_skip == -1) {
		$GroupNuclStarts[$num_groups] = $line_nucl_pos;
	    } elsif ($nucl_pos + $recent_skip != $line_nucl_pos) {
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

    # We'll generate a bunch of processes to make things go down sooooo smooth
    $num_cpus = Min($num_cpus,scalar(@OutFileList));
    my $threadID = SpawnProcesses($num_cpus);

    # Off to the file mines with you!
    my $start_file_index =  $threadID    * (scalar(@OutFileList)/$num_cpus);
    my $end_file_index   = ($threadID+1) * (scalar(@OutFileList)/$num_cpus);
    $end_file_index = scalar(@OutFileList) if ($threadID == $num_cpus-1);

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
	    if ($chr ne 'Incomplete' && $chr ne 'Chimeric'
		&& $Chrs{$chr} > $top_chr_count) {
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
		if ($chr eq 'Incomplete' || $chr eq 'Chimeric') {
		    while (length($chr) < 11) {
			$chr = $chr.' ';
		    }
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
	
    }
    
    # Take a breather -- you've earned it!
    if ($threadID) { exit(0); }
    while (wait() != -1) {}
    
}








# EOF











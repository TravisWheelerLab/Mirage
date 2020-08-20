#!/usr/bin/env perl
#
# MultiMSA.pl - Alex Nord - 2016
#
# ABOUT: In the mirage pipeline, this script is used to convert the exon-
#        aware hits produced by 'Quilter' into intra-species gene-family
#        MSAs.  It does this by scanning the nucleotide positions for each
#        hit within a gene family and inserting the corresponding (character
#        and sequence name) pair into a hash.  The keys to this hash, when
#        sorted, represent the columns of the hash.  A little extra work is
#        used to track splice sites.
#
use warnings;
use strict;
use POSIX;
use Cwd;

# YUCKITY YUCK YUCK
sub GetThisDir { my $lib = $0; $lib =~ s/\/MultiMSA\.pl$//; return $lib; }
use lib GetThisDir();
use BureaucracyMirage;
use DisplayProgress;


sub WriteMSAToFile;


sub GetNextExonRange;   # Identify the start and end points of the next exon (nucleotide values)
sub UniteARFSegments;   # When we have adjacent ARF segments, group them by sequence
sub MatrixRecurse;      # Used by UniteARFSegments during clustering
sub MinorClean;         # Clean up obvious minor errors in the alignment
sub QuickAlign;         # Used to resolve indels
sub ResolveExon;        # When an error has been detected in an exon, try to adjust using gap tricks
sub CheckColumnProfile; # Build a profile of an MSA column



if (@ARGV < 1) { die "\n  USAGE:  ./MapsToMSAs.pl [SpeciesMSAs/species/seqs]\n\n"; }


# We should be able to work entirely from the directory with Quilter output
my $dirname = ConfirmDirectory($ARGV[0]);

# Let's figure out what our species' name is
$dirname =~ /\/([^\/]+)\/seqs\/$/;
my $species = $1;

# How many threads do we want to use?
my $ThreadGuide = OpenInputFile($dirname.'Thread-Guide');
my $num_cpus = <$ThreadGuide>;
close($ThreadGuide);
$num_cpus =~ /Num CPUs\: (\d+)/;
$num_cpus = $1;

# Spawn your friends and get to work!
my $threadID = SpawnProcesses($num_cpus);

# Ugh, do we *really* have to reopen that stinky old Thread-Guide?
$ThreadGuide = OpenInputFile($dirname.'Thread-Guide');
my $tg_line = <$ThreadGuide>; # Num CPUs line
while ($tg_line = <$ThreadGuide>) {

    # We only want the finest sequences, hand-selected for this thread
    next if ($tg_line !~ /(\d+) (\S+)/);
    my $req_thread = $1;
    my $gene = $2;

    next if ($req_thread != $threadID);

    my $mapfname = $dirname.$gene.'.quilter.out';
    my $seqfname = $dirname.$gene.'.fa';

    # Read in the sequences, both in their original case-sensitive forms
    # and an all-uppercase format for internal use.
    my @OrigSeqs;
    my @Seqs;
    my @SeqNames;
    my @SeqLengths;
    my %SeqNameToIndex;
    my $num_seqs = -1;
    my $SeqFile = OpenInputFile($seqfname);
    while (my $line = <$SeqFile>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	if ($line =~ /\>(\d+)/) {

	    $SeqNames[++$num_seqs] = $1;
	    $SeqLengths[$num_seqs] = 0;
	    $OrigSeqs[$num_seqs] = '';
	    $Seqs[$num_seqs] = '';

	    # Because things go poorly if this is a zero, we need to
	    # have the index be one higher than the actual index.
	    $SeqNameToIndex{$SeqNames[$num_seqs]} = $num_seqs+1;

	} else {

	    $SeqLengths[$num_seqs] += length($line);
	    $OrigSeqs[$num_seqs] = $OrigSeqs[$num_seqs].$line;
	    $Seqs[$num_seqs] = $Seqs[$num_seqs].uc($line);

	}
	
    }
    close($SeqFile);

    # If there isn't a mapping file, then we write out the whole family
    # as unmapped.
    if (!(-s $mapfname)) {
	my $missf = OpenOutputFile($dirname.$gene.'.misses');
	foreach my $seqname (@SeqNames) {
	    print $missf "$seqname: Unmapped\n";
	}
	close($missf);
	next;
    }

    # Next up, read in the mappings from our mapping file
    my $MapFile = OpenInputFile($mapfname);
    my %Mapping;
    my @Unmapped;

    # What's the canonical chromosome for this gene family?
    my $canon_chr = <$MapFile>;
    $canon_chr =~ /^Canonical Chromosome\: (\S+)/;
    $canon_chr = $1;

    # Are we on the reverse strand?
    my $revcomp = 0;
    $revcomp = 1 if ($canon_chr =~ /\[revcomp\]/);

    # We'll need to be prepared for some sequences to not have usable mappings
    my $num_mapped = 0;

    # Read in those dang mappings!
    while (my $line = <$MapFile>) {

	$line =~ s/\n|\r//g;
	next if (!$line);

	$line =~ /Sequence ID\: (\d+)/;
	my $seqname = $1;
	my $seq_id  = $SeqNameToIndex{$seqname}-1;

	$line = <$MapFile>;
	$line =~ /Map Method \: (\S+)/;
	my $map_method = $1;

	if ($map_method eq 'Unmapped') {
	    push(@Unmapped,$seqname.': Unmapped');
	    next;
	}

	$line = <$MapFile>;
	$line =~ /Chromosome \: (\S+)/;
	my $chr = $1;

	# If we have a match between the chromosomes, we're incorporating this sequence
	# right into our MSA!
	my @SeqChars;
	if ($chr eq $canon_chr) {
	    $SeqNames[$num_mapped] = $SeqNames[$seq_id];
	    $SeqLengths[$num_mapped] = $SeqLengths[$seq_id];
	    $OrigSeqs[$num_mapped] = $OrigSeqs[$seq_id];
	    $Seqs[$num_mapped] = $Seqs[$seq_id];
	    @SeqChars = split(//,$Seqs[$seq_id]);
	    $seq_id = $num_mapped;
	    $num_mapped++;
	} else {
	    push(@Unmapped,$seqname,': Noncanonical Chromosome');
	}

	$line = <$MapFile>;
	$line =~ /Num Exons  \: (\d+)/;
	my $num_exons = $1;

	# Read in the actual mappings! Note that we don't really need
	# to use the first line, since that just gives metadata about
	# the exon.
	my $seq_pos = 0;
	for (my $i=0; $i<$num_exons; $i++) {

	    $line = <$MapFile>;
	    $line = <$MapFile>;
	    next if ($chr ne $canon_chr);

	    foreach my $map_pos (split(/\,/,$line)) {
		my $entry = $seq_id.':'.$SeqChars[$seq_pos++];
		if ($Mapping{$map_pos}) { $Mapping{$map_pos} = $Mapping{$map_pos}.','.$entry; }
		else                    { $Mapping{$map_pos} = $entry;                        }
	    }
	    
	}
	
	
    }
    close($MapFile);

    # Next up, if we had any unmapped sequences, we'll sing about them reverently in
    # the form of a little dirge called file I/O
    if (scalar(@Unmapped)) {
	my $MissFile = OpenOutputFile($dirname.$gene.'.misses');
	foreach my $unmapped_seq (@Unmapped) {
	    print $MissFile "$unmapped_seq\n";
	}
	close($MissFile);
    }

    # Did this entire family miss? That would make me very sad :(
    $num_seqs = $num_mapped;
    next if ($num_seqs == 0);

    # Excellent! Let's kick things off by drafting our spliced MSA!
    # First up, get your list of coordinates ready
    my @MapCoords;
    if ($revcomp) { @MapCoords = sort { $b <=> $a } keys %Mapping; }
    else          { @MapCoords = sort { $a <=> $b } keys %Mapping; }

    # Now we'll actually draft that MSA! We'll also store the mapping coordinates
    # in a top-secret final row, but don't go telling about it
    my @MSA;
    my $msa_len = 0;
    my $last_coord = -1;
    my @ExonStarts;
    my @ExonEnds;
    my $num_exons = -1;
    my %IndelColumns;
    foreach my $coord (@MapCoords) {

	# Is this the start of a new exon?
	if (abs($coord-$last_coord)>4) {
	    $ExonEnds[$num_exons] = $msa_len-1 if ($num_exons >= 0);
	    for (my $i=0; $i<=$num_seqs; $i++) {
		$MSA[$i][$msa_len] = '*';
	    }
	    $msa_len++;
	    $ExonStarts[++$num_exons] = $msa_len;
	}
	$last_coord = $coord;

	# Initialize a list to represent the column
	my @Column;
	for (my $i=0; $i<$num_seqs; $i++) {
	    $Column[$i] = 0;
	}

	# Now we'll check that coordinate to see who's representing this column
	# NOTE that we need to be careful about possible indels
	foreach my $entry (split(/\,/,$Mapping{$coord})) {
	    $entry =~ /(\d+)\:(\S)/;
	    my $seq_id = $1;
	    my $char = $2;
	    if ($Column[$seq_id]) {
		$Column[$seq_id] = $Column[$seq_id].$char;
		$IndelColumns{$msa_len} = 1;
	    } else {
		$Column[$seq_id] = $char;
	    }
	}

	# MSA time, baby!
	for (my $i=0; $i<$num_seqs; $i++) {
	    $MSA[$i][$msa_len] = $Column[$i];
	}
	$MSA[$num_seqs][$msa_len] = $coord;
	$msa_len++;
	
    }
    $ExonEnds[$num_exons] = $msa_len-1;
    $num_exons++;



    # AND FINALLY
    my $outfname = $dirname.$gene.'.afa';
    WriteMSAToFile(\@MSA,\@SeqNames,\@OrigSeqs,$num_seqs,$msa_len,$outfname);
    
}

close($ThreadGuide);

if ($threadID) { exit(0); }
while (wait() != -1) {}


# Now that it's just you 'n' me, Master Thread, let's concatenate all of our
# missed sequences into a single file.
my $Dir = OpenDirectory($dirname);
my @MissFiles;
while (my $fname = readdir($Dir)) {
    if ($fname =~ /\.misses$/) {
	push(@MissFiles,$dirname.$fname);
    }
}
closedir($Dir);

# Oh, how nice it would be not to enter this conditional...
if (scalar(@MissFiles)) {
    my $BigMissFile = OpenOutputFile($dirname.'../'.$species.'.misses');
    foreach my $missfname (sort @MissFiles) {
	my $missf =  OpenInputFile($missfname);
	while (my $line = <$missf>) {
	    $line =~ s/\n|\r//g;
	    next if (!$line);
	    print $BigMissFile "$line\n";
	}
	close($missf);
    }
    close($BigMissFile);
}



# CALLIN' IT
1;








#
# OLD STUFF BELOW
#

# NOTE: This stuff is just so that we can run the new script without getting
#       a billion warnings.
if (0) {
    my ($i,$j,$k,$chr_hits,);
    my @PositionIndex;
    my %MSA;
    my @ProteinLens;
    
    
    my @FinalMSA = ();
    for ($i=0; $i<$chr_hits; $i++) { $FinalMSA[$i][0] = '*'; }
    my $msa_len  = 1;

    # In order to log ARF positions, we need to keep track of where we are in
    # a given sequence.
    my @ARFNameField;
    my @ContentPositions;
    for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) { 
	$ARFNameField[$seq_id] = 0; 
	$ContentPositions[$seq_id] = 0;
    }

    # If our last segment had ID'd a subset of our sequences as non-ARFs,
    # we'll want that labeling to cary through into the next segment.
    my $last_arf_seqs = 'X'; # Since 0 could be a valid last sequence to see an ARF in

    my $pos_index = 0;
    while ($pos_index < scalar(@PositionIndex)) {

	# Figuring out the range of indices in our PositionIndex
	# corresponding to the exon we're currently building.
	my ($next_index,$segment_ref) = GetNextExonRange($pos_index,\@PositionIndex,\%MSA);
	my @SegmentList = @{$segment_ref};

	# We need to know where segments with multiple frames
	# occur adjacent to one another, so that if we need to
	# stitch sequences together (i.e., organize them by sequence) 
	# we can.
	my @MultiFrameSegs;
	my @MultiFrameStarts;
	my @MultiFrameEnds;
	my $num_segs = 0;
	
	# Because we may need to combine ARFs within some segments, we
	# need a temporary place to record the positions of ARFs.
	my @TempARFNameField;
	for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
	    for (my $seg_id = 0; $seg_id < scalar(@SegmentList); $seg_id++) {
		$TempARFNameField[$seq_id][$seg_id] = 0;
	    }
	}
	
	foreach my $segment (@SegmentList) {

	    my @SegmentParts = split(/\,/,$segment);
	    my $start_index  = $SegmentParts[0];
	    my $end_index    = $SegmentParts[1];
	    my $num_frames   = $SegmentParts[2];

	    # If there are multiple frames, but none start/end the protein, we'll do a quick
	    # check to see which frame is the most common and designate that the non-ARF.
	    #
	    # NOTE: This is our ad-hoc selection -- later on, we'll do a check to see if there
	    # are any sequences that begin or end in a DCE region, and if so we'll call them the
	    # ARF (since we suspect they're providing signal for NMD).
	    #
	    my $non_arf_frame = 0;
	    my $non_arf_frame_size = scalar(split(/\,/,$MSA{$PositionIndex[$start_index]}));
	    for (my $frame=1; $frame<$num_frames; $frame++) {
		my $frame_size = scalar(split(/\,/,$MSA{$PositionIndex[$start_index+$frame]}));
		if ($frame_size > $non_arf_frame_size) {
		    $non_arf_frame = $frame;
		    $non_arf_frame_size = $frame_size;
		}
	    }
	    
	    # It's probably easiest to just do this here, rather than on a
	    # conditional
	    my @ContentStarts;
	    for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
		# Note that we do a '+1' for non-CS enumeration
		$ContentStarts[$seq_id] = $ContentPositions[$seq_id]+1;
	    }


	    # If there are multiple frames, we'd better take note
	    if ($num_frames > 1) {
		push(@MultiFrameSegs,$num_segs);
		push(@MultiFrameStarts,$msa_len);
	    }

	    # UGH, this is a naming convention nightmare...
	    my %SeqsToFrames;

	    # I mean, pretty slick, right? (this is how we handle ARFs)
	    for (my $frame=0; $frame<$num_frames; $frame++) {

		# We benefit from knowing which sequences we've seen in this frame
		my @SeqsInFrame;
		for (my $seq_id=0; $seq_id<$chr_hits; $seq_id++) { $SeqsInFrame[$seq_id] = 0; }

		for (my $index=$start_index+$frame; $index<$end_index; $index+=$num_frames) {

		    my $genome_pos   = $PositionIndex[$index];
		    my @SeqsAndChars = split(/\,/,$MSA{$genome_pos});
		    my @Seqs;
		    my @Chars;
		    
		    # If we have a longest entry >1 then we use our quick N-W like
		    # thing to generate an alignment for the whole column.
		    my @Column;
		    for ($i=0; $i<$chr_hits; $i++) { $Column[$i] = '-'; }

		    # We'll need to check the lengths of the characters at
		    # this entry to see if we're in the land of indels
		    my $longest_entry_len = 0;
		    my $longest_entry_seq = -1;
		    for ($i=0; $i<scalar(@SeqsAndChars); $i++) {

			$SeqsAndChars[$i] =~ /(\d+)\:(\S+)/;
			my $seq_id = $1;
			my $chars  = $2;

			$SeqsInFrame[$seq_id] = 1;
			$SeqsToFrames{$seq_id} = $frame + 1; # Need '+1' to avoid missing frame zero usage...

			# If there are multiple characters in this entry, we'll
			# see if there's a more creative solution to this puzzle
			# (namely, pushing all but the leftmost over to the next
			# position).
			if (length($chars) > 1 && $index+$num_frames < $end_index) {

			    $chars =~ /^(\S)(\S)/;
			    my $first_char  = $1;
			    my $second_char = uc($2);

			    # If the second character is what more than half of the
			    # entries at the upcoming position have, then we'll move
			    # our indel content that-a-way
			    my $check_pos    = $PositionIndex[$index+$num_frames];
			    my $check_entry  = uc($MSA{$check_pos});
			    my $num_elements = scalar(split(/\,/,$check_entry));
			    my $num_matches  = scalar(split(/\:$second_char/,$check_entry));
			    if ($num_matches > $num_elements/2) {

				$chars =~ /^\S(\S+)$/;
				my $next_entry = $1;

				# We'll need to be careful, in case this sequence
				# already has an entry at that position.
				if ($MSA{$check_pos} =~ /^$seq_id\:/) {
				    $MSA{$check_pos} =~ s/^$seq_id\:/$seq_id\:$next_entry/;
				} elsif ($MSA{$check_pos} =~ /\,$seq_id\:/) {
				    $MSA{$check_pos} =~ s/\,$seq_id\:/\,$seq_id\:$next_entry/;
				} else {
				    $MSA{$check_pos} = $MSA{$check_pos}.','.$seq_id.':'.$next_entry;
				}

				$chars = $first_char;

			    }
			}

			# We want to know what the longest entry was in case we need to
			# do some pseudo-alignment to make us happy with this set of
			# columns.
			if (length($chars) > $longest_entry_len) {
			    $longest_entry_len = length($chars);
			    $longest_entry_seq = $seq_id;
			}

			# Add this set of characters to the column
			$Column[$seq_id] = $chars;
			$ContentPositions[$seq_id] += length($chars);

		    }

		    # If we had a column with more than one character in it, then
		    # we'll need to do a quick little bit of alignment work.
		    if ($longest_entry_len > 1) {
			my $column_ref = QuickAlign(\@Column,$longest_entry_seq,$longest_entry_len,$chr_hits);
			@Column = @{$column_ref};
		    }

		    # And, finally, load up the MSA!
		    for ($i=0; $i<$chr_hits; $i++) {
			my @RowItems = split(//,$Column[$i]);
			for ($j=0; $j<$longest_entry_len; $j++) {
			    $FinalMSA[$i][$msa_len+$j] = $RowItems[$j];
			}
		    }

		    # Finally (even more finally than the last thing we did finally),
		    # increment the length of the MSA.
		    $msa_len += $longest_entry_len;

		}

	    }

	    # We want to prioritize designating exons that start/end proteins as
	    # ARFs, so we'll first do a check to see whether there are multiple
	    # reading frames and if one starts off the sequence.
	    
	    # If we just wrapped up an ARF frame, then we'll record how much 
	    # each of the sequences in this frame moved.
	    
	    # NEW APPROACH
	    if ($num_frames > 1) {

		# If we encountered an ARF in the last segment, we'll want to
		# ensure that a contiguous DCE uses a consistent labeling.
		my $precleared = 0;
		if ($last_arf_seqs !~ /X/) {

		    # Turn this into a list, and see if there are any frames without
		    # a member from this list...
		    my @ClearFrames;
		    for (my $frame=0; $frame<$num_frames; $frame++) { $ClearFrames[$frame] = 1; }
		    foreach my $seq_id (split(/\,/,$last_arf_seqs)) {
			if ($SeqsToFrames{$seq_id}) {
			    $ClearFrames[$SeqsToFrames{$seq_id}-1] = 0;
			}
		    }

		    # If our current non-ARF frame looks good, we're solid
		    if ($ClearFrames[$non_arf_frame]) {
			$precleared = 1;
		    } else {
			# Just pick one of the clear frames (if there is one!)
			for (my $frame=0; $frame<$num_frames; $frame++) {
			    if ($ClearFrames[$frame]) {
				$non_arf_frame = $frame;
				$precleared = 1;
				last;
			    }
			}
		    }

		} 

		# If we couldn't get a clear non-ARF from previous knowledge, we'll check
		# to see which of these sequences look non-ARF-y
		if (!$precleared) {

		    # First, we'll organize sequences into these three categories
		    my @TrueInteriors;
		    my @CoversStartExon;
		    my @CoversEndExon;
		    
		    foreach my $seq_id (keys %SeqsToFrames) {
			# Recall that these are from '1..seq_len'
			my $start_pos = $ContentStarts[$seq_id];
			my $end_pos   = $ContentPositions[$seq_id];
			if ($start_pos > 1 && $end_pos < $ProteinLens[$seq_id]) {
			    push(@TrueInteriors,$seq_id);
			} else {
			    # We'll allow sequences to be both starting and ending, since
			    # that's conceivable.
			    if ($start_pos == 1) { push(@CoversStartExon,$seq_id); }
			    if ($end_pos == $ProteinLens[$seq_id]) { push(@CoversEndExon,$seq_id); }
			}
		    }
		    
		    # If we don't have any DCEs that start or end their sequences, we defer back
		    # to our silly 'non_arf_frame' designation.
		    if (scalar(@CoversStartExon) + scalar(@CoversEndExon) > 0) {
			
			# Check if we have a frame that doesn't touch any starting / ending exons
			my @ClearFrames;
			for (my $frame=0; $frame<$num_frames; $frame++) { $ClearFrames[$frame] = 1; }
			
			foreach my $seq_id (@CoversStartExon) { $ClearFrames[$SeqsToFrames{$seq_id}-1] = 0; }
			foreach my $seq_id (@CoversEndExon)   { $ClearFrames[$SeqsToFrames{$seq_id}-1] = 0; }
			
			# If our original non-ARF-frame is clear, stick with it.  Otherwise, pick any safe frame.
			if ($ClearFrames[$non_arf_frame] == 0) {
			    my $safe_non_arf = -1;
			    for (my $frame=0; $frame<$num_frames; $frame++) {
				if ($ClearFrames[$frame]) {
				    $safe_non_arf = $frame;
				}
			    }
			    if ($safe_non_arf != -1) { $non_arf_frame = $safe_non_arf; }
			}
			
		    }
		}
	    }

	    # This is always a good place to start forgetting the last segment
	    $last_arf_seqs = 'X';
	    
	    # Now we can actually do our labeling!
	    foreach my $seq_id (keys %SeqsToFrames) {
		if ($SeqsToFrames{$seq_id}-1 != $non_arf_frame) {
		    # To avoid weirdness associated with SPALN sticking single amino acids in places,
		    # we'll skip annotating single amino ARFs
		    if ($ContentStarts[$seq_id] < $ContentPositions[$seq_id]) {
			$TempARFNameField[$seq_id][$num_segs] = $ContentStarts[$seq_id].'..'.$ContentPositions[$seq_id];
			if ($last_arf_seqs =~ 'X') { $last_arf_seqs = $seq_id;                    }
			else                       { $last_arf_seqs = $last_arf_seqs.','.$seq_id; }
		    }
		} 
	    }


	    # Cool! That wraps up that segment, but if we saw multiple reading frames
	    # then we're going to need to take note of where the sequence ended.
	    # Note that we'll use the [start,end) format for easy loop-writing.
	    if ($num_frames > 1) {
		push(@MultiFrameEnds,$msa_len);
	    }
	    
	    # 'Nother seggo done
	    $num_segs++;

	}

	# Before we call it wraps on this exon, we'll make sure that we don't have any
	# runs of multi-frame segments (i.e., 2 -> 3 -> 2 stuff) to potentially gussy
	# up.
	#
	# NOTE:  While the above describes the original use of this code (where we
	#        conditioned the stuff after computing 'seg_run_end' on
	#        seg_run_len > 1, we've switched to just running it anyways, so that
	#        (hopefully) all ARFs look as visually intuitive as possible.
	#
	my $seg_run_start = 0;
	my $num_multi_frames = scalar(@MultiFrameSegs);
	while ($seg_run_start < $num_multi_frames) {

	    # Find the starting segment and the length of this multi-frame run
	    my $seg_run_len = 1;
	    while ($seg_run_start + $seg_run_len < $num_multi_frames 
		   && $MultiFrameSegs[$seg_run_start + $seg_run_len] == $MultiFrameSegs[$seg_run_start] + $seg_run_len) {
		$seg_run_len++;
	    }
	    
	    # What's the ending segment?
	    my $seg_run_end = $seg_run_start + ($seg_run_len - 1);

	    # In case this isn't the first segment in this exon, we'll try to
	    # keep our first cluster in line with the end of the last segment
	    # (or, if this is the first segment, same-ish thing, but with the
	    # last). 
	    my @LineUpWithSeqs;
	    my $line_up = 0;
	    if ($seg_run_start) {
		$line_up = 1;
		my $last_seg_end = $MultiFrameStarts[$seg_run_start]-1;
		for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
		    if ($FinalMSA[$seq_id][$last_seg_end] =~ /[A-Za-z]/) {
			$LineUpWithSeqs[$seq_id] = 1;
		    } else {
			$LineUpWithSeqs[$seq_id] = 0;
		    }
		}
	    } elsif ($MultiFrameSegs[$seg_run_end]+1 < $num_segs) {
		$line_up = 2;
		my $next_seg_start = $MultiFrameEnds[$seg_run_end]+1;
		for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
		    if ($FinalMSA[$seq_id][$next_seg_start] =~ /[A-Za-z]/) {
			$LineUpWithSeqs[$seq_id] = 1;
		    } else {
			$LineUpWithSeqs[$seq_id] = 0;
		    }
		}
	    }

	    # Correct the MSA
	    my $final_msa_ref = UniteARFSegments(\@FinalMSA,$chr_hits,$MultiFrameStarts[$seg_run_start],$MultiFrameEnds[$seg_run_end],\@LineUpWithSeqs,$line_up);
	    @FinalMSA         = @{$final_msa_ref};
	    
	    # Correct our ARF annotation
	    for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
		my $new_annotation;
		for (my $scan1 = $seg_run_start; $scan1 <= $seg_run_end; $scan1++) {
		    if ($TempARFNameField[$seq_id][$scan1] =~ /^(\d+)\.\./) {
			$new_annotation = $1.'..';
			for (my $scan2 = $seg_run_end; $scan2 >= $seg_run_start; $scan2--) {
			    if ($TempARFNameField[$seq_id][$scan2] =~ /\.\.(\d+)$/) {
				$new_annotation = $new_annotation.$1;
				last;
			    }
			    
			}
			$TempARFNameField[$seq_id][$seg_run_start] = $new_annotation;
			for (my $scan2 = $seg_run_start + 1; $scan2 <= $seg_run_end; $scan2++) {
			    $TempARFNameField[$seq_id][$scan2] = 0;
			}
		    }
		}
	    }
		
	    $seg_run_start += $seg_run_len;

	}
	
	for (my $seq_id=0; $seq_id<$chr_hits; $seq_id++) {
	    for (my $seg_id=0; $seg_id<$num_segs; $seg_id++) {
		if ($TempARFNameField[$seq_id][$seg_id]) {
		    if ($ARFNameField[$seq_id]) {
			$ARFNameField[$seq_id] =~ s/^ARF\:/ARFs\:/;
			$ARFNameField[$seq_id] = $ARFNameField[$seq_id].','.$TempARFNameField[$seq_id][$seg_id];
		    } else {
			$ARFNameField[$seq_id] = 'ARF:'.$TempARFNameField[$seq_id][$seg_id];
		    }
		}
	    }
	}

	# EXON DONE BABY!!!!
	for($i=0; $i<$chr_hits; $i++) { $FinalMSA[$i][$msa_len] = '*'; }
	$msa_len++;

	$pos_index = $next_index;

    }


    #
    #  WOWEE!  All done with the first sketch of the final MSA!
    #


    # Now we can actually assemble our final MSA.  The MSA and 'Disagreements' that
    # we get back from this function are what get stored as the intra-species alignment
    # for this gene family.
    #
    my @Disagreements = ();
    my $final_len = $msa_len;
    if ($chr_hits > 1) {
	my ($FinalMSARef,$DisagreementsRef,$tmp_len) = MinorClean(\@FinalMSA,$chr_hits,$msa_len);
	@FinalMSA = @{$FinalMSARef};
	@Disagreements = @{$DisagreementsRef};
	$final_len = $tmp_len;
    }


    # See if we can connect any ARFs together
    #
    # One thing to consider about this is that it would otherwise preserve clarity
    # of exon boundaries... Probably still better in the long-run to concatenate these
    # multi-exon ARFs, but just a thought.
    #
    for (my $seq_id=0; $seq_id<$chr_hits; $seq_id++) {
	if ($ARFNameField[$seq_id] =~ /\:(.*)$/) {

	    my $arf_range_str = $1;
	    my @ARFRangeList  = split(/\,/,$arf_range_str);
	    $ARFRangeList[0] =~ /^(\d+)\.\.(\d+)$/;

	    my $arf_start_coord = $1;
	    my $arf_end_coord   = $2;
	    my $final_arf_str   = $arf_start_coord.'..';

	    for (my $arf_range=1; $arf_range<scalar(@ARFRangeList); $arf_range++) {
		$ARFRangeList[$arf_range] =~ /^(\d+)\.\.(\d+)$/;
		$arf_start_coord = $1;
		my $next_arf_end_coord = $2;

		if ($arf_start_coord != $arf_end_coord+1) {
		    $final_arf_str = $final_arf_str.$arf_end_coord.','.$arf_start_coord.'..';
		}
		$arf_end_coord = $next_arf_end_coord;
	    }
	    
	    $final_arf_str = $final_arf_str.$arf_end_coord;

	    # We now have the full comma-separated coordinate list, now
	    # we just need the 'ARF(s)' bit...
	    @ARFRangeList = split(/\,/,$final_arf_str);
	    if (scalar(@ARFRangeList) > 1) { $final_arf_str = 'ARFs:'.$final_arf_str; }
	    else                           { $final_arf_str = 'ARF:'.$final_arf_str;  }
	    $ARFNameField[$seq_id] = $final_arf_str;

	}
    }
    
}


# Mad phresh stylez 4 lyfe
1;






#########################
#########################
###                   ###
###    SUBROUTINES    ###
###                   ###
#########################
#########################








################################################################
#
#  Function: ConvertToOrigSeqs
#
sub ConvertToOrigSeqs
{
    my $msa_ref = shift;
    my $orig_seqs_ref = shift;
    my $num_seqs = shift;

    my @MSA = @{$msa_ref};
    my @OrigSeqs = @{$orig_seqs_ref};

    for (my $i=0; $i<$num_seqs; $i++) {

	my @Seq = split(//,$OrigSeqs[$i]);
	my $seq_len = length($OrigSeqs[$i]);

	my $seq_pos = 0;
	my $msa_pos = 0;
	while ($seq_pos < $seq_len) {
	    if ($MSA[$i][$msa_pos] && $MSA[$i][$msa_pos] =~ /[A-Za-z]/) {
		$MSA[$i][$msa_pos] = $Seq[$seq_pos++];
	    }
	    $msa_pos++;
	}

    }

    return \@MSA;
    
}








################################################################
#
#  Function: WriteMSAToFile
#
sub WriteMSAToFile
{
    my $msa_ref = shift;
    my $seqnames_ref = shift;
    my $orig_seqs_ref = shift;
    my $num_seqs = shift;
    my $msa_len = shift;
    my $fname = shift;

    my @MSA = @{$msa_ref};
    my @OrigSeqs = @{$orig_seqs_ref};
    $msa_ref = ConvertToOrigSeqs(\@MSA,\@OrigSeqs,$num_seqs);
    @MSA = @{$msa_ref};
    
    my @SeqNames = @{$seqnames_ref};

    my $outf = OpenOutputFile($fname);
    for (my $i=0; $i<$num_seqs; $i++) {
	print $outf ">$SeqNames[$i]\n";
	for (my $j=0; $j<$msa_len; $j++) {
	    if ($MSA[$i][$j]) { print $outf "$MSA[$i][$j]"; }
	    else              { print $outf '-';            }
	    print $outf "\n" if (($j+1) % 60 == 0);
	}
	print $outf "\n" if ($msa_len % 60);
	print $outf "\n";
    }
    close($outf);
    
}







################################################################
#
# FUNCTION: GetNextExonRange
#
sub GetNextExonRange
{
    my $start_index   = shift;
    my $position_ref  = shift;
    my $msa_hash_ref  = shift;
    my @PositionIndex = @{$position_ref};
    my %MSA           = %{$msa_hash_ref};

    #  NOTE:  We need access to the MSA hash because we need to track
    #         which sequences are involved in which reading frames of
    #         specific segments.  This way, if an indel pushes a sequence
    #         from one reading frame into another we can record the two
    #         segments as separate (and thus avoid weirdness).

    # First, we just want to know how far we go from the start of this
    # exon before we hit an intron
    my $next_index = $start_index;
    my $last_nucl  = $PositionIndex[$start_index];
    $next_index++;
    while ($next_index < scalar(@PositionIndex) && abs($PositionIndex[$next_index]-$last_nucl)<4) {
	$last_nucl = $PositionIndex[$next_index];
	$next_index++;
    }

    # Next, we'll segment the exon based on any ARF-looking stuff.
    # Each segment will have (1.) a start position, (2.) an end position,
    # and (3.) a number of reading frames.
    my @SegmentList;
    my $index = $start_index;
    while ($index < $next_index) {

	my $segment_start = $index;
	
	# In case we end up with a length-1 segment
	$index++;
	if ($index == $next_index) { 
	    push(@SegmentList,$segment_start.','.$index.',1');
	    last;
	}

	my $last_nucl    = $PositionIndex[$index-1];
	my $current_nucl = $PositionIndex[$index];

	# Peeking around super carefully for arfs
	if (abs($last_nucl-$current_nucl)==3) {

	    # If the next entry is 3 nucleotides away, we're in ARF-free territory, baby!
	    $index++;
	    while ($index < $next_index && abs($last_nucl-$current_nucl)==3) {
		$last_nucl = $current_nucl;
		$current_nucl = $PositionIndex[$index];
		$index++;
	    }

	    # In case we appear to have entered into ARF territory, we want to be sure that
	    # we don't include the first character of the ARF with the non-ARF stuff
	    $index-- if (abs($last_nucl-$current_nucl) != 3);

	    push(@SegmentList,$segment_start.','.$index.',1');

	} else {

	    # Uh-oh!  Better start keeping track of which reading frame each sequence
	    # belongs to...
	    my %SeqToFrame;
	    foreach my $seq_char_pair (split(/\,/,$MSA{$PositionIndex[$index]})) {
		$seq_char_pair =~ /^(\d+)\:/;
		my $seq_id = $1;
		$SeqToFrame{$seq_id} = 1;		
	    }

	    my $frame_num     = 2;
	    my $indel_trouble = 0;

	    # Moving up in the MSA world -- note that we're keeping the current_nucl
	    # and last_nucl values a tiny bit desync-ed from the index...
	    $index++;

	    # Could this be... a... TRIPLE ARF?!
	    if ($index+1 < $next_index && abs($PositionIndex[$index+1]-$last_nucl)==3) {

		# IMPOSSIBLE!!!
		while ($index < $next_index && abs($last_nucl-$current_nucl)==1) {

		    # Looking out for trouble
		    foreach my $seq_char_pair (split(/\,/,$MSA{$PositionIndex[$index]})) {
			$seq_char_pair =~ /^(\d+)\:/;
			my $seq_id = $1;
			if (!$SeqToFrame{$seq_id}) {
			    $SeqToFrame{$seq_id} = $frame_num;
			} elsif ($SeqToFrame{$seq_id} != $frame_num) {
			    $indel_trouble = 1;
			}
		    }
		    last if ($indel_trouble);

		    if ($frame_num < 3) { $frame_num++; }
		    else                { $frame_num=1; }

		    $last_nucl = $current_nucl;
		    $current_nucl = $PositionIndex[$index];
		    $index++;

		}

		# In case the ARF doesn't cleanly line up with the end of the
		# exon, we won't count the first character we saw indicating that
		# the reading frame pattern changed.
		$index-- if (abs($last_nucl-$current_nucl) != 1); 

		push(@SegmentList,$segment_start.','.$index.',3');

	    } else {

		# Nope, just the regular ol' kind (not that there's an ARF out there
		# that doesn't blow my socks off).
		#
		# NOTE:  We need to be careful for transitions into triple ARFs by
		#        tracking the pattern of 1->2->1...
		#
		my $expected_jump = 2;
		if (abs($last_nucl-$current_nucl)==1) { $expected_jump = 1; }

		while ($index < $next_index && abs($last_nucl-$current_nucl) == $expected_jump) {

		    # Looking out for trouble
		    foreach my $seq_char_pair (split(/\,/,$MSA{$PositionIndex[$index]})) {
			$seq_char_pair =~ /^(\d+)\:/;
			my $seq_id = $1;
			if (!$SeqToFrame{$seq_id}) {
			    $SeqToFrame{$seq_id} = $frame_num;
			} elsif ($SeqToFrame{$seq_id} != $frame_num) {
			    $indel_trouble = 1;
			}
		    }
		    last if ($indel_trouble);

		    if ($frame_num == 1) { $frame_num=2; }
		    else                 { $frame_num=1; }

		    $last_nucl = $current_nucl;
		    $current_nucl = $PositionIndex[$index];
		    if ($expected_jump == 2) { $expected_jump = 1; }
		    else                     { $expected_jump = 2; }
		    $index++;
		}

		# In case the ARF doesn't cleanly line up with the end of the
		# exon, we won't count the first character we saw indicating that
		# the reading frame pattern changed.
		$index-- if (abs($last_nucl-$current_nucl) != $expected_jump || $indel_trouble);

		push(@SegmentList,$segment_start.','.$index.',2');

	    }
	    
	}
	
    }

    return ($next_index,\@SegmentList);
    
}






################################################################
#
#  FUNCTION:  UniteARFSegments
#
sub UniteARFSegments
{
    my $msa_ref       = shift;
    my $num_seqs      = shift;
    my $segment_start = shift;
    my $segment_end   = shift;
    my $line_up_seqs  = shift;
    my $line_up       = shift;

    my @MSA = @{$msa_ref};

    # If we have information about the sequences that made up the content of
    # the final column of the last segment than we give priority to whichever
    # cluster has the most overlaps with those sequences.
    my @LineUpWithSeqs;
    if ($line_up) { @LineUpWithSeqs = @{$line_up_seqs}; }


    # The way we're going to figure out which sequences need to
    # be grouped together is by generating a CoOccurrence matrix
    # that records which sequences share at least one column,
    # along with an array tracking which sequences don't have
    # any non-gap content.

    
    # Initialization - note that we're only initializing one half,
    # so that for [i][j] we only care about situations where j < i
    my @Cooccur;
    my @NotAllGaps;
    for (my $i=0; $i<$num_seqs; $i++) { # Need 0 for NotAllGaps
	$NotAllGaps[$i] = 0;
	for (my $j=0; $j<$i; $j++) {
	    $Cooccur[$i][$j] = 0;
	}
    }

    # Party time (Scanning)!
    for (my $j=$segment_start; $j<$segment_end; $j++) {

	# Build up our column
	my @Column;
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($MSA[$i][$j] =~ /[A-Za-z]/) {
		$NotAllGaps[$i] = 1;
		push(@Column,$i);
	    }
	}

	# If we saw multiple sequences in this column, record all sequences
	# that co-occur
	if (scalar(@Column) > 1) {
	    for (my $j=0; $j<scalar(@Column)-1; $j++) {
		for (my $i=$j+1; $i<scalar(@Column); $i++) {
		    $Cooccur[$Column[$i]][$Column[$j]] = 1 
		}
	    }
	}

    }

    # Prepare for clustering, you dastardly sequences!
    my @Clustered;
    my @CheckEm;
    my $num_to_check = 0;
    for (my $i=0; $i<$num_seqs; $i++) {
	$Clustered[$i] = 0;
	if ($NotAllGaps[$i]) {
	    $CheckEm[$num_to_check] = $i;
	    $num_to_check++;
	}
    }

    # And, cluster!
    #
    # We start with the highest number because this lets us scan a
    # row from 0 to i-1 without having to keep track of i explicitly
    my $num_clusters = 0;
    for (my $check=$num_to_check-1; $check>=0; $check--) {

	# Skip anyone we've already clustered
	next if ($Clustered[$CheckEm[$check]]);

	# We'll need to cluster at least this sequence
	$num_clusters++;

	# Figure out the set of sequences that belong to this cluster
	my ($cluster_list,) = MatrixRecurse(\@Cooccur,$CheckEm[$check],$num_seqs);
	foreach my $i (split(/\,/,$cluster_list)) {
	    $Clustered[$i] = $num_clusters;
	}

    }


    # Terrific!  Now we know (a.) which sequences ought to be clustered
    # together, and (b.) which sequences we don't really care about (i.e.,
    # are all-gaps).
    #
    # Now we can pass through the "Clustered" array num_clusters times and
    # construct our replacement alignment of the segments we're interested
    # in.

    
    # Before we fully jump in, however, we'll check to see if we want to
    # privilege any of the clusters over the others
    my $privileged_cluster  = 0;
    if ($line_up) {
	my $privileged_overlaps = 0;
	for (my $cluster = 1; $cluster <= $num_clusters; $cluster++) {
	    my $cluster_score = 0;
	    for (my $i=0; $i<$num_seqs; $i++) {
		if ($Clustered[$i] == $cluster && $LineUpWithSeqs[$i]) { 
		    $cluster_score++;
		}
	    }
	    if ($cluster_score > $privileged_overlaps) {
		$privileged_cluster  = $cluster;
		$privileged_overlaps = $cluster_score;
	    }
	}
    }

    # Construct a cluster order, taking privilege into account
    my @ClusterOrder;
    if ($line_up == 1) { push(@ClusterOrder,$privileged_cluster); }
    for (my $cluster = 1; $cluster <= $num_clusters; $cluster++) {
	next if ($cluster == $privileged_cluster);
	push(@ClusterOrder,$cluster);
    }
    if ($line_up == 2) { push(@ClusterOrder,$privileged_cluster);}

    # Initialize our replacement alignment matrix
    my @Replacement;
    my $replacement_len = 0;
    foreach my $cluster (@ClusterOrder) {

	# Which sequences are in this cluster?
	my @SeqsInCluster;
	my $num_seqs_in_cluster = 0;
	for (my $i=0; $i<$num_seqs; $i++) {
	    if ($Clustered[$i] == $cluster) {
		$SeqsInCluster[$num_seqs_in_cluster] = $i;
		$num_seqs_in_cluster++;
	    }
	}

	# Now we run through the relevant part of the MSA,
	# and any column with one of the sequences in this
	# cluster gets pushed onto our Replacement
	for (my $j=$segment_start; $j<$segment_end; $j++) {

	    # Does this column belong to this cluster?
	    my $add_this_col = 0;
	    for (my $i=0; $i<$num_seqs_in_cluster; $i++) {
		if ($MSA[$SeqsInCluster[$i]][$j] =~ /[A-Za-z]/) {
		    $add_this_col = 1;
		    last;
		}
	    }

	    # If it belongs, stick it on, baby
	    if ($add_this_col) {
		for (my $i=0; $i<$num_seqs; $i++) {
		    $Replacement[$i][$replacement_len] = $MSA[$i][$j];
		}
		$replacement_len++;
	    }

	}

    }

    # Hooray!  Now, for our final step, we replace the relevant part of the
    # MSA with our replacement
    my $pos = 0;
    for (my $j=$segment_start; $j<$segment_end; $j++) {
	for (my $i=0; $i<$num_seqs; $i++) {
	    $MSA[$i][$j] = $Replacement[$i][$pos];
	} 
	$pos++;
    }

    # Return that sucker!
    return \@MSA;

}





################################################################
#
#  FUNCTION:  MatrixRecurse
#
sub MatrixRecurse
{
    
    my $matrix_ref = shift;
    my $seq_id     = shift;
    my $max_dim    = shift;

    return ('0',$matrix_ref) if ($seq_id == 0);

    my @M = @{$matrix_ref};

    my %RecurseOn;
    
    for (my $j=$seq_id-1; $j>=0; $j--) {
	if ($M[$seq_id][$j]) {
	    $M[$seq_id][$j] = 0;
	    $RecurseOn{$j}  = 1;
	}
    }

    for (my $i=1; $i<$max_dim; $i++) {
	if ($M[$i][$seq_id]) {
	    $M[$i][$seq_id] = 0;
	    $RecurseOn{$i}  = 1;
	}
    }
    
    my $return_string = "$seq_id";
    foreach my $next_seq_id (sort {$b <=> $a} keys %RecurseOn) {
	my ($next_string,$new_matrix_ref) = MatrixRecurse(\@M,$next_seq_id,$max_dim);
	$return_string = $return_string.','.$next_string;
	@M = @{$new_matrix_ref};
    }
    
    return ($return_string,\@M);

}





################################################################
#
#  FUNCTION:  MinorClean
#
sub MinorClean
{
    my $msa_ref  = shift;
    my $num_seqs = shift;
    my $msa_len  = shift;

    my @MSA = @{$msa_ref};

    # It'll be useful for later components of 'cleanup' to know where
    # our introns are
    my @IntronPositions;

    # We'll begin by doing a forward sweep looking for disagreements,
    # and seeing if we can resolve them by shifting characters across
    # gaps.
    # Note that we won't create new columns -- just consider shifting
    # around pieces that have obvious landing zones.
    my $exon_start;
    my $exon_end;
    my $exon_len;
    my $exon_conflicts = 0;
    my $exon_agreements = 0;
    my @TopChars;
    for (my $j=0; $j<$msa_len; $j++) {
	
	# If we hit an intron marker we look back and consider fixing this exon
	if ($MSA[0][$j] eq '*') {

	    push(@IntronPositions,$j);
	    $exon_end = $j;

	    # A somewhat arbitrary selection, but whatever...
	    # NOTE:  Used to be 'if ($exon_conflicts > 4) {' if that's important...
	    if ($exon_conflicts) {

		# We know this exon is having some issues -- let's see if we can
		# help out!
		$msa_ref = ResolveExon(\@MSA,$num_seqs,$exon_start,$exon_end);
		@MSA = @{$msa_ref};

	    }
	    
	    # Rest our stats for the next round
	    $exon_start = $j+1;
	    $exon_conflicts = 0;
	    $exon_agreements = 0;
	    $exon_len = 0;	    
	    next;
	}
	
	# Otherwise, tally whether we have an agreement or disagreement and
	# record the top character
	my ($num_chars,$top_char,$top_char_count) = CheckColumnProfile(\@MSA,$num_seqs,$j,-1,0);
	if ($num_chars > 1) { $exon_conflicts++;  }
	else                { $exon_agreements++; }
	$TopChars[$exon_len] = $top_char;
	$exon_len++;

    }


    # FINALLY do a pass to see which columns are still disagreeing
    my @Disagreements;
    my $final_len = 0;
    for (my $j=0; $j<$msa_len; $j++) {
	my ($num_chars,$top_char,$top_char_count) = CheckColumnProfile(\@MSA,$num_seqs,$j,-1,0);
	if ($num_chars) {
	    for (my $i=0; $i<$num_seqs; $i++) { $MSA[$i][$final_len] = $MSA[$i][$j]; }
	    $final_len++;
	}
	push(@Disagreements,$final_len) if ($num_chars > 1);
    }

    return(\@MSA,\@Disagreements,$final_len);

}





################################################################
#
# FUNCTION: ResolveExon
#
sub ResolveExon
{
    my $msa_ref  = shift;
    my $num_seqs = shift;
    my $start    = shift;
    my $end      = shift;

    my @MSA = @{$msa_ref};

    # For each sequence, figure out if it's having trouble,
    # and (if it is) see whether the run of characters having
    # the trouble are encountering a problematic gap.
    for (my $i=0; $i<$num_seqs; $i++) {

	# We need to recompute the profile for this exon without
	# this sequence included
	my @TopChar;
	for (my $j=$start; $j<$end; $j++) {
	    my ($num_chars,$top_char,$top_char_count) = CheckColumnProfile(\@MSA,$num_seqs,$j,$i,0);
	    push(@TopChar,$top_char); 
	}

	# Tracking regions within the exon
	my $region_start = 0;
	my $region_end   = 0;
	my $num_probs  = 0;

	#
	#  NOTE:  We'll start things off with a smarter approach where we check
	#    the ends of gaps and try to make long jumps to accomodate our mismatched
	#    characters.  If that doesn't work, we'll check if there are any single-
	#    position shifts (into gaps) that will cover our booties.
	#

	###############
	#             #
	#  BIG JUMPS  #
	#             #
	###############
	
	#
	# We'll start with a forward pass...
	#
	for (my $j=0; $j<$end-$start; $j++) {

	    # How we access the MSA
	    my $pos = $start+$j;

	    # Is this a gap character?
	    if ($MSA[$i][$pos] eq '-' || $j == ($end-$start)-1) {

		$region_end = $j;
		$region_end-- unless ($MSA[$i][$pos] ne '-');

		# If we saw any problems, we'll check whether we can resolve
		# these problems by moving characters to either side of a gap
		if ($num_probs) {

		    my $back_check = $region_start-1;
		    $back_check-- while ($back_check >= 0 && $MSA[$i][$start+$back_check] eq '-');
		    $back_check++; # Bringing us back from the brink

		    if ($back_check != $region_start) {

			# We recompute the TopChar as we move along, since it's liable to change
			my ($num_chars,$top_char,$top_char_count) 
			    = CheckColumnProfile(\@MSA,$num_seqs,$start+$back_check,$i,$MSA[$i][$start+$region_start]);
			$TopChar[$back_check] = $top_char;

			while ($TopChar[$back_check] eq uc($MSA[$i][$start+$region_start]) && $region_start <= $region_end) {

			    $MSA[$i][$start+$back_check] = $MSA[$i][$start+$region_start];
			    $MSA[$i][$start+$region_start] = '-';
			    $region_start++;
			    $back_check++;

			    ($num_chars,$top_char,$top_char_count) 
				= CheckColumnProfile(\@MSA,$num_seqs,$start+$back_check,$i,$MSA[$i][$start+$region_start]);
			    $TopChar[$back_check] = $top_char;

			}
		    }

		}
		
		# Reset the 'region'
		$region_start = $j+1;
		$num_probs    = 0;

	    } elsif (uc($MSA[$i][$pos]) ne $TopChar[$j]) {
		$num_probs++;
	    }

	}

	#
	# ... followed by a reverse pass
	#
	# Note that region_start and region_end follow the logic of rev. comp.,
	# in all of its obnoxiousness.
	#
	for (my $j=($end-$start)-1; $j>=0; $j--) {

	    # How we access the MSA
	    my $pos = $start+$j;

	    # Is this a gap character?
	    if ($MSA[$i][$pos] eq '-' || $j==0) {

		$region_start = $j;
		$region_start++ unless ($MSA[$i][$pos] ne '-');

		# If we saw any problems, we'll check whether we can resolve
		# these problems by moving characters to either side of a gap
		if ($num_probs) {

		    my $fwd_check = $region_end+1;
		    $fwd_check++ while ($fwd_check <= $end-$start && $MSA[$i][$start+$fwd_check] eq '-');
		    $fwd_check--; # Bringing us back from the brink
		    
		    if ($fwd_check != $region_end) {
			
			# We recompute the TopChar as we move along, since it's liable to change
			my ($num_chars,$top_char,$top_char_count)
			    = CheckColumnProfile(\@MSA,$num_seqs,$start+$fwd_check,$i,$MSA[$i][$start+$region_end]);
			$TopChar[$fwd_check] = $top_char;
			    
			while ($TopChar[$fwd_check] eq uc($MSA[$i][$start+$region_end]) && $region_start <= $region_end) {
				
			    $MSA[$i][$start+$fwd_check] = $MSA[$i][$start+$region_end];
			    $MSA[$i][$start+$region_end] = '-';
			    $region_end--;
			    $fwd_check--;

			    ($num_chars,$top_char,$top_char_count) 
				= CheckColumnProfile(\@MSA,$num_seqs,$start+$fwd_check,$i,$MSA[$i][$start+$region_end]);
			    $TopChar[$fwd_check] = $top_char;

			}
		    }

		}
		
		# Reset the 'region'
		$region_end = $j-1;
		$num_probs  = 0;

	    } elsif (uc($MSA[$i][$pos]) ne $TopChar[$j]) {
		$num_probs++;
	    }

	}

    }

    ################
    #              #
    #  LIL SCOOTS  #
    #              #
    ################

    #
    #  Here, what we'll do is just go sequence-by-sequence, and if
    #  a sequence has a run of mismatches with any other sequences,
    #  then we'll see if that run of mismatches can be resolved by a
    #  single shift of the mismatching block to the right or left
    #  into a gap (so that it now matches to *some* sequence).
    #
    #  NOTE:  This is, theoretically, a bad approach, but in practice
    #  will cover a species of bug where we actually need to shift the
    #  majority sequence to line up with the minority sequence.
    #
    #  EXAMPLE:
    #              ABCD        ABCD
    #              BCD-        -BCD
    #              BCD-   ==>  -BCD
    #              BCD-        -BCD
    #              BCD-        -BCD
    #


    # We only need to do one pass -- find a run of mismatches, check for
    # single moves that would resolve the whole run, and make it do the
    # thing.

    for (my $i=0; $i<$num_seqs; $i++) {

	my $j=$start;
	my $mismatch_run_start = 0;

	while ($j < $end) {

	    my $mismatch = 0;
	    if ($MSA[$i][$j] ne '-') {
		for (my $k=0; $k<$num_seqs; $k++) {
		    $mismatch = 1 if ($MSA[$k][$j] ne '-' && lc($MSA[$k][$j]) ne lc($MSA[$i][$j]));
		    last if ($mismatch);
		}
	    }
	    
	    if ($mismatch_run_start == 0 && $mismatch) {
		$mismatch_run_start = $j;	
	    } 

	    # Note:  Need an extra little catch, just in case the mismatch extends
	    # to the very end of the exon.
	    if (($mismatch_run_start && $mismatch == 0) || ($mismatch_run_start && $j == $end-1)) {
		
		# OOOWEE, we got ourselves a mismatch run over here!

		# Look backwards
		my $moved_em_back = 0;
		if ($MSA[$i][$mismatch_run_start-1] eq '-') {
		    
		    for (my $x=0; $x<$num_seqs; $x++) {
			
			next if ($x == $i);
			
			my $match = 1;
			for (my $y=$mismatch_run_start; $y<$j; $y++) {
			    if (uc($MSA[$x][$y-1]) ne uc($MSA[$i][$y])) {
				$match = 0;
				last;
			    }
			}
			
			# OOOOOOOOWEEEEEEEE! We're on the move!
			if ($match) {
			    for (my $y=$mismatch_run_start; $y<$j; $y++) {
				$MSA[$i][$y-1] = $MSA[$i][$y];
			    }
			    $MSA[$i][$j-1] = '-';
			    $moved_em_back = 1;
			    last;
			}
			
		    }
		    
		}
		# Stop looking backwards
		
		# And start looking forwards
		if (!$moved_em_back && $MSA[$i][$j] eq '-') {

		    for (my $x=0; $x<$num_seqs; $x++) {
			
			next if ($x == $i);
			
			my $match = 1;
			for (my $y=$mismatch_run_start; $y<$j; $y++) {
			    if (uc($MSA[$x][$y+1]) ne uc($MSA[$i][$y])) {
				$match = 0;
				last;
			    }
			}
			
			# OOOOOOOOWEEEEEEEE! We're on the move!
			if ($match) {
			    for (my $y=$j; $y>$mismatch_run_start; $y--) {
				$MSA[$i][$y] = $MSA[$i][$y-1];
			    }
			    $MSA[$i][$mismatch_run_start] = '-';
			    last;
			}
			
		    }
		    
		}
		# Jeez, stop looking already
		
	    }

	    $mismatch_run_start = 0 if ($mismatch == 0);
	    $j++;

	}

    }
    
    
    # Pass back the (hopefully) resolved MSA
    return \@MSA;

}





################################################################
#
# FUNCTION: QuickAlign
#
sub QuickAlign
{
    my ($i,$j,$k,$x,$y);
    
    my $match_score    = 3;
    my $mismatch_score = 0;
    my $gap_cost       = -1;

    my $ColumnRef = shift;
    my @Column    = @{$ColumnRef};
    
    my $max_col_ID   = shift;
    my $max_col_size = shift; # <-- This shouldn't ever be larger than a handful of bases

    my $num_seqs = shift;


    # We'll want to have the largest multi-character entry
    # on hand.
    #
    my @BigSeq = split('',$Column[$max_col_ID]);
    
    
    # We can initialize the NW matrix here and just
    # recycle the space as needed.
    #
    my @NWMatrix = ();
    $NWMatrix[0][0] = 0;
    my $all_gaps = '';    
    for ($i=0; $i<=$max_col_size; $i++) {

	# Since we'll be counting one over the max entry length
	# we don't add when i==0.
	if ($i) { $all_gaps = $all_gaps.'-'; }

	# Set the edges, leave the rest to the actual run
	$NWMatrix[0][$i] = $i*$gap_cost;
	$NWMatrix[$i][0] = $i*$gap_cost;
	
    }


    # Now we can just cruise along, trying our darndest
    # to find something of an alignment...
    #
    my @FinalColumn = ();
    for ($i=0; $i<$num_seqs; $i++) {
	
	if ($Column[$i] eq '-') { #------------------------# Just a gapper
	    
	    $FinalColumn[$i] = $all_gaps;
	    
	} elsif (length($Column[$i]) == length($Column[$max_col_ID])) { #--------------------# Max-length seq
	    
	    $FinalColumn[$i] = $Column[$i];
	    
	} else { #-----------------------------------------# Pseudo Needleman-Wunsch
	    
	    my @CompSeq = split('',$Column[$i]);
	    my $complen = @CompSeq;

	    # We only allow matches and HORIZONTAL motion -- gaps in CompSeq
	    my $score;
	    for ($x=0; $x<@BigSeq; $x++) {

		for ($y=0; $y<@CompSeq && $y<=$x; $y++) {

		    if (lc($BigSeq[$x]) eq lc($CompSeq[$y])) { $score = $match_score;    }
		    else                                     { $score = $mismatch_score; }

		    $NWMatrix[$x+1][$y+1] = MAX($NWMatrix[$x][$y]+$score,
						$NWMatrix[$x+1][$y]+$score+$gap_cost);
		    
		}
	    }

	    # Traceback -- here's where we get the aligned sequence
	    my $align_seq = '';
	    $x = $max_col_size-1;
	    $y = $complen-1;

	    while ($x && $x > $y) {

		if (lc($BigSeq[$x]) eq lc($CompSeq[$y])) { $score = $match_score;    }
		else                                     { $score = $mismatch_score; }
		
		if ($NWMatrix[$x+1][$y+1] == $NWMatrix[$x][$y]+$score) {
		    $align_seq = $CompSeq[$y].$align_seq;
		    $y--;
		} else {
		    $align_seq = '-'.$align_seq;
		}

		$x--;

	    }

	    # Walking along the diagonal to (0,0), if applicable
	    #
	    while ($x == $y && $x >= 0) {
		$align_seq = $CompSeq[$y].$align_seq;
		$x--;
		$y--;
	    }

	    # Walking along the edge to (0,0), if applicable
	    #
	    while ($x >= 0) {
		$align_seq = '-'.$align_seq;
		$x--;
	    }

	    # Store that aligned seq!
	    #
	    $FinalColumn[$i] = $align_seq;
	
	}
	
    }

    
    # Radical! Now we can pass back reference to the adjusted column
    #
    return \@FinalColumn;
	
}







################################################################
#
#  FUNCTION:  CheckColumnProfile
#
sub CheckColumnProfile
{
    my $msa_ref    = shift;
    my $num_seqs   = shift;
    my $column_num = shift;
    my $ignore_id  = shift;
    my $preference = shift; # If there's a tie and we have a preferred winner

    my @MSA = @{$msa_ref};
    my %Column;
    for (my $i=0; $i<$num_seqs; $i++) {
	next if ($i == $ignore_id);
	my $char = uc($MSA[$i][$column_num]);
	if ($char ne '-') {
	    if ($Column{$char}) { $Column{$char}++; }
	    else                { $Column{$char}=1; }
	}
    }

    # What was the most popular character?
    my $top_char = 'x';
    my $num_chars = 0;
    my $top_char_count = 0;
    foreach my $char (keys %Column) {
	$num_chars++;
	my $char_count = $Column{$char};
	if ($char_count > $top_char_count) {
	    $top_char = $char;
	    $top_char_count = $char_count;
	} elsif ($preference && $char_count == $top_char_count && $char eq $preference) {
	    $top_char = $char;
	}
    }

    # Return how many characters there were in this column,
    # what the most popular character was, and how many
    # times we saw the most popular character
    return ($num_chars,$top_char,$top_char_count);

}







#######################  END OF FILE  ##########################









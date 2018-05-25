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


sub MIN;                # Get the minimum of two values
sub MAX;                # Get the maximum of two values
sub GetNextExonRange;   # Identify the start and end points of the next exon (nucleotide values)
sub MinorClean;         # Clean up obvious minor errors in the alignment
sub QuickAlign;         # Used to resolve indels
sub ResolveExon;        # When an error has been detected in an exon, try to adjust using gap tricks
sub CheckColumnProfile; # Build a profile of an MSA column


if (@ARGV < 2) { die "\n  USAGE:  perl  MultiMSA.pl  <Quilter output>  <Protein database>  [Opt.s]\n\n"; }

my ($i,$j,$k);

# Marks the ends of exons
my $specialChar = '*';

my $FinalMisses = 'MultiMSA.misses.out';
my $printConsensus = 0;
my $verbose  = 0;
my $specific = 0;
my $outputfolder = 'multilignments';
my $CPUs = 2;
my $canonical_species;
my $mask_arfs = 1;

$i = 2;
while ($i < @ARGV) {
    if ($ARGV[$i] eq '-species') {
	$i++;
	$specific = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-folder') {
	$i++;
	$outputfolder = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-misses') {
	$i++;
	$FinalMisses = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-stack-arfs') {
	$mask_arfs = 0;
    } elsif ($ARGV[$i] eq '-v') {
	$verbose = 1;
    } elsif ($ARGV[$i] eq '-n') {
	$i++;
	$CPUs = int($ARGV[$i]);
	if ($CPUs <= 0) {
	    print "\n  Unsupported number of CPUs requested ($CPUs)\n";
	    print "  Reverting to default (2)\n\n";
	    $CPUs = 2;
	}
    } else {
	print "\n  Unrecognized option '$ARGV[$i]' ignored\n\n";
    }
    $i++;
}

my %Hits;
my %Nums;

if (-e $outputfolder) {
    $i = 2;
    while (-e $outputfolder.'_'.$i) {
	$i++;
    }
    $outputfolder = $outputfolder.'_'.$i;
}
system("mkdir $outputfolder");

open(my $infile,'<',$ARGV[0]) || die "\n  Failed to open input file '$ARGV[0]'\n\n";
my %ChrsByFam;
while (!eof($infile)) {

    my $line = <$infile>;

    if ($line =~ /Isoform ID/) {

	my $ID_line       = $line;
	my $method_line   = <$infile>;
	my $chr_name_line = <$infile>;
	my $pos_line      = <$infile>;
    
	# (1.) Strip down ID line into components 
	$ID_line =~ s/\r|\n//g;
	if (!($ID_line =~ /Isoform ID \: ([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|(\S*\|)?([^\|\s]+)\s*$/)) {
	    die "\n  ERROR:  Unexpected ID line:  $ID_line\n\n";
	}
	
	# NOTE: "GN:" ARTIFACT EXCISION
	my $db_entry  = $1.'|'.$2.'|'.$3.'|'.$4.'|';
	$db_entry = $db_entry.$5 if ($5);
	$db_entry = $db_entry.$6;
	my $orig_group = $6;
	my $group_id   = lc($orig_group);
	
	# Record the species
	$canonical_species = lc($3) if (!$canonical_species);
	
	# (2.) Rip out chromosome name
	$chr_name_line =~ s/\n|\r//g;
	if (!($chr_name_line =~ /Chromosome \: ([\w\-\,\.\/\(\)\'\[\]]+)/)) {
	    die "\n  ERROR:  Unexpected Chromosome line:  $chr_name_line\n\n";
	}
	my $chr_name = $1;
	
	# (3.) Rip out how we matched to this chromosome
	$pos_line =~ s/\n|\r//g;
	if (!($pos_line =~ /Match Pos\.s: (.*+)/)) {
	    die "\n  ERROR:  Unexpected Match Pos.s line:  $pos_line\n\n";
	}
	$pos_line = $1;
	
	# Stuff our new entry down the chimney
	my $hash_entry = $chr_name.'#'.$db_entry.'#'.$pos_line;
	if ($Hits{$group_id}) {
	    $hash_entry = $Hits{$group_id}.'&'.$hash_entry;
	    $Nums{$group_id}++;
	} else {
	    $Nums{$group_id}=1;
	}
	$Hits{$group_id} = $hash_entry;

	if ($ChrsByFam{$group_id}) { $ChrsByFam{$group_id} = $ChrsByFam{$group_id}.'|'.$chr_name; }
	else                       { $ChrsByFam{$group_id} = $chr_name;                           }
	
    }
    
}
close($infile);


# First progress message
my $progmsg;
if ($canonical_species) {
    $progmsg = '  MultiMSA.pl: Preparing to generate MSAs for '.$canonical_species;
} else {
    $progmsg = '  MultiMSA.pl: Preparing to generate MSAs';
}
while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
print "$progmsg\r" if (!$verbose);


# Naming convention for progress file
my $progressbase = 'MultiMSA.thread_progress.';


# We'll be breaking up the dataset based on gene names
my @AllGroupIDs  = sort keys %Hits;
my $TotNumGroups = @AllGroupIDs;
exit(0) if (!$TotNumGroups);
if ($TotNumGroups < $CPUs) { $CPUs = $TotNumGroups; }
my $portion;
if ($CPUs) { $portion = $TotNumGroups / $CPUs; }
else       { exit(0);                          } # decline gracefully


# Split into the requested number of (default 2) processes
my $processes = 1;
my $ThreadID  = 0; # In case we only have one process
my $pid;
while ($processes < $CPUs) {
    
    # Create new thread
    if ( $pid = fork ) {
	# Bail if we fail (nice rhyme!)
	if (not defined $pid) { die "\n\tFork failed\n\n"; }    
	# Giving the 'master' thread an ID of '0'
        $ThreadID = 0;
    } else {
        $ThreadID = $processes;
        last;
    }

    $processes++;

}


# Status printing stuff
my $progtime = time();
my $progfilename = $progressbase.$ThreadID;
my $num_complete = 0;


# Open up the file for writing misses (i.e., multi-chromosome
# hitting sequence groups) to.
# These files will all get appended into the same missfile
# generated by Quilter
my $missbase = 'MultiMSA.misses.';
my $missfilename = $missbase.$ThreadID;
open(my $missfile,'>',$missfilename) ||  die "\n  ERROR:  Failed to open file '$missfilename'\n\n";


# Now we can generate our MSAs (sorted for progress vis.)
my $startpoint = $ThreadID * $portion;
my $endpoint   = ($ThreadID+1) * $portion;
if ($ThreadID == $CPUs-1) { $endpoint = $TotNumGroups; }
foreach my $group_index ($startpoint..$endpoint-1) {
    
    my $group_id = $AllGroupIDs[$group_index];
    print "  $group_id\n" if ($verbose);
    
    my $outfilename = $outputfolder.'/'.$group_id;
    if ($specific) {
	$outfilename = $outfilename.'_'.$specific;
    }
    $outfilename = $outfilename.'.afa';

    my $numhits = $Nums{$group_id};

    # Figure out the dominant chromosome for this family
    my @ChrList = split(/\|/,$ChrsByFam{$group_id});
    my %TopChr;
    foreach my $chr (@ChrList) {
	if ($TopChr{$chr}) { $TopChr{$chr}++; }
	else               { $TopChr{$chr}=1; }
    }
    my $curr_chr = 0;
    my $curr_chr_counts = -1;
    foreach my $chr (keys %TopChr) {
	if ($TopChr{$chr} > $curr_chr_counts) {
	    $curr_chr = $chr;
	    $curr_chr_counts = $TopChr{$chr};
	}
    }

    my $revcomp;
    if ($curr_chr =~ /\[revcomp\]/) { $revcomp = 1; }
    else                            { $revcomp = 0; }

    
    my $ChrName;
    my @GeneNames;
    my @IsoNames;
    my @Species;
    my @DBEntries;
    my @PosStrings;
    my @ProteinSeqs;
    my @IsoIDs;
    my @ExtraInfo;
    my @GroupField; # Because there's the possibility of case varying
    
    my ($NuclStart,$NuclEnd);
    my $hash_entry = $Hits{$group_id};
    my @each_entry = split('&',$hash_entry);

    $i = 0;
    my $chr_hits = 0;
    foreach $i (0..$numhits-1) {

	my @this_entry  = split('#',$each_entry[$i]);
	$ChrName        = $this_entry[0];

	if ($ChrName ne $curr_chr) {

	    # Write this guy out for future alignment using translated
	    # Needleman-Wunsch
	    print $missfile "$this_entry[1]\n";	    
	    next;

	}

	$PosStrings[$chr_hits] = $this_entry[2];
	$DBEntries[$chr_hits]  = $this_entry[1];
	
	# Some entries might have additional data fields before the
	# group identifier, so we'll have to pay attention to how
	# many fields there are.
	if ($DBEntries[$chr_hits]  =~ /([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|(\S*\|)([^\|]+)\s*$/) {
	    $GeneNames[$chr_hits]  = $1;
	    $IsoNames[$chr_hits]   = $2;
	    $Species[$chr_hits]    = $3;
	    $IsoIDs[$chr_hits]     = $4;
	    $ExtraInfo[$chr_hits]  = $5;
	    $GroupField[$chr_hits] = $6;
	} elsif ($DBEntries[$chr_hits]  =~ /([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\s*$/) {
	    $GeneNames[$chr_hits]  = $1;
	    $IsoNames[$chr_hits]   = $2;
	    $Species[$chr_hits]    = $3;
	    $IsoIDs[$chr_hits]     = $4;
	    $ExtraInfo[$chr_hits]  = 0;
	    $GroupField[$chr_hits] = $5;
	} else {
	    die "\n  ERROR:  Failed to parse sequence name '$DBEntries[$chr_hits]'\n\n";
	}
	
	# Get a hold of the proteins and record them as strings
	my $eslsfetchCmd = "esl-sfetch $ARGV[1] \"$DBEntries[$chr_hits]\" \|";
	open(my $eslinput,$eslsfetchCmd) || die "\n  esl-sfetch failed to grab $DBEntries[$chr_hits] from $ARGV[1]\n\n";
	
	# Eat header line <- could be used for sanity check
	my $line = <$eslinput>;
	
	# Read in the protein sequence
	my $Protein = '';
	while ($line = <$eslinput>) {
	    $line =~ s/\n|\r//g;
	    next if (!$line);
	    $Protein = $Protein.$line;
	}
	close($eslinput);
	
	# Store the protein
	$Protein =~ s/\s//g;
	$ProteinSeqs[$chr_hits] = $Protein;
	$chr_hits++;	    
	
    }
    
    # Hash each codon to a designated nucleotide position
    my %MSA;
    my %MultiAminos; # Indicates where we've seen indels (multiple aminos on one codon)
    foreach $i (0..$chr_hits-1) {

	# Converting the string of positions into individual codon centers
	my @Codons = split(/\,/,$PosStrings[$i]);
	my @Seq    = split(//,$ProteinSeqs[$i]);

	# Walking through the position indices, filling in our hash of
	# chromosome positions to amino acids
	my $amino_pos = 0;
	my $last_codon = -1;
	foreach my $codon_pos (@Codons) {

	    next if ($codon_pos =~ /\*/);

	    my $amino = $Seq[$amino_pos];
	    $amino_pos++;

	    if ($MSA{$codon_pos}) {

		if ($last_codon == $codon_pos) { 

		    # Indel ahoy!
		    $MSA{$codon_pos} = $MSA{$codon_pos}.$amino;
		    $MultiAminos{$codon_pos} = 1;

		} else {
		    $MSA{$codon_pos} = $MSA{$codon_pos}.','.$i.':'.$amino; 
		}

	    } else {
		$MSA{$codon_pos} = $i.':'.$amino;
	    }

	    $last_codon = $codon_pos;

	}

    }

    
    # Generate a sorted list of our hash entries (positions in the table)
    my @PositionIndex = keys %MSA;
    if ($revcomp) { @PositionIndex = sort {$b <=> $a} @PositionIndex; }
    else          { @PositionIndex = sort {$a <=> $b} @PositionIndex; }



    ###############################################################################
    #                                                                             #
    #  From here, we're doing the actual MSA assembly work for this gene family.  #
    #  A couple things worth mentioning: (1.) The 'FinalMSA' isn't actually the   #
    #  final MSA at the end of this 'foreach' loop, since a little bit of work    #
    #  is required to make sure that any insertions get handled correctly.        #
    #                                                                             #
    ###############################################################################


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

    my $pos_index = 0;
    while ($pos_index < scalar(@PositionIndex)) {

	# Figuring out the range of indices in our PositionIndex
	# corresponding to the exon we're currently building.
	my ($next_index,$segment_ref) = GetNextExonRange($pos_index,\@PositionIndex);
	my @SegmentList = @{$segment_ref};

	foreach my $segment (@SegmentList) {

	    my @SegmentParts = split(/\,/,$segment);
	    my $start_index = $SegmentParts[0];
	    my $end_index   = $SegmentParts[1];
	    my $num_frames  = $SegmentParts[2];

	    # If there are multiple frames, we'll do a quick scan to find out which
	    # frame is the most common and designate that the non-ARF.
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
		    $msa_len += $longest_entry_len;

		}

		# If we just wrapped up an ARF frame, then we'll record how much 
		# each of the sequences in this frame moved.
		if ($frame != $non_arf_frame) {
		    for (my $seq_id=0; $seq_id<$chr_hits; $seq_id++) {
			if ($SeqsInFrame[$seq_id]) {
			    # To avoid weirdness associated with SPALN sticking single amino acids in places,
			    # we'll skip annotating single amino ARFs
			    if ($ContentStarts[$seq_id] < $ContentPositions[$seq_id]) {
				my $new_arf_entry = $ContentStarts[$seq_id].'..'.$ContentPositions[$seq_id];
				if ($ARFNameField[$seq_id]) {
				    $ARFNameField[$seq_id] =~ s/^ARF\:/ARFs\:/;
				    $ARFNameField[$seq_id] = $ARFNameField[$seq_id].','.$new_arf_entry;
				} else {
				    $ARFNameField[$seq_id] = 'ARF:'.$new_arf_entry;
				}
			    }
			}
		    }
		}

	    }

	}

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


    # Now that we have our clean and tidy MSA we can go
    # ahead and write it out!
    open(my $outfile,'>',$outfilename);
    for ($i=0; $i<$chr_hits; $i++) {

	$IsoNames[$i] =~ s/\s//g;
	if ($ExtraInfo[$i]) {
	    if ($ARFNameField[$i]) {
		print $outfile ">$GeneNames[$i]|$IsoNames[$i]|$Species[$i]|$IsoIDs[$i]|$ExtraInfo[$i]$ARFNameField[$i]|$GroupField[$i]\n";
	    } else {
		print $outfile ">$GeneNames[$i]|$IsoNames[$i]|$Species[$i]|$IsoIDs[$i]|$ExtraInfo[$i]$GroupField[$i]\n";
	    }
	} else {
	    if ($ARFNameField[$i]) {
		print $outfile ">$GeneNames[$i]|$IsoNames[$i]|$Species[$i]|$IsoIDs[$i]|$ARFNameField[$i]|$GroupField[$i]\n";
	    } else {
		print $outfile ">$GeneNames[$i]|$IsoNames[$i]|$Species[$i]|$IsoIDs[$i]|$GroupField[$i]\n";
	    }
	}
	    
	$j = 0;
	$k = 0;
	while ($j < $final_len) {
	    
	    print $outfile "$FinalMSA[$i][$j]";

	    $k++;
	    $j++;
	    
	    if ($k == 50) {
		print $outfile "\n";
		$k = 0;
	    }

	}

	print $outfile "\n";
	print $outfile "\n" if ($k);
	
    }
    close($outfile);
    

    ## If we had any disagreements let the user know
    #if (@Disagreements) {
    #
    #my $numdisagreements = @Disagreements;
    #my $disfilename = $outfilename;
    #
    #$disfilename =~ s/\.afa/\_ARFs/;
    ##print "  > WARNING: $numdisagreements positions did not reach unanimous consensus (see $disfilename)\n" if ($verbose);
    #
    #open(my $disfile,'>',$disfilename);
    #print $disfile "$numdisagreements mismatched sites\n";
    #print $disfile "$Disagreements[0]-";
    #foreach $j (1..$numdisagreements-1) {
    #if ($Disagreements[$j] != $Disagreements[$j-1]+1) {
    #print $disfile "$Disagreements[$j-1],$Disagreements[$j]-";
    #}
    #}
    #print $disfile "$Disagreements[$numdisagreements-1]\n";
    #close($disfile);
    #
    #}


    # Status reporting stuff:
    if (!$verbose) {
	$num_complete++;
	if (time()-$progtime > 2) {
	    
	    if ($ThreadID) {
		open(my $progfile,'>',$progfilename);
		print $progfile "$num_complete\n";
		close($progfile);
		
	    } else {
		
		# Generate call to ProgressTimer
		my $progressCmd = './src/ProgressTimer.pl '.$progressbase.' '.$CPUs;
		$progressCmd = $progressCmd.' '.$num_complete.' |';
		
		# Make the call and read the results
		open(my $progfile,$progressCmd);
		my $overall_progress = readline($progfile);
		close($progfile);
		$overall_progress =~ s/\n|\r//g;
		
		# Pull together the message
		my $progmsg = '  MultiMSA.pl: '.$overall_progress.' gene families aligned ('.$canonical_species.')';
		while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
		print "$progmsg\r";
		
	    }
	    
	    # Set new prog timer
	    $progtime = time();
	    
	}
    }
}


# Close up individual missfiles
close($missfile);


# de-fork
if ($ThreadID) { exit(0); }
while (wait() != -1) {}


# Clean up progress files and concatenate all missfiles into 'FinalMisses'
for ($i = 0; $i < $CPUs; $i++) {

    $progfilename = $progressbase.$i;
    if (-e $progfilename && system("rm $progfilename")) { die "\n  Failed to remove file $progfilename\n\n"; }

    $missfilename = $missbase.$i;
    if (-e $missfilename) {
	system("cat $missfilename >> $FinalMisses");
	system("rm $missfilename");
    }
    
}


# Final report
$progmsg = '  MultiMSA.pl: '.$canonical_species.' complete';
while (length($progmsg) < 63) { $progmsg = $progmsg.' '; }
print "$progmsg\r" if (!$verbose);


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
# FUNCTION: MAX
#
sub MAX
{
    my $a = shift;
    my $b = shift;
    return $a if ($a >= $b);
    return $b;
}





################################################################
#
# FUNCTION: MIN
#
sub MIN
{
    my $a = shift;
    my $b = shift;
    return $a if ($a <= $b);
    return $b;
}





################################################################
#
# FUNCTION: GetNextExonRange
#
sub GetNextExonRange
{
    my $start_index   = shift;
    my $position_ref  = shift;
    my @PositionIndex = @{$position_ref};

    # First, we just want to know how far we go from the start of this
    # exon before we hit an intron
    my $next_index = $start_index;
    my $last_nucl  = $PositionIndex[$start_index];
    $next_index++;
    while ($next_index < scalar(@PositionIndex) && abs($PositionIndex[$next_index]-$last_nucl)<10) {
	$last_nucl = $PositionIndex[$next_index];
	$next_index++;
    }

    # Next, we'll segment the exon based on any ARF-looking stuff.
    # Each segment will have (1.) a start position, (2.) a length,
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
	    push(@SegmentList,$segment_start.','.$index.',1');

	} else {

	    # Could this be... a... TRIPLE ARF?!
	    if ($index+1 < $next_index && abs($PositionIndex[$index+1]-$last_nucl)==2) {

		# IMPOSSIBLE!!!
		$index++;
		while ($index < $next_index && abs($last_nucl-$current_nucl)==1) {
		    $last_nucl = $current_nucl;
		    $current_nucl = $PositionIndex[$index];
		    $index++;
		}
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

		$index++;
		while ($index < $next_index && abs($last_nucl-$current_nucl) == $expected_jump) {
		    $last_nucl = $current_nucl;
		    $current_nucl = $PositionIndex[$index];
		    if ($expected_jump == 2) { $expected_jump = 1; }
		    else                     { $expected_jump = 2; }
		    $index++;
		}

		# In case the ARF doesn't cleanly line up with the end of the
		# exon, we won't count the first character we saw indicating that
		# the reading frame pattern changed.
		$index-- if (abs($last_nucl-$current_nucl) != $expected_jump); 

		push(@SegmentList,$segment_start.','.$index.',2');

	    }
	    
	}
	
    }

    return ($next_index,\@SegmentList);
    
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
	    if ($exon_conflicts > 4) {

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
    for (my $j=0; $j<$msa_len; $j++) {
	my ($num_chars,$top_char,$top_char_count) = CheckColumnProfile(\@MSA,$num_seqs,$j,-1,0);
	push(@Disagreements,$j) if ($num_chars > 1);
    }

    return(\@MSA,\@Disagreements,$msa_len);

}





################################################################
#
# FUNCTION: SortOutExon
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









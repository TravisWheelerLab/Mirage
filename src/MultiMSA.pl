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


sub MIN;         # Get the minimum of two values
sub MAX;         # Get the maximum of two values
sub AssembleMSA; # Pull together the final MSA
sub SimpleClean; # Very basic cleanup, geared towards correcting little SPALN mistakes
sub QuickAlign;  # Used to resolve indels


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
    
    my $ChrName;
    my $revcomp;
    my @GeneNames;
    my @IsoNames;
    my @Species;
    my @DBEntries;
    my @PosStrings;
    my @ProteinSeqs;
    my @IsoIDs;
    my @ExtraInfo;
    my @GroupField; # Because there's the possibility of case varying
    
    my ($curr_chr,$NuclStart,$NuclEnd);
    my $hash_entry = $Hits{$group_id};
    my @each_entry = split('&',$hash_entry);

    my $i = 0;
    my $chr_hits = 0;
    foreach $i (0..$numhits-1) {

	my $alt_chr = 0;
	
	my @this_entry  = split('#',$each_entry[$i]);
	$ChrName        = $this_entry[0];

	# Verifying that this group's isoforms all hit to the same chromosome
	if ($i == 0) { 

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
	    
	    $curr_chr = $ChrName;
	    if ($ChrName =~ /\[revcomp\]/) { $revcomp = 1; }
	    else                           { $revcomp = 0; }

	} elsif ($ChrName ne $curr_chr) {

	    # Write this guy out for future alignment using translated
	    # Needleman-Wunsch
	    print $missfile "$this_entry[1]\n";	    
	    $alt_chr = 1;

	} else {

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
	    
	}

	# No need to keep playing around with this one
	next if ($alt_chr);

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
	$ProteinSeqs[$chr_hits] = $Protein;
	$chr_hits++;	    
	
    }
    
    # It's worth remarking that even though we're hashing positions
    # we still need relative positioning, so it makes sense to use
    # ExonStarts and ExonEnds.

    # Now that we've got the overall min/max and each hit's string, we
    # can walk along, filling in an array for each of the strings.
    my %MSA;
    my %IntronTracker;
    foreach $i (0..$chr_hits-1) {

	# Converting the string of positions into individual codon centers
	my @Codons = split(',',$PosStrings[$i]);

	$ProteinSeqs[$i] =~ s/\s//g;
	my @Seq = split('',$ProteinSeqs[$i]);
		
	$j = 0; # How far along the series of codons we are
	$k = 0; # Position in Seq

	my $next_pos;
	my $last_write;
	my $last_read;

	# Walking through the position indices
	while ($j < @Codons) {

	    $next_pos = $Codons[$j];
	    
	    if ($next_pos =~ /\*/) {

		# We don't want to double-down on splice boundaries
		if ($Codons[$j-1] ne '*') {

		    # Which direction are we going?
		    my $next_intron = $Codons[$j-1];
		    if ($revcomp) { $next_intron -= 3; }
		    else          { $next_intron += 3; }

		    # Well, better record the dang thing!
		    $IntronTracker{$next_intron} = 1;
		    
		}

	    } else {
		
		# If we're walking through the shadow of the valley of indels
		# take special notice.
		if ($last_read && abs($last_read-$next_pos) < 3) {


		    # We'll stack the next character on top of this
		    # one and deal with it later.
		    $MSA{$last_write} = $MSA{$last_write}.$Seq[$k];


		} else { #===============[ normal entry ]================#

		    # Standard placement, adding onto existing position
		    if ($MSA{$next_pos}) {
			$MSA{$next_pos} = $MSA{$next_pos}.','.$i.':'.$Seq[$k];
			$last_write = $next_pos;

			
			# Special considerations for weird cases (like mouse OBSCN) where
			# we have regions that don't quite align.  Note that we do a
			# little checking to make sure we aren't doubling down
		    } elsif ($MSA{$next_pos+1} && $MSA{$next_pos+1} !~ /\,$i\:/) {
			$MSA{$next_pos+1} = $MSA{$next_pos+1}.','.$i.':'.$Seq[$k];
			$last_write = $next_pos+1;
			
		    } elsif ($MSA{$next_pos-1} && $MSA{$next_pos-1} !~ /\,$i\:/) {
			$MSA{$next_pos-1} = $MSA{$next_pos-1}.','.$i.':'.$Seq[$k];
			$last_write = $next_pos-1;

			
			# New placement -- note that we need to start with a comma for the
			# preceding 'elsif's to work in every case (we'll trim this later)
		    } else {
			$MSA{$next_pos} = ','.$i.':'.$Seq[$k];
			$last_write = $next_pos;
		    }
		    
		}
	    
		
		$k++;

		# NOTE: This catch shouldn't be necessary...
		#last if ($k >= @Seq);
		
	    }

	    # Record what we just worked with as what we last saw
	    $last_read = $next_pos if ($next_pos ne '*');

	    # Next stop, next codon (but hopefully not next stop codon!)
	    $j++;

	}

    }

    
    # Generate a sorted list of our hash entries (positions in the table)
    my @PositionIndex = keys %MSA;
    my @IntronIndex   = keys %IntronTracker;
    if ($revcomp) { 
	@PositionIndex = sort {$b <=> $a} @PositionIndex; 
	@IntronIndex   = sort {$b <=> $a} @IntronIndex;
    } else { 
	@PositionIndex = sort {$a <=> $b} @PositionIndex; 
	@IntronIndex   = sort {$a <=> $b} @IntronIndex;
    }



    ###############################################################################
    #                                                                             #
    #  From here, we're doing the actual MSA assembly work for this gene family.  #
    #  A couple things worth mentioning: (1.) The 'FinalMSA' isn't actually the   #
    #  final MSA at the end of this 'foreach' loop, since a little bit of work    #
    #  is required to make sure that any insertions get handled correctly.        #
    #                                                                             #
    ###############################################################################

    
    
    # Run along the MSA hash putting together a rough sketch of the alignment
    #
    my @FinalMSA = (); # Psssst... This isn't actually the final MSA, yet...
    my $msa_len  = 1;

    # Prime the MSA with a column of splice-site markers
    #
    for ($i=0; $i<$chr_hits; $i++) { $FinalMSA[$i][0] = '*'; }

    # We should be able to fill out a pretty darn nice MSA from all this!
    #
    my $aa_index     = 0;
    my $intron_index = 0;
    my $next_intron_pos = $IntronIndex[$intron_index];
    while ($aa_index < @PositionIndex) {

	# Where are we in the genome?
	my $genome_pos = $PositionIndex[$aa_index];
	$aa_index++;

	# Have we hit an intron boundary 
	if ($intron_index < @IntronIndex 
	    && ((!$revcomp && $next_intron_pos < $genome_pos)
		||($revcomp && $next_intron_pos > $genome_pos))) {

	    # If there's an imminent boundary, we'll hold off until that one
	    $intron_index++;
	    if ($intron_index < @IntronIndex) {
		
		if ((!$revcomp && $IntronIndex[$intron_index-1]+10 < $IntronIndex[$intron_index])
		    || ($revcomp && $IntronIndex[$intron_index-1]-10 > $IntronIndex[$intron_index])) {
		    
		    for ($i=0; $i<$chr_hits; $i++) { $FinalMSA[$i][$msa_len] = '*'; }
		    $msa_len++;
		    
		}
		
		$next_intron_pos = $IntronIndex[$intron_index];

	    } else {

		for ($i=0; $i<$chr_hits; $i++) { $FinalMSA[$i][$msa_len] = '*'; }
		$msa_len++;

	    }
	    
	}

	# Next position up!  We need to do some careful work if there's a stack
	# (i.e., indel indicator) associated with this position.
	my @Seqs;
	my @Chars;

	my $cleaned_entry = $MSA{$genome_pos}; # We need to chop off the leading comma first
	$cleaned_entry =~ s/^\,//;
	my @Entries = split(',',$cleaned_entry);

	my $num_entries = @Entries;
	my @MultiChars;
	my $num_multis = 0;
	for ($i=0; $i<$num_entries; $i++) {

	    # Split each of the items into the sequence ID number and
	    # its character at this position.
	    my @Entry  = split(':',$Entries[$i]);
	    $Seqs[$i]  = $Entry[0];
	    $Chars[$i] = $Entry[1];

	    # Check for multiple characters being stacked at this position
	    if (length($Chars[$i]) > 1) {
		push(@MultiChars,$i);
		$num_multis++;
	    }
	    
	}

	# Is this an easy one (i.e., everyone is only one character)?
	if ($num_multis == 0) {

	    # If we're in the seqs, we get a present, otherwise
	    # a gap
	    for ($i=0; $i<$chr_hits; $i++)    { $FinalMSA[$i][$msa_len]        = '-';        }
	    for ($i=0; $i<$num_entries; $i++) { $FinalMSA[$Seqs[$i]][$msa_len] = $Chars[$i]; }
	    $msa_len++;
	    
	} else { # Dangerous territory!  Mind your insertions!

	    # If this isn't the last position then we can (try to) get sneaky
	    if ($aa_index < @PositionIndex) {

		# Sample the acceptable characters at this site
		$cleaned_entry = $MSA{$PositionIndex[$aa_index]};
		$cleaned_entry =~ s/^\,//;
		my @NextEntries = split(',',$cleaned_entry);
		my %AcceptedChars;
		for ($i=0; $i<@NextEntries; $i++) {
		    my @NextEntry = split(':',$NextEntries[$i]);
		    if (length($NextEntry[1]) == 1) { 
			$AcceptedChars{lc($NextEntry[1])} = 1; 
		    }
		}
		
		# Grab the next genome position, if there is one, and see
		# if we can align there (cascadery) <- for each insert-y seq.
		my @UnresolvedMultis;
		my $num_unresolved=0;
		for ($i=0; $i<$num_multis; $i++) {

		    $j = $MultiChars[$i];

		    my @PosStr = split('',$Chars[$j]);
		    my $second = $PosStr[1];

		    # If this character can be used in the next position, do some
		    # swaggadelic cascadery!
		    if ($AcceptedChars{lc($second)}) {

			my $remainder = '';
			for ($k=1; $k<@PosStr; $k++) {$remainder=$remainder.$PosStr[$k];}

			$Chars[$j] = $PosStr[0];
			
			# We decide what to do based next based on whether this
			# sequence occurs in the next amino-acid associated site.
			my $next_entry = $MSA{$PositionIndex[$aa_index]};
			if ($next_entry =~ /\,$Seqs[$j]\:/) {
			    $next_entry =~ s/\,$Seqs[$j]\:/\,$Seqs[$j]\:$remainder/;
			} else {
			    $next_entry = $next_entry.','.$Seqs[$j].':'.$remainder;
			}
			$MSA{$PositionIndex[$aa_index]} = $next_entry;
			    
		    } else {

			# Very sad :'(
			push(@UnresolvedMultis,$j);
			$num_unresolved++;
			
		    }

		}

		# If we have an Multi-character things that didn't get resolved
		# we'll have to take care of them here using our 'QuickAlign'
		if ($num_unresolved) {

		    # We need to know who among the unresolved is longest
		    my $max_unres_id  = $Seqs[$UnresolvedMultis[0]];
		    my $max_unres_len = length($Chars[$UnresolvedMultis[0]]);
		    for ($i=1; $i<$num_unresolved; $i++) {
			my $comp_len = length($Chars[$UnresolvedMultis[$i]]);
			if ($comp_len > $max_unres_len) {
			    $max_unres_len = $comp_len;
			    $max_unres_id  = $Seqs[$UnresolvedMultis[$i]];
			}
		    }

		    # Let's load up a column with every sequence.  We need to
		    # start by loading every position with gaps in case our
		    # cascading has put things out of order
		    my @Column;
		    for ($i=0; $i<$chr_hits; $i++)    { $Column[$i]        = '-';        }
		    for ($i=0; $i<$num_entries; $i++) { $Column[$Seqs[$i]] = $Chars[$i]; }

		    # Pseudo-align those suckers!
		    my $columnref = QuickAlign(\@Column,$max_unres_id,$max_unres_len,$chr_hits);
		    @Column = @{$columnref};

		    # And load up that FinalMSA with all sorts of tasty treats!
		    for ($i=0; $i<$chr_hits; $i++) {
			my @Stuff = split('',$Column[$i]);
			for ($j=0; $j<$max_unres_len; $j++) {
			    $FinalMSA[$i][$msa_len+$j] = $Stuff[$j];
			}
		    }
		    $msa_len += $max_unres_len;

		} else { # Total resolution is a beautiful thing!

		    # It's just like not having a multi-character pile-up!
		    for ($i=0; $i<$chr_hits; $i++)    { $FinalMSA[$i][$msa_len]        = '-';        }
		    for ($i=0; $i<$num_entries; $i++) { $FinalMSA[$Seqs[$i]][$msa_len] = $Chars[$i]; }
		    $msa_len++;		    
		    
		}
		
	    } else { # The last genome position with associated amino acids

		# Do our standard fallback pseudo-alignment
		# We need to know who among the unresolved is longest
		my $max_multi_id  = $Seqs[$MultiChars[0]];
		my $max_multi_len = length($Chars[$MultiChars[0]]);
		for ($i=1; $i<$num_multis; $i++) {
		    my $comp_len = length($Chars[$MultiChars[$i]]);
		    if ($comp_len > $max_multi_len) {
			$max_multi_len = $comp_len;
			$max_multi_id  = $Seqs[$MultiChars[$i]];
		    }
		}
		
		# Let's load up a column with every sequence.  We need to
		# start by loading every position with gaps in case our
		# cascading has put things out of order
		my @Column;
		for ($i=0; $i<$chr_hits; $i++)    { $Column[$i]       = '-';        }
		for ($i=0; $i<$num_entries; $i++) { $Column[$Seqs[$i]]= $Chars[$i]; }
		
		# Pseudo-align those suckers!
		my $columnref = QuickAlign(\@Column,$max_multi_id,$max_multi_len,$chr_hits);
		@Column = @{$columnref};
		
		# And load up that FinalMSA with all sorts of tasty treats!
		for ($i=0; $i<$chr_hits; $i++) {
		    my @Stuff = split('',$Column[$i]);
		    for ($j=0; $j<$max_multi_len; $j++) {
			$FinalMSA[$i][$msa_len+$j] = $Stuff[$j];
		    }
		}
		$msa_len += $max_multi_len;
		
	    }

	}
	
    }

    # Wrap things up with a column of splice-site markers
    for ($i=0; $i<$chr_hits; $i++) { $FinalMSA[$i][$msa_len] = '*'; }
    $msa_len++;



    
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
	my ($FinalMSARef,$DisagreementsRef,$tmp_len) = AssembleMSA(\@FinalMSA,$chr_hits,$msa_len);
	@FinalMSA = @{$FinalMSARef};
	@Disagreements = @{$DisagreementsRef};
	$final_len = $tmp_len;
    }


    # Unless we want visible ARF-candidates, we'll run through the disagreements to
    # get ARFs into their own groups.
    #
    #    BEGIN ARF-SPLITTING HERE
    #   +========================+
    #
    my @ARFNameField;
    for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) { $ARFNameField[$seq_id] = 0; }
    if ($mask_arfs && scalar(@Disagreements)) {
	
	# Run through the disagreement position indices and check for 
	# any 5+ runs of contiguous disagreements (to account for spurious
	# matches in ARFs).
	#
	my @ARFIndices;
	my $dis_index = 0;
	while ($dis_index < scalar(@Disagreements)) {
	    
	    my $runner = $dis_index+1;
	    $runner++ while ($runner < scalar(@Disagreements) && abs($Disagreements[$runner] - $Disagreements[$runner-1]) <= 4);

	    if ($runner - $dis_index > 5) { push(@ARFIndices,$Disagreements[$dis_index].'-'.$Disagreements[$runner-1]); }
	    
	    $dis_index = $runner;

	}

	# If we have ARF indices, we'll need to make some adjustments...
	# Included in this are considerations for how much we'll be increasing
	# the length of the MSA as we go along.
	#
	if (scalar(@ARFIndices)) {

	    my @ContentPosition;
	    for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) { 
		$ContentPosition[$seq_id] = 0; 
	    }

	    my @ARFMSA;
	    my $arfmsalen = 0;
	    my $orig_pos  = 0;
	    my $last_end  = 0;
	    foreach my $msa_arf_pair (@ARFIndices) {
		
		#
		#  In this loop we scan through each of the ARF indices (in the "original final" MSA)
		#  and [1.] get caught up (in our ARFMSA) to the starting point, [2.] figure out if
		#  we need to resolve (i.e., break apart) sequences to make the ARF not an abomination
		#  for our alignment, and [3.] break it up if we decide that's what makes sense.
		#

		my @ARFPair = split(/\-/,$msa_arf_pair);
		my $start   = $ARFPair[0];
		my $end     = $ARFPair[1]+1; # NOTE that ARFPair[1] is the INTERNAL end point!
		my $arf_len = $end-$start;
		$last_end   = $end;

		while ($orig_pos < $start) {
		    for (my $seq_id=0; $seq_id<$chr_hits; $seq_id++) {
			$ARFMSA[$seq_id][$arfmsalen] = $FinalMSA[$seq_id][$orig_pos];
			$ContentPosition[$seq_id]++ if ($FinalMSA[$seq_id][$orig_pos] =~ /[A-Za-z]/);
		    }
		    $orig_pos++;
		    $arfmsalen++;
		}
		
		my @ContentStart;
		for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) { 
		    $ContentStart[$seq_id] = $ContentPosition[$seq_id]; 
		}

		# We'll cluster sequences by finding out (1.) which sequence in the
		# ARF-candidate region is most common, and (2.) dividing between
		# sequences that are like this (>50% ID) and those that are < 50%.
		#
		my %CandidateSeqHash;
		my @CandidateSeqs;
		my @Intronic;
		for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
		    my $letter_count = 0;
		    my $nonletter_count = 0;
		    my $seq_str = '';
		    for (my $pos = $start; $pos < $end; $pos++) {
			my $next_char = uc($FinalMSA[$seq_id][$pos]);
			$seq_str = $seq_str.$next_char;
			if ($next_char =~ /[A-Za-z]/) {
			    $ContentPosition[$seq_id]++;
			    $Intronic[$pos-$start] = 0;
			    $letter_count++;
			} else {
			    $nonletter_count++;
			    if ($next_char eq '*')  { $Intronic[$pos-$start] = 1; }
			}
		    }
		    
		    # Any sequences that don't have any sequence here will just
		    # go to the most common (i.e., non-ARF-candidate) sequence.
		    #
		    if ($nonletter_count < $letter_count) { 
			$CandidateSeqs[$seq_id] = $seq_str;
			if ($CandidateSeqHash{$seq_str}) { $CandidateSeqHash{$seq_str}++; }
			else                             { $CandidateSeqHash{$seq_str}=1; }
		    } else {
			$CandidateSeqs[$seq_id] = 0;
		    }
		}

		# Get a list of intron positions (if there were any)
		#
		my @IntronPositions;
		for (my $pos = 0; $pos<$arf_len; $pos++) { 
		    if ($Intronic[$pos]) { push(@IntronPositions,$pos); } 
		}
		push(@IntronPositions,$arf_len);
		
		# Figure out which of the candidates is the canonical non-ARF
		#
		my $top_seq   = '!';
		my $top_count = 0;
		foreach my $cand_seq (keys %CandidateSeqHash) {
		    my $cand_count = $CandidateSeqHash{$cand_seq};
		    if ($cand_count > $top_count) {
			$top_count = $cand_count;
			$top_seq = $cand_seq;
		    }
		}

		# Not sure why this would happen, but whatever...
		#
		next if ($top_seq eq '!');
		
		# Alright -- who's aligned with the canonical non-ARF?
		#
		my @IsNonARF;
		my @ReferenceChars = split(//,$top_seq);
		my $ref_non_gaps = 0;
		foreach my $ref_char (@ReferenceChars) {
		    $ref_non_gaps++ if ($ref_char =~ /[a-zA-Z]/);
		}

		my $num_arf_seqs = 0; # In case we just want to leave this alone...
		for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
		    if ($CandidateSeqs[$seq_id]) {

			my $id_positions = 0;
			my @CandidateChars = split(//,$CandidateSeqs[$seq_id]);
			for (my $pos = 0; $pos < $arf_len; $pos++) {
			    next if ($ReferenceChars[$pos] !~ /[a-zA-Z]/);
			    $id_positions++ if ($CandidateChars[$pos] eq $ReferenceChars[$pos]);
			}

			if ($id_positions / $ref_non_gaps > 0.75) {
			    $IsNonARF[$seq_id] = 1;
			} else {
			    $IsNonARF[$seq_id] = 0;
			    $num_arf_seqs++;
			}

		    } else {

			# Gappy or low-content sequences will be noted as
			# non-ARFs.
			$IsNonARF[$seq_id] = 1;
			
		    }
		}

		# If we saw at least one strong arf candidate then we adjust the 
		# MSA -- otherwise leave it alone!
		#
		# Note that by virtue of how we tabulate things this pseudo-arf
		# will automatically be added to the alignment if we skip this
		# conditional.
		#
		if ($num_arf_seqs) {

		    for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {

			my $non_arf  = $IsNonARF[$seq_id];
			my $full_str = '';
			my $char_str = '';
			my $gap_str  = '';
			my $pos      = 0;
			my $exon_num = 0;
			while ($pos < $arf_len) {
			    if ($pos == $IntronPositions[$exon_num]) {
				$exon_num++;
				if ($non_arf) { $full_str = $full_str.'*'.$char_str.$gap_str; }
				else          { $full_str = $full_str.'*'.$gap_str.$char_str; }
				$char_str = '';
				$gap_str  = '';
			    } else {
				$char_str = $char_str.$FinalMSA[$seq_id][$start+$pos];
				$gap_str  = $gap_str.'-';
			    }
			    $pos++;
			}

			# Record the altered sequence!  Also, record the range if we're the ARF
			#
			if ($non_arf) { 
			    if (scalar(@IntronPositions) > 1) { $full_str = $full_str.'*'.$char_str.$gap_str; }
			    else                              { $full_str = $char_str.$gap_str;               }
			} else { 
			    if (scalar(@IntronPositions) > 1) { $full_str = $full_str.'*'.$gap_str.$char_str; }
			    else                              { $full_str = $gap_str.$char_str;               }

			    # NOTE: We're going to use the terrible counting scheme where
			    #       we start counting at 1, instead of 0 (which is why we
			    #       use "ContentStart[$seq_id]+1")
			    #
			    if ($ARFNameField[$seq_id]) {
				$ARFNameField[$seq_id] =~ s/ARF\:/ARFs\:/;
				$ARFNameField[$seq_id] = $ARFNameField[$seq_id].','.($ContentStart[$seq_id]+1).'..'.$ContentPosition[$seq_id];
			    } else {
				$ARFNameField[$seq_id] = 'ARF:'.($ContentStart[$seq_id]+1).'..'.$ContentPosition[$seq_id];
			    }

			}

			# NOW we can add in the ARF-y bit to the MSA
			#
			$pos = $arfmsalen;	
			foreach my $char (split(//,$full_str)) { $ARFMSA[$seq_id][$pos++] = $char; }

			# Shifting our 'big-picture' of the MSA forward
			if ($seq_id == $chr_hits-1) { $arfmsalen = $pos; }

		    }
		    
		} else {

		    # Just a straightforward copy from 'FinalMSA'
		    #
		    while ($start < $end) {
			for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {			    
			    $ARFMSA[$seq_id][$arfmsalen] = $FinalMSA[$seq_id][$start];
			}
			$start++;
			$arfmsalen++;
		    }
		    
		}

		$orig_pos = $end;

	    }


	    # We'll need to really quickly incorporate everything after the last
	    # ARF candidate.
	    #
	    while ($last_end < $final_len) {
		for (my $seq_id = 0; $seq_id < $chr_hits; $seq_id++) {
		    $ARFMSA[$seq_id][$arfmsalen] = $FinalMSA[$seq_id][$last_end];
		}
		$arfmsalen++;
		$last_end++;
	    }

	    # Swap ARFMSA into the FinalMSA variables
	    #
	    # NOTE:  It's been observed that there are some pathological ARFs that
	    #        are recognized as having inserted AAs relative to the genome,
	    #        so we need to do ANOTHER check for all-gap columns
	    #
	    $final_len = 0;
	    for (my $j=0; $j<$arfmsalen; $j++) {
		my $not_all_gaps = 0;
		for (my $i=0; $i<$chr_hits; $i++) {
		    $FinalMSA[$i][$final_len] = $ARFMSA[$i][$j];
		    $not_all_gaps = 1 if ($ARFMSA[$i][$j] ne '-');
		}
		$final_len += $not_all_gaps;
	    }
	    
	}

    }
    #
    #     END OF ARF-SPLITTING
    #    +====================+
    #


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
    

    # If we had any disagreements let the user know
    if (@Disagreements) {

	my $numdisagreements = @Disagreements;
	my $disfilename = $outfilename;

	$disfilename =~ s/\.afa/\_ARFs/;
	#print "  > WARNING: $numdisagreements positions did not reach unanimous consensus (see $disfilename)\n" if ($verbose);

	open(my $disfile,'>',$disfilename);
	print $disfile "$numdisagreements mismatched sites\n";
	print $disfile "$Disagreements[0]-";
	foreach $j (1..$numdisagreements-1) {
	    if ($Disagreements[$j] != $Disagreements[$j-1]+1) {
		print $disfile "$Disagreements[$j-1],$Disagreements[$j]-";
	    }
	}
	print $disfile "$Disagreements[$numdisagreements-1]\n";
	close($disfile);

    }


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
# FUNCTION: AssembleMSA
#
sub AssembleMSA
{

    my ($i,$j,$k);

    
    # Grab what we have to work with, so far
    #
    my $InputMSAref = shift;
    my @InputMSA = @{$InputMSAref};
    

    # The dimensions of the input MSA
    #
    my $num_seqs  = shift;
    my $input_len = shift;


    # The final MSA -- same number of seqs (duh), but possibly
    # a different length, depending on how things clean up...
    #
    my @FinalMSA;
    my $final_len = 0;


    # Walk along the input MSA, taking care of business
    #
    for ($j=0; $j<$input_len; $j++) {

	
	my @Column = ();
	my $max_col_ID   = -1;
	my $max_col_size = 0;

	
	# Get a quick look at this input column
	#
	for ($i=0; $i<$num_seqs; $i++) {
	    
	    my $next = $InputMSA[$i][$j];
	    push(@Column,$next);
	    
	    if (length($next) > $max_col_size) {
		$max_col_size = length($next);
		$max_col_ID   = $i;
	    }
	    
	}

	
	# Any multiple-character entries to resolve?
	#
	if ($max_col_size > 1) {

	    # Do some quick NW-esque alignment
	    #
	    my $ColumnRef = QuickAlign(\@Column,$max_col_ID,$max_col_size,$num_seqs);
	    @Column = @{$ColumnRef};

	}

	
	# Store everything right proper!
	#
	for ($i=0; $i<$num_seqs; $i++) {

	    my @Row = split('',$Column[$i]);
	    for ($k=0; $k<$max_col_size; $k++) { $FinalMSA[$i][$final_len+$k] = $Row[$k];}
	    
	}

	
	# And, finally, increment the length of the final sequence
	#
	$final_len += $max_col_size;
	
	
    }

    
    # Pass this over to a cleanup function that checks for obvious spacing goofs.
    #
    # This function gives us back everything that we need, so we can just forward
    # the results.
    #
    return SimpleClean(\@FinalMSA,$num_seqs,$final_len);


}






################################################################
#
# FUNCTION: SimpleClean
#
sub SimpleClean
{

    # Grab the MSA and its associated metadata
    #
    my $msa_ref = shift;
    my @MSA = @{$msa_ref};

    my $num_seqs = shift;
    my $msa_len  = shift;

    
    # The first thing we'll do is build a 'consensus sequence'
    # of all sites with unanimous agreement about their contents
    # (that is, all columns that store the same character) as
    # well as all columns that disagree somewhere (not counting
    # gaps).
    #
    my @ConsensusSeq;
    my @Disagreements;
    for (my $i=0; $i<$msa_len; $i++) {

	# Not interested in splice sites
	#
	if ($MSA[0][$i] eq '*') {
	    push(@ConsensusSeq,'*');
	    next;
	}

	# Run through each sequence, checking if we're in agreement
	#
	my $disagreement  = 0;
	my $consensuschar = 0;
	for (my $j=0; $j<$num_seqs; $j++) {

	    if ($MSA[$j][$i] =~ /\w/) {

		if (!$consensuschar) { 
		    $consensuschar = uc($MSA[$j][$i]);
		} elsif (uc($MSA[$j][$i]) ne $consensuschar) {
		    $disagreement = 1;
		    last;
		}

	    }

	}

	# How are we looking? Did we find a disagreement?
	#
	if ($disagreement) {
	    push(@ConsensusSeq,'#');
	    push(@Disagreements,$i);
	} else {
	    push(@ConsensusSeq,$consensuschar);    
	}

    }


    # If there weren't any disagreements we can just skip over this
    # whole error-correction stuff.
    #
    if (@Disagreements) {

	# Now that we have a list disagreeing positions, we can go through
	# and check whether there are any characters that we can swap
	# around to make those pesky disagreements disappear.
	#
	my @Unresolved;
	for (my $dis_num=0; $dis_num<scalar(@Disagreements); $dis_num++) {
	            
	    my $i = $Disagreements[$dis_num];
	            
	    # If there's a disagreement about the next column over (left) we
	    # don't consider placement.  Or if it's a splice site.
	    #
	    if ($ConsensusSeq[$i-1] eq '#') {
		push(@Unresolved,$i);
		next;
	    }
	    
	    # Iterate through each position checking if we can shift it
	    # safely into a consensus column
	    #
	    for (my $j=0; $j<$num_seqs; $j++) {
		
		# Check if the position to the left meets the following
		# criteria:
		#
		#     OPTION 1:
		#     1. Gap character
		#     2. Part of consensus column with this character
		#
		#     OPTION 2:
		#     1. Splice site marker immediately after a gap
		#     2. The gap before the splice site meets (1).
		#
		if ($MSA[$j][$i-1] =~ /\W/) {
		    if ($MSA[$j][$i-1] eq '-' && uc($MSA[$j][$i]) eq $ConsensusSeq[$i-1]) {
			$MSA[$j][$i-1] = $MSA[$j][$i];
			$MSA[$j][$i]   = '-';
		    } elsif ($MSA[$j][$i-1] eq '*' && $i-2 >= 0 
			     && $MSA[$j][$i-2] eq '-' &&  uc($MSA[$j][$i]) eq $ConsensusSeq[$i-2]) {
			$MSA[$j][$i-2] = $MSA[$j][$i];
			$MSA[$j][$i]   = '-';
		    }
		}
		
	    }
	            
	    # Now that we've corrected this column (as best we can) see if we've
	    # resolved what it should contain...
	    #
	    my $consensuschar = 0;
	    my $disagreement  = 0;
	    for (my $j=0; $j<$num_seqs; $j++) {
		
		next if ($MSA[$j][$i] =~ /\W/);
		
		if (!$consensuschar) {
		    $consensuschar = uc($MSA[$j][$i]);
		} elsif (uc($MSA[$j][$i]) ne $consensuschar) {
		    $disagreement = 1;
		    last;
		}
		
	    }
	        
	    # If there isn't a consensus character, that means that this column
	    # has become an all-gap column (and will need to be cleared out).
	    #
	    # If there's still disagreement we'll try to catch it on the flipside
	    # (i.e., when we run through in reverse checking right neighbors).
	    #
	    if (!$consensuschar) {
		$ConsensusSeq[$i] = '-';
		
	    } elsif ($disagreement) {
		push(@Unresolved,$i);
		
	    } else {
		$ConsensusSeq[$i] = $consensuschar;
		
	    }

	}
	
	
	# We're really cruising now!  Let's check the reverse...
	#
	for (my $dis_num=scalar(@Unresolved)-1; $dis_num>=0; $dis_num--) {
	            
	    my $i = $Unresolved[$dis_num];
	    
	    # If there's a disagreement about the next column over (right) we
	    # don't consider placement.  Again, we also don't mess with splice sites.
	    #
	    next if ($ConsensusSeq[$i+1] eq '#');
	    
	    # Iterate through each position checking if we can shift it
	    # safely into a consensus column
	    #
	    for (my $j=0; $j<$num_seqs; $j++) {
		
		# Check if the position to the left meets the criteria from before.
		#
		if ($MSA[$j][$i+1] =~ /\W/) {
		    if ($MSA[$j][$i+1] eq '-' && uc($MSA[$j][$i]) eq $ConsensusSeq[$i+1]) {
			$MSA[$j][$i+1] = $MSA[$j][$i];
			$MSA[$j][$i]   = '-';
		    } elsif ($MSA[$j][$i+1] eq '*' && $i+2 < $msa_len 
			     && $MSA[$j][$i+2] eq '-' &&  uc($MSA[$j][$i]) eq $ConsensusSeq[$i+2]) {
			$MSA[$j][$i+2] = $MSA[$j][$i];
			$MSA[$j][$i]   = '-';
		    }
		}
		
	    }
	            
	            
	    # We don't need to worry about locating unresolved disagreements
	    # (two passes is sufficient), but we benefit from knowing whether
	    # we've converted this into an all-gap-characters column.
	    #
	    my $consensuschar = 0;
	    my $disagreement  = 0;
	    for (my $j=0; $j<$num_seqs; $j++) {
		
		next if ($MSA[$j][$i] eq '-');
		
		if (!$consensuschar) {
		    $consensuschar = uc($MSA[$j][$i]);
		} elsif (uc($MSA[$j][$i]) ne $consensuschar) {
		    $disagreement = 1;
		    last;
		}
		
	    }
	            
	    # If there isn't a consensus character, that means that this column
	    # has become an all-gap column (and will need to be cleared out).
	    #
	    # If there's still disagreement we'll try to catch it on the flipside
	    # (i.e., when we run through in reverse checking right neighbors).
	    #
	    if (!$consensuschar) { 
		$ConsensusSeq[$i] = '-'; 
	    } elsif (!$disagreement) {
		$ConsensusSeq[$i] = $consensuschar;
	    }
	    
	}
    }
    
    
    # Actually, for one last sneak, let's see if we can't clean up any
    # jigsaw-type weirdness...
    #
    for (my $i=2; $i<$msa_len; $i++) {
	
	# Per usual, don't worry about those silly old splice sites
	#
	next if ($MSA[0][$i] eq '*');
	
	# Check if we can give it a little scoot
	#
	my $switcheroo = 0;
	for (my $j=0; $j<$num_seqs; $j++) {
	    if ($ConsensusSeq[$i-1] eq uc($MSA[$j][$i]) && $MSA[$j][$i-1] eq '-' && $MSA[$j][$i-2] =~ /\w/) {
		$MSA[$j][$i-1] = $MSA[$j][$i];
		$MSA[$j][$i]   = '-';
		$switcheroo = 1;
	    }
	}
	
	# If we switched this column around, maybe we just created an
	# all-gap column, which will need to be fixed...
	#
	if ($switcheroo) {
	    for (my $j=0; $j<$num_seqs; $j++) {
		if ($MSA[$j][$i] ne '-') {
		    $switcheroo = 0;
		    last;
		}
	    }
	    if ($switcheroo) { 
		$ConsensusSeq[$i] = '-';
	    }
	}
    }


    # EXCELLENT!
    #
    # Let's go ahead and make a final MSA without any all-gap columns
    # (or double-'*' columns).
    #
    my @FinalMSA;
    my @FinalDisagreements;
    my $final_len = 0;
    my $after_gap = 0;
    for (my $i=0; $i<$msa_len; $i++) {

	# Is this a splice-site boundary?
	#
	if ($ConsensusSeq[$i] eq '*') {
	        
	    # If we're immediately following another splice boundary,
	    # just forget about this column...
	    #
	    if (!$after_gap) {
		
		# We can just tear through this column.
		#
		for (my $j=0; $j<$num_seqs; $j++) { $FinalMSA[$j][$final_len] = '*'; }

		# One small step for this multiple sequence alignment.
		# One giant leap for...  Let me get back to you.
		#
		$final_len++;

	    }

	    # Splice-column observed!
	    #
	    $after_gap = 1;

	} elsif ($ConsensusSeq[$i] ne '-') {

	    # Looks like we're just looking at a regular old column.
	    # Boring!
	    #
	    $after_gap = 0;
	        
	    # Copy-paste (plagarism alert!)
	    #
	    for (my $j=0; $j<$num_seqs; $j++) { $FinalMSA[$j][$final_len] = $MSA[$j][$i]; }
	        
	    # Were there any disagreements at this column?
	    #
	    push(@FinalDisagreements,$i) if ($ConsensusSeq[$i] eq '#');
	        
	    # Moving forward with our lives
	    #
	    $final_len++;

	}

    }


    # Hope this helps!  Take care, now!
    #
    return (\@FinalMSA,\@FinalDisagreements,$final_len);

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














#######################  END OF FILE  ##########################









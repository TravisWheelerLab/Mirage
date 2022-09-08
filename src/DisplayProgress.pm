#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

# So every function doesn't need to loop up for length...
sub PrintProgress;
sub ClearProgress;

# Functions for setting global variables
sub InitMirageProgressVars;
sub InitQuilterProgressVars;
sub InitMapsToMSAsProgressVars;

# Many functions for specific outs-put
sub DispProgMirage;
sub DispProgQuilter;
sub DispProgMapsToMSAs;

# Before we get going, let's set a standardized character limit for status messages.
# Note that we go with a fairly large number to try to get at least one line of
# text-wrap. NEVERMIND, that causes the carriage-return to effectively work as a
# newline.
my $DispProg_line_len = 90;

# Variables that we'll want to have access to
my $DispProg_species;
my $DispProg_dirname;
my $DispProg_cpus;

my @QuilterFM2Genes;
my @QuilterB2SGenes;
my @QuilterCleanedGenes;
my @MapsToMSAsGenes;
my @MultiSeqNWGenes;



###############################################################
#
#  Functions: PrintProgress & ClearProgress
#
sub PrintProgress
{
    my $str = shift;
    $str = $str.' ' while (length($str) < $DispProg_line_len);
    print "$str\r";
}
sub ClearProgress { PrintProgress(' '); }




###############################################################
#
#  Function: InitMirageProgressVars
#
sub InitMirageProgressVars
{
    my $progress_dirname = shift;
    my $num_cpus = shift;

    # Record the location of the progress alcove
    $DispProg_dirname = $progress_dirname;

    # Copy over the number of cpus (requested -- still want to be careful!)
    $DispProg_cpus = $num_cpus;

    # Let 'em know we're initialized!
    DispProgMirage('init');    
}




###############################################################
#
#  Function: InitQuilterProgressVars
#
sub InitQuilterProgressVars
{
    my $seqdirname = shift;
    my $num_cpus = shift;

    # Reverse engineer the location of the progress alcove and species name
    # from the sequence directory name
    $seqdirname =~ /\/([^\/]+)\/seqs\/$/;
    $DispProg_species = $1;

    $seqdirname =~ s/[^\/]+\/$//; # seqs -> species
    $seqdirname =~ s/[^\/]+\/$//; # species -> SpeciesMSAs
    $seqdirname =~ s/[^\/]+\/$//; # SpeciesMSAs -> MirageResults    
    $seqdirname = $seqdirname.'.progress/';
    $DispProg_dirname = $seqdirname;

    # Copy over the number of cpus (requested -- still want to be careful!)
    $DispProg_cpus = $num_cpus;

    # Zero-out the gene completion counters
    for (my $i=0; $i<$num_cpus; $i++) {
	$QuilterFM2Genes[$i] = 0;
	$QuilterB2SGenes[$i] = 0;
	$QuilterCleanedGenes[$i] = 0;
    }

    # Let 'em know we're initialized!
    DispProgQuilter('init');
    
}



###############################################################
#
#  Function: InitMapsToMSAsProgressVars
#
sub InitMapsToMSAsProgressVars
{
    my $seqdirname = shift;
    my $num_cpus = shift;

    # Reverse engineer the location of the progress alcove and species name
    # from the sequence directory name
    $seqdirname =~ /\/([^\/]+)\/seqs\/$/;
    $DispProg_species = $1;

    $seqdirname =~ s/[^\/]+\/$//; # seqs -> species
    $seqdirname =~ s/[^\/]+\/$//; # species -> SpeciesMSAs
    $seqdirname =~ s/[^\/]+\/$//; # SpeciesMSAs -> MirageResults    
    $seqdirname = $seqdirname.'.progress/';
    $DispProg_dirname = $seqdirname;

    # Copy over the number of cpus (requested -- still want to be careful!)
    $DispProg_cpus = $num_cpus;

    # Zero-out the gene completion counters
    for (my $i=0; $i<$num_cpus; $i++) {
	$MapsToMSAsGenes[$i] = 0;
    }

    # Let 'em know we're initialized!
    DispProgMapsToMSAs('init');
    
}




###############################################################
#
#  Function: DispProgMirage
#
sub DispProgMirage
{
    my $data_str = shift;
    my @Data = split(/\|/,$data_str);

    # Field 0 is always the part of the program we're working on
    my $part = $Data[0];

    # During the MultiSeqNW-y stages of the pipeline, we might be threaded
    my $threadID = 0;
    $threadID = $Data[1] if (scalar(@Data) > 1);
    
    my $status = "  Mirage: ";
    if ($part eq 'init') {

	$status = $status."Checking files and performing setup";
	
    } elsif ($part eq 'db-speciation') {
	
	$status = $status."Dividing protein database according to species";
	
    } elsif ($part eq 'msnw-init') {

	for (my $i=0; $i<$DispProg_cpus; $i++) {
	    $MultiSeqNWGenes[$i] = 0;
	}

	$DispProg_species = $Data[2];
	if ($DispProg_species eq 'FINAL') {
	    $status = $status."Preparing to generate interspecies alignments";
	} else {
	    $status = $status."Preparing to join unmapped sequences to $DispProg_species MSAs";
	}

    } elsif ($part eq 'misc-ali') {

	$status = $status."Aligning sequences belonging to species without genomes";
	
    } elsif ($part eq 'msnw-loop') {

	# Make a wisdom saving throw (to avoid excessive printing)
	if (rand() < 0.3 && !($threadID==0 && $MultiSeqNWGenes[0]==0)) {
	    return;
	}

	# Write out how many genes you've completed to a secret file!
	my $genes_completed = $Data[2];
	my $outfbase = $DispProg_dirname.$DispProg_species.'.MSNW.';
	open(my $outf,'>',$outfbase.$threadID);
	print $outf "$genes_completed\n";
	close($outf);

	# If you're the master, count the total number of completed genes and
	# compose a status report
	if (!$threadID) {
	    for (my $i=1; $i<$DispProg_cpus; $i++) {
		my $outfname = $outfbase.$i;
		if (-e $outfname) {
		    open(my $inf,'<',$outfname);
		    my $thread_genes = <$inf>;
		    close($inf);
		    if ($thread_genes && $thread_genes =~ /(\d+)/) {
			$MultiSeqNWGenes[$i] = $1;
		    }
		    $genes_completed += $MultiSeqNWGenes[$i];
		}
	    }
	    $status = $status."$genes_completed MSAs constructed";
	    $status =~ s/ 1 MSAs / 1 MSA /;
	}

    }
    
    PrintProgress($status) if (!$threadID);
    
}






###############################################################
#
#  Function: DispProgQuilter
#
sub DispProgQuilter
{
    my $data_str = shift;
    my @Data = split(/\|/,$data_str);

    # The first part of our data is always going to be the part of the
    # program we're working on right now
    my $part = $Data[0];

    # The second piece will be the thread ID (if there is a second piece!)
    my $threadID = 0;    
    $threadID    = $Data[1] if (scalar(@Data) > 1);

    # In any case, we'll be using this classic line
    my $status = "  Quilter ($DispProg_species): ";

    # Now, what part are we working through currently?
    if ($part eq 'init') {

	$status = $status."Preparing to map protein sequences to genome";
	
    } elsif ($part eq 'parsing-gtf') {

	$status = $status."Loading GTF data from file";

    } elsif ($part eq 'fm2') {

	# Make a dexterity saving throw (to avoid excessive printing)
	if (rand() < 0.3 && !($threadID==0 && $QuilterFM2Genes[0]==0)) {
	    return;
	}

	# Write out how many genes you've completed to a secret file!
	my $genes_completed = $Data[2];
	my $outfbase = $DispProg_dirname.$DispProg_species.'.Quilter.fm2.';
	open(my $outf,'>',$outfbase.$threadID);
	print $outf "$genes_completed\n";
	close($outf);

	# If you're the master, count the total number of completed genes and
	# compose a status report
	if (!$threadID) {
	    for (my $i=1; $i<$DispProg_cpus; $i++) {
		my $outfname = $outfbase.$i;
		if (-e $outfname) {
		    open(my $inf,'<',$outfname);
		    my $thread_genes = <$inf>;
		    close($inf);
		    if ($thread_genes && $thread_genes =~ /(\d+)/) {
			$QuilterFM2Genes[$i] = $1;
		    }
		    $genes_completed += $QuilterFM2Genes[$i];
		}
	    }
	    $status = $status."$genes_completed genes examined using GTF coordinates";
	    $status =~ s/ 1 genes / 1 gene /;
	}

    } elsif ($part eq 'blatprep') {

	$status = $status."Preparing files for BLAT";
	
    } elsif ($part eq 'blatrunning') {

	$status = $status."Running BLAT";
	
    } elsif ($part eq 'blat2spaln') {

	# Make an insight check (to avoid excessive printing)
	if (rand() < 0.3 && !($threadID==0 && $QuilterB2SGenes[0]==0)) {
	    return;
	}

	# Write out how many genes you've completed to a different secret file!
	my $genes_completed = $Data[2];
	my $outfbase = $DispProg_dirname.$DispProg_species.'.Quilter.b2s.';
	open(my $outf,'>',$outfbase.$threadID);
	print $outf "$genes_completed\n";
	close($outf);

	# If you're the master, count the total number of completed genes and
	# compose a status report
	if (!$threadID) {
	    for (my $i=1; $i<$DispProg_cpus; $i++) {
		my $outfname = $outfbase.$i;
		if (-e $outfname) {
		    open(my $inf,'<',$outfname);
		    my $thread_genes = <$inf>;
		    close($inf);
		    if ($thread_genes && $thread_genes =~ /(\d+)/) {
			$QuilterB2SGenes[$i] = $1;
		    }
		    $genes_completed += $QuilterB2SGenes[$i];
		}
	    }
	    $status = $status."$genes_completed genes examined using BLAT coordinates";
	    $status =~ s/ 1 genes / 1 gene /;
	}

    } elsif ($part eq 'cleanup') {

	# Make a stealth check (to avoid excessive printing)
	if (rand() < 0.3 && !($threadID==0 && $QuilterCleanedGenes[0]==0)) {
	    return;
	}

	# Write out how many genes you've completed to a different secret file!
	my $genes_completed = $Data[2];
	my $outfbase = $DispProg_dirname.$DispProg_species.'.Quilter.cleanup.';
	open(my $outf,'>',$outfbase.$threadID);
	print $outf "$genes_completed\n";
	close($outf);

	# If you're the master, count the total number of completed genes and
	# compose a status report
	if (!$threadID) {
	    for (my $i=1; $i<$DispProg_cpus; $i++) {
		my $outfname = $outfbase.$i;
		if (-e $outfname) {
		    open(my $inf,'<',$outfname);
		    my $thread_genes = <$inf>;
		    close($inf);
		    if ($thread_genes && $thread_genes =~ /(\d+)/) {
			$QuilterCleanedGenes[$i] = $1;
		    }
		    $genes_completed += $QuilterCleanedGenes[$i];
		}
	    }
	    $status = $status."file cleanup complete for $genes_completed genes";
	    $status =~ s/ 1 genes/ 1 gene/;
	}

    }

    # Only the master gets to report progress (lucky duck!)
    PrintProgress($status) if (!$threadID);
    
}






###############################################################
#
#  Function: DispProgMapsToMSAs
#
sub DispProgMapsToMSAs
{
    my $data_str = shift;
    my @Data = split(/\|/,$data_str);

    my $part = $Data[0];

    # Do we have a thread ID?
    my $threadID = 0;
    $threadID = $Data[1] if (scalar(@Data) > 1);

    # Gotta love this one
    my $status = "  MapsToMSAs ($DispProg_species): ";

    # Do a lil' bit o' work
    if ($part eq 'init') {

	$status = $status."Preparing to generate alignments using genome mappings";
	
    } elsif ($part eq 'aligning') {

	# Make an athletics check (to avoid excessive printing)
	if (rand() < 0.3 && !($threadID==0 && $MapsToMSAsGenes[0]==0)) {
	    return;
	}

	# Write out how many genes you've completed to a secret file!
	my $genes_completed = $Data[2];
	my $outfbase = $DispProg_dirname.$DispProg_species.'.MapsToMSAs.';
	open(my $outf,'>',$outfbase.$threadID);
	print $outf "$genes_completed\n";
	close($outf);

	# If you're the master, count the total number of completed genes and
	# compose a status report
	if (!$threadID) {
	    for (my $i=1; $i<$DispProg_cpus; $i++) {
		my $outfname = $outfbase.$i;
		if (-e $outfname) {
		    open(my $inf,'<',$outfname);
		    my $thread_genes = <$inf>;
		    close($inf);
		    if ($thread_genes && $thread_genes =~ /(\d+)/) {
			$MapsToMSAsGenes[$i] = $1;
		    }
		    $genes_completed += $MapsToMSAsGenes[$i];
		}
	    }
	    $status = $status."$genes_completed genes aligned using genome mappings";
	    $status =~ s/ 1 genes / 1 gene /;
	}
	
    }

    PrintProgress($status) if (!$threadID);
    
}






1; # EOF

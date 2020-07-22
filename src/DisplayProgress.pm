#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

# So every function doesn't need to loop up for length...
sub PrintProgress;
sub ClearProgress;

# Functions for setting global variables
sub InitQuilterProgressVars;

# Many functions for specific outs-put
sub ProgressMirageInit;
sub ProgressMirageQuilter;
sub ProgressMirageMultiMSA;
sub ProgressMirageMultiSeqNW;
sub ProgressMirageCleanup;

# Before we get going, let's set a standardized character limit for status messages
my $DispProg_line_len = 120;

# Variables that we'll want to have access to
my $DispProg_species;
my $DispProg_dirname;
my $DispProg_cpus;





###############################################################
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
    
}



###############################################################
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
sub ProgressMirageInit
{
    # What part of the initialization are we immersed in?
    my $part = shift;
    my $status;

    if ($part eq 'db-speciation') {
	$status = "Dividing protein database according to species";
    }
    
    PrintProgress($status);
}






###############################################################
#
sub ProgressMirageQuilter
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
    if ($part eq 'parsing-gtf') {

	$status = $status."Loading GTF data from file";

    } elsif ($part eq 'fm2') {

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
		    if ($thread_genes =~ /(\d+)/) {
			$genes_completed += $1;
		    }
		}
	    }
	    $status = $status."$genes_completed genes examined using GTF coordinates";
	}

    } elsif ($part eq 'blatprep') {

	$status = $status."Preparing files for BLAT";
	
    } elsif ($part eq 'blatrunning') {

	$status = $status."Running BLAT";
	
    } elsif ($part eq 'blat2spaln') {

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
		    if ($thread_genes =~ /(\d+)/) {
			$genes_completed += $1;
		    }
		}
	    }
	    $status = $status."$genes_completed genes examined using BLAT coordinates";
	}

    }

    # Only the master gets to report progress (lucky duck!)
    PrintProgress($status) if (!$threadID);
    
}






###############################################################
#
sub ProgressMirageMultiMSA
{

}






###############################################################
#
sub ProgressMirageMultiSeqNW
{

}






###############################################################
#
sub MirageCleanup
{

}





1; # EOF

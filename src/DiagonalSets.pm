#
# DiagonalSet.pm - Alex Nord - 2016
#
# USED IN: Quilter.pl
#
# ABOUT: This perl module contains the DiagonalSet object prototype, along
#        with a function used to sort the contents of the object based on
#        its 'ProtStarts' array.
#
#
use warnings;
use strict;
use POSIX;


###########################################################################
#
# DiagonalSet - This object is used to track all segments of highly similar
#               sequence between a protein and a chromosome.  Because we're
#               doing translated search to build this, we need to track what
#               offset (0, 1, or 2) we're using and whether we encountered a
#               stop codon.
#
package DiagonalSet;

sub New
{
    my $class = shift;
    my $self  = {};
    
    $self->{ChrName}            = shift;
    $self->{NumHits}            = 0;
    $self->{NumTerminals}       = 0;
    $self->{ProtStrings}        = ();
    $self->{ProtStarts}         = ();
    $self->{ProtEnds}           = ();
    $self->{StartOffsetLengths} = ();
    $self->{StartOffsets}       = ();
    $self->{NuclStrings}        = ();
    $self->{NuclStarts}         = ();
    $self->{NuclEnds}           = ();
    $self->{EndOffsets}         = ();
    $self->{EndOffsetLengths}   = ();
    $self->{Scores}             = ();
    $self->{StopCodon}          = ();
    $self->{IsTerminal}         = ();
    
    bless $self,$class;
    return $self;
}




###########################################################################
#
# FUNCTION NAME:  SortHitsByStart
#
# ABOUT: SortHitsByStart simply sorts all the array fields of a DiagonalSet
#        based on their starting positions in the protein (ProtStarts).
#        The sorting is determined using a simple implementation of 
#        mergesort that works indirectly using an index, which may then
#        be applied to each of the fields in the DiagonalSet.
#
sub SortHitsByField
{
    my ($i,$j,$k);

    # The set we're sorting
    my $self    = shift;
    my $numHits = $self->{NumHits};

    # What field are we sorting on?
    my $field    = shift;
    my @fieldRef;
    if ($field eq 'ProtStarts') {
	foreach $i (0..$numHits-1) {
	    $fieldRef[$i] = ${$self->{ProtStarts}}[$i];
	}
    } elsif ($field eq 'Length') {
	foreach $i (0..$numHits-1) {
	    $fieldRef[$i] = ${$self->{ProtEnds}}[$i] - ${$self->{ProtStarts}}[$i];
	}
    } else {
	die "\n\tUnrecognized field option '$field'\n\n";
    }

    # Do we want the field to be sorted in ascending order?
    my $ascending = shift;
    
    # Create an index array.  This is used so we don't jostle things around
    # during the sorting step.
    my @index;
    foreach $i (0..$numHits-1) { $index[$i] = $i; }
    
    # Run mergesort on the set, using ProtStarts to determine order.
    my ($start,$end,$placer);
    my @temp;
    my $groupSize = 1;
    while ($groupSize < $numHits) {
	foreach my $groupID (0..POSIX::ceil(($numHits/(2*$groupSize)))-1) {
	    $start = $groupID * (2 * $groupSize);
	    $end   = $start + (2 * $groupSize) - 1;
	    $end   = $numHits-1 if ($end >= $numHits);
	    foreach $i (0..$end-$start) {
		$temp[$i] = $index[$start+$i];
	    }
	    $i = 0;
	    $j = $groupSize;
	    $placer = $start;
	    while ($i < $groupSize && $j <= $end-$start) {
		if ($fieldRef[$temp[$i]] < $fieldRef[$temp[$j]]) {
		    $index[$placer] = $temp[$i];
		    $i++;
		} else {
		    $index[$placer] = $temp[$j];
		    $j++;
		}
		$placer++;
	    }
	    if ($i < $groupSize) {
		while ($i < $groupSize) {
		    $index[$placer] = $temp[$i];
		    $i++;
		    $placer++;
		}
	    } else {
		while ($j <= $end-$start) {
		    $index[$placer] = $temp[$j];
		    $j++;
		    $placer++;
		}
	    }
	}
	$groupSize *= 2;
    }
    
    $self->OrganizeOnIndex(\@index,$ascending);

}


sub OrganizeOnIndex
{
    my ($i,$j,$k);

    my $self = shift;
    my $numHits = $self->{NumHits};

    my $indexRef = shift;
    my @index    = @{$indexRef};
    
    my $ascending = shift;

    my @temp;

    # Now that the sorting is complete, we just need to re-arrange all 
    # of the attributes based on the index we've constructed.

    if ($ascending) {

	# ProtStrings
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{ProtStrings}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{ProtStrings}}[$i] = $temp[$i];         }
	
	# ProtStarts
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{ProtStarts}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{ProtStarts}}[$i] = $temp[$i];         }
	
	# ProtEnds
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{ProtEnds}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{ProtEnds}}[$i] = $temp[$i];         }
	
	# StartOffsets
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{StartOffsets}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{StartOffsets}}[$i] = $temp[$i];         }
	
	# StartOffsetLengths
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{StartOffsetLengths}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{StartOffsetLengths}}[$i] = $temp[$i];         }
	
	# NuclStrings
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{NuclStrings}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{NuclStrings}}[$i] = $temp[$i];         }
	
	# NuclStarts
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{NuclStarts}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{NuclStarts}}[$i] = $temp[$i];         }
	
	# NuclEnds
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{NuclEnds}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{NuclEnds}}[$i] = $temp[$i];         }
	
	# EndOffsets
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{EndOffsets}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{EndOffsets}}[$i] = $temp[$i];         }
	
	# EndOffsetLengths
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{EndOffsetLengths}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{EndOffsetLengths}}[$i] = $temp[$i];         }
	
	# Scores
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{Scores}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{Scores}}[$i] = $temp[$i];         }
	
	# StopCodon
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{StopCodon}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{StopCodon}}[$i] = $temp[$i];         }
	
	# IsTerminal
	foreach $i (0..$numHits-1) { $temp[$i] = ${$self->{IsTerminal}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{IsTerminal}}[$i] = $temp[$i];         }
	
    } else { # descending order

	# ProtStrings
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{ProtStrings}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{ProtStrings}}[$i] = $temp[$i];                      }
	
	# ProtStarts
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{ProtStarts}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{ProtStarts}}[$i] = $temp[$i];                      }
	
	# ProtEnds
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{ProtEnds}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{ProtEnds}}[$i] = $temp[$i];                      }
	
	# StartOffsets
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{StartOffsets}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{StartOffsets}}[$i] = $temp[$i];                      }
	
	# StartOffsetLengths
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{StartOffsetLengths}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{StartOffsetLengths}}[$i] = $temp[$i];                      }
	
	# NuclStrings
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{NuclStrings}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{NuclStrings}}[$i] = $temp[$i];                      }
	
	# NuclStarts
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{NuclStarts}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{NuclStarts}}[$i] = $temp[$i];                      }
	
	# NuclEnds
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{NuclEnds}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{NuclEnds}}[$i] = $temp[$i];                      }
	
	# EndOffsets
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{EndOffsets}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{EndOffsets}}[$i] = $temp[$i];                      }
	
	# EndOffsetLengths
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{EndOffsetLengths}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{EndOffsetLengths}}[$i] = $temp[$i];                      }
	
	# Scores
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{Scores}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{Scores}}[$i] = $temp[$i];                      }
	
	# StopCodon
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{StopCodon}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{StopCodon}}[$i] = $temp[$i];                      }
	
	# IsTerminal
	foreach $i (0..$numHits-1) { $temp[($numHits-$i)-1] = ${$self->{IsTerminal}}[$index[$i]]; }
	foreach $i (0..$numHits-1) { ${$self->{IsTerminal}}[$i] = $temp[$i];                      }
	
    }

    #### DONE! ####
}


1;

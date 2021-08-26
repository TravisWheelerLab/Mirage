# Mirage

Multiple-sequence IsofoRm Alignment tool Guided by Exon boundaries.

## About

Mirage produces multiple sequence alignments of proteoforms by mapping sequences
back to their constitutive exons on the genome and then reducing those mappings
to alignments that preserve the proteoforms' exonic structures.

## Dependencies

Mirage has no build dependencies beyond the C standard library, but it does have
several runtime dependencies. These are described below. Executables should
exist somewhere in the user's `PATH`.

**blat** - <http://www.kentinformatics.com/index.html>

**hsi** - <https://github.com/TravisWheelerLab/hsi>

**spaln2** - <http://www.genome.ist.i.kyoto-u.ac.jp/~aln_user/spaln/>

**tblastn** - <https://blast.ncbi.nlm.nih.gov/Blast.cgi>

Appropriate versions of these are included in the `dependencies` directory of
the repository. The `SETUP.pl` script will build and install them.

## Build

Building Mirage requires a modern C compiler. Use the `CC` environment variable
to control the compiler used. Set the `WORK_DIR` environment variable to control
where the binaries end up, the default is `dist`.

Then, building is as simple as running `make`.

## Installation

Once the build is complete, run `make install` to copy binaries. Use the
`PREFIX` environment variable to control the install location. The default
location is `/usr/local`.


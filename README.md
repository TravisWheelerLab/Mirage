# Mirage

Multiple-sequence IsofoRm Alignment tool Guided by Exon boundaries.

## Dependencies

Mirage has no build dependencies beyond the C standard library, but it does have
several runtime dependencies. These are listed below. The build process will
include suitable versions of each in the build output, so there is no need to
install them manually.

**blat** - <http://www.kentinformatics.com/index.html>

**hsi** - <https://github.com/TravisWheelerLab/hsi>

**spaln2** - <http://www.genome.ist.i.kyoto-u.ac.jp/~aln_user/spaln/>

**tblastn** - <https://blast.ncbi.nlm.nih.gov/Blast.cgi>

## Build

We use CMake to build Mirage. The following commands will place everything
needed to run Mirage in the `build/` directory.

```
cmake .
make
```

## Installation

 The software can be "installed" simply by adding the `build/` directory to the
 `PATH`. Optionally, it can be relocated first (for example, to `/opt/`).

## Usage

TODO

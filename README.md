# Mirage

Multiple-sequence IsofoRm Alignment tool Guided by Exon boundaries.

## Build

Building Mirage requires a modern C compiler. Use the `CC` environment variable
to control the compiler used. Set the `WORK_DIR` environment variable to control
where the binaries end up, the default is `dist`.

Then, building is as simple as running `make`.

## Installation

Once the build is complete, run `make install` to copy binaries. Use the
`PREFIX` environment variable to control the install location. The default
location is `/usr/local`.


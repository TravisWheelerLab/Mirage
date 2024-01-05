# Mirage

Multiple-sequence IsofoRm Alignment tool Guided by Exon boundaries.

## Citation

[*Mirage2*’s high-quality spliced protein-to-genome mappings produce accurate multiple-sequence alignments of isoforms](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0285225).
Nord et al., 2023.


## About

Mirage produces multiple sequence alignments of proteoforms by mapping sequences
back to their constitutive exons on the genome and then reducing those mappings
to alignments that preserve the proteoforms' exonic structures.

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

```
./mirage2  [options]  proteins.fa  species-guide.txt
```

Mirage requires two inputs: a FASTA-formatted file of protein sequences and a
"species guide" file listing the paths to genome and, if available, GTF index
locations for each species.


**Protein Sequence File**

In order for Mirage to organize proteoform sequences according to their gene
families and species, we have developed the following |-delineated naming
convention for sequences in the protein sequence file:

```
>species|gene-family|id
```

The `species` and `gene-family` fields provide the names of the species and
gene family that each sequence belongs to, and `id` is a unique identifier
for the sequence within that species-family pair.
Additional information can be stored in the name field as a comment,
signified as anything content following a '#.'
For example, the following is a valid sequence name entry:

```
>human|obscn|isoform_1  # obscurin, primary isoform, Swiss Prot accession Q5VST9
```

Alternatively, protein sequences can be named following UniProt conventions,
where Mirage looks to the contents of the `OS` and `GN` fields to recognize the
sequence's species and gene family:

```
>sp|Q5VST9_iso1|OBSCN_HUMAN Obscurin OS=Homo_sapiens OX=9606 GN=OBSCN PE=1 SV=3
```

Because the simplified Mirage naming convention and the UniProt convention both
incorporate a triple of |-separated fields, it is critical to preserve the `OS`
and `GN` fields in sequences intended to be parsed under the UniProt convention.
In the above example, removing those fields would cause Mirage to mistakenly
identify the sequence as belonging to a species named 'sp' and a gene family
named 'Q5VST9_iso1.'


**Species Guide File**

The species guide provides whitespace-delineated mappings of species to their
local genome and (optionally) GTF file locations.
A '-' in a species' GTF field indicates that a GTF index has not been provided.
Species that are not listed in the species guide file will be integrated into
alignments using a Needleman-Wunsch-style dynamic programming method rather than
*via* genome mapping.
Finally, a Newick-formatted species tree can be provided on the first line of
the species guide file to guide the order by which intra-species alignments are
merged.
The following example is a valid species guide file:

```
(human,(mouse,rat)),magpie
human  ~/genomes/human.fa  ~/gtfs/human.gtf
mouse  ~/genomes/mouse.fa  -
magpie ~/genomes/magpie.fa -
rat    ~/genomes/rat.fa    ~/gtfs/rat.gtf
```

In order to facilitate genome and GTF downloading and organization, the
'DownloadGenomicData.pl' script can be used to download the most up-to-date
genomes and associated GTF indices for species represented in the
UCSC Genome Browser.
The input to this script is a text file where each line lists a desired
species, for example:

```
human
mouse
golden eagle
```

The script can simply be called with the following command:

```
./build/DownloadGenomicData.pl  species-list.txt
```


**Options**

`-cpus n`

Run Mirage using n cpu cores.

`-outdirname my-results`

Set the name of the Mirage results directory to be 'my-results.'

`--time`

Print robust timing data

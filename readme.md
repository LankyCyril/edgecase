edgeCase
========

*edgeCase* is a framework for extraction and interpretation of telomeric reads
from long-read single-molecule whole genome sequencing datasets. Associated
preprint: https://www.biorxiv.org/content/10.1101/2020.01.31.929307v1

![densityplot_sample](assets/densityplot-haplotypes.png?raw=true "densityplot example")


## Installation

#### With Conda (preferred):

```{sh}
$ git clone https://github.com/LankyCyril/edgecase
$ cd edgecase
$ conda env create --name edgecase --file environment.yaml
$ conda activate edgecase
$ ./edgecase
```

#### By manually installing dependencies:

```{sh}
$ git clone https://github.com/LankyCyril/edgecase
$ cd edgecase
$ pip install numpy scipy scikit-learn statsmodels numba
$ pip install pandas matplotlib seaborn tqdm regex pysam
$ ./edgecase
```


## Input data and formats

### The extended reference genome and BAM files

*edgeCase* works with SAM/BAM files aligned to a reference that is annotated
with known subtelomeric regions and uses reads anchored to the outermost ends
of subtelomeres (5' on the *p* arm, 3' on the *q* arm). For a BAM file
*dataset.bam* aligned to *ref.fa*, it needs several files:

```
dataset.bam.bai: a BAI index; create with "samtools index dataset.bam"
ref.fa.fai:      a FAI index; create with "samtools faidx ref.fa"
ref.fa.ecx:      an index containing annotations of subtelomere-telomere boundareis
```

*ref.fa.ecx*, a.k.a. the edgeCase indeX, describes anchors of interest in the
reference genome; the format is based on the BED format. Usable "flag" values
*have* to be among 4096 (hard mask), 8192 (fork), 16384 (telomeric tract). Two
examples of ECX files can be found in the "assets" subdirectory.

Specifically, as described in the bioRxiv preprint, the human reference can be
constructed from the hg38/GRCh38 reference genome and subtelomeric assemblies
published by [Stong et al., 2014](https://dx.doi.org/10.1101%2Fgr.166983.113).
To generate this reference, which we call "extended", or *hg38ext*, run
`tools/generate-hg38ext.py --remote > hg38ext.fa`.


### Custom SAM flags

*edgeCase* extends the zoo of SAM flags with four of its own. The full table of
flag names:  
(also see https://broadinstitute.github.io/picard/explain-flags.html)

name               | value | hex value | comment
-------------------|-------|-----------|----------------------------------------------
paired             | 1     | 0x0001    | SAM specification flag
mapped_proper_pair | 2     | 0x0002    | SAM specification flag
unmapped           | 4     | 0x0004    | SAM specification flag
mate_unmapped      | 8     | 0x0008    | SAM specification flag
rev                | 16    | 0x0010    | SAM specification flag
mate_rev           | 32    | 0x0020    | SAM specification flag
1stmate            | 64    | 0x0040    | SAM specification flag
2ndmate            | 128   | 0x0080    | SAM specification flag
secondary          | 256   | 0x0100    | SAM specification flag
qcfail             | 512   | 0x0200    | SAM specification flag
pcrdup             | 1024  | 0x0400    | SAM specification flag
supp               | 2048  | 0x0800    | SAM specification flag
ucsc_mask_anchor   | 4096  | 0x1000    | edgeCase-specific flag; added during pipeline
fork               | 8192  | 0x2000    | edgeCase-specific flag; added during pipeline
tract_anchor       | 16384 | 0x4000    | edgeCase-specific flag; added during pipeline
is_q               | 32768 | 0x8000    | edgeCase-specific flag; added during pipeline

*Note:* All edgeCase routines that allow flag filtering recognize both the
numeric flag format (such as 3844) and the "human-readable" format such as "rev"
or "is_q|paired". Combinations are also understood, for example, "3844|is_q".


## The edgeCase pipeline

```
usage: ./edgecase [-h] {tailpuller,tailchopper,repeatfinder,kmerscanner,levenshtein,densityplot} ...

positional arguments:
    tailpuller                select overhanging reads
    tailchopper               get overhanging heads/tails of reads
    repeatfinder              discover enriched repeats in candidate sequences
    kmerscanner               perform scan of known kmers/motifs
    levenshtein               cluster reads by edit distance
    densityplot               visualize densities of candidate reads
```


### tailpuller

```
usage: ./edgecase tailpuller --index X [options] bam > sam

positional arguments:
  bam                        name of input BAM/SAM file

optional arguments:
  -x X, --index X            location of the reference .ecx index (REQUIRED)
  -f f, --flags f            process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g        process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F      process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q      process only entries with MAPQ >= Q (default: 0)
  -m M, --max-read-length M  max read length to consider when selecting lookup regions (default: None)
```

Outputs a subset SAM file that contains only the reads that overhang anchors
defined in the ECX. If the read overhangs the mask anchor, the 4096 SAM flag is
added; for forks, 8192 is added; for telomeric tracts, 16384.  
For reads on the *q* arm (i.e., on the 3' end), the 32768 flag is added (see
above for the full list and the explanation of flags).

**NB**: these flags are unused in the SAM specification and should not clash with
anything. `samtools view` can correctly subset using these flags.

Suggestions:
* use `-F 3844` to skip secondary, supplementary and QC-fail alignments;
* pipe the output through `samtools view -bh -` to compress on the fly;
* supplying `--max-read-length` drastically improves wall time if reads are
significantly shorter than chromosomes.


### tailchopper

```
usage: ./edgecase tailchopper --index X [options] bam > fasta

positional arguments:
  bam                    name of input BAM/SAM file

optional arguments:
  -x X, --index X        location of the reference .ecx index (REQUIRED)
  -f f, --flags f        process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g    process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F  process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q  process only entries with MAPQ >= Q (default: 0)
  -t ?, --target ?       an ECX flag (cut relative to reference) or 'cigar' (default: tract_anchor)
```

Truncates reads in the tailpuller file either to soft/hard-clipped ends (when
--target is "cigar"), or to sequences extending past given anchor (when
--target is "tract_anchor", "fork", or "ucsc_mask_anchor").

**NB**: outputs a SAM file with unmapped reads (sets the 0x0004 bit in the
flag), but *retains the original mapping position*; do *not* use this value for
downstream analyses unless you know exactly what you are after.


### repeatfinder

```
usage: ./edgecase repeatfinder [options] sequencefile > tsv

positional arguments:
  sequencefile                name of input SAM/BAM/FASTA/FASTQ file

optional arguments:
  -f f, --flags f             process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g         process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F       process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q       process only entries with MAPQ >= Q (default: 0)
  --fmt ?                     format of input file(s) (default: sam)
  -m ?, --min-k ?             smallest target repeat length (default: 4)
  -M ?, --max-k ?             largest target repeat length (default: 16)
  -r R, --min-repeats R       minimum number of consecutive repeats (default: 2)
  -n ?, --max-motifs ?        maximum number of motifs to report (default: None)
  -P ?, --max-p-adjusted ?    cutoff adjusted p-value (default: 0.05)
  --jellyfish ?               jellyfish binary (unless in $PATH) (default: None)
  --jellyfish-hash-size ?     jellyfish initial hash size (default: 2G)
  -j J, --jobs J              number of jellyfish jobs (parallel threads) (default: 1)
```

Performs Fisher's exact tests on *k*-mer counts to identify significantly
enriched repeating motifs of lengths from *--min-k* to *--max-k* in the input
file.  
Relies on [jellyfish](http://www.genome.umd.edu/jellyfish.html) to count
*k*-mers. If *edgeCase* has been installed with the Conda method (by creating
an environment from *environment.yaml*), *jellyfish* is already installed and no
special action is needed. Otherwise, it needs to be installed manually and, if
not in $PATH, supplied with the *--jellyfish* option.


### kmerscanner

```
usage: ./edgecase kmerscanner --motif-file M [options] bam > dat

positional arguments:
  bam                    name of input SAM/BAM file

optional arguments:
  --motif-file M         file with repeated motif sequences (REQUIRED)
  -f f, --flags f        process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g    process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F  process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q  process only entries with MAPQ >= Q (default: 0)
  -w W, --window-size W  size of the rolling window (default: 100)
  --head-test H          length of head to use for density filter (if specified) (default: None)
  --tail-test T          length of tail to use for density filter (if specified) (default: None)
  -c C, --cutoff C       use hard cutoff for density (default: None)
  -j J, --jobs J         number of jobs to run in parallel (default: 1)
```

In a rolling window along each read in a BAM file, calculates densities of given
motifs and outputs a DAT file.  
Optionally filters input by terminal density (outputs data only for reads
exceeding density cutoff). By default, outputs data for all input reads.  
*--motif-file* is usually the output of *repeatfinder*, but can be an arbitrary
tab-separated file where the first field of each line is a motif (except for
lines starting with "#" which are treated as comments).

**NB**: it is possible to run *kmerscanner* on an entire WGS BAM and look for
reads that pass tests encoded by *--head-test*, *--tail-test*, and *--cutoff*,
but this use case is experimental and is discouraged.  Generally, *kmerscanner*
is intended to be used downstream of *tailpuller/tailchopper* and *repeatfinder*
to calculate densities of identified motifs in telomeric candidate reads. In
this case, *--head-test*, *--tail-test*, and *--cutoff* should be omitted.


### levenshtein

```
usage: ./edgecase levenshtein [options] sequencedata > tsv

positional arguments:
  sequencedata           name of input BAM/SAM file or directory with precomputed distances

optional arguments:
  -f f, --flags f        process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g    process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F  process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q  process only entries with MAPQ >= Q (default: 0)
  --kmerscanner-file ?   kmerscanner file (optional, for use with --output-dir)
  --min-cluster-size ?   minimum cluster size to consider (default: 5)
  -o ?, --output-dir ?   output directory for clustermaps and per-haplotype SAM files (default: None)
```

For each chromosome arm in the BAM/SAM file, clusters the reads that align there
by their relative pairwise edit distance (Levenshtein distance) and decides on
the best number of clusters by maximizing the Bayesian information criterion.
If more than one cluster is identified, performs a Mann-Whitney U one-tailed
test on all intra-cluster distances vs. all inter-cluster distances.  
If the input is a directory with precomputed matrices (files matching mask
*${chromosome_name}-matrix.tsv*), uses these values to cluster and compute
*p*-values (skips the actual step of distance computation).  
If *--kmerscanner-file* is provided, generates kmerscanner files for read
clusters (haplotypes) on each arm where more than one such cluster is
detected.

**NB**: this is an experimental module, and the maximum number of outliers is
hard-coded as 1. This worked for the datasets analyzed in the bioRxiv preprint,
but the number of outliers may have to be adjusted for other datasets.  
**NB**: this algorithm scales quadratically with the number of input reads and
is computationally infeasible for large datasets.


### densityplot

```
usage: ./edgecase densityplot --index X [options] dat > pdf

positional arguments:
  dat                    input density file

optional arguments:
  -x X, --index X        location of the reference .ecx index (REQUIRED)
  -f f, --flags f        process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g    process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F  process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q  process only entries with MAPQ >= Q (default: 0)
  -z, --gzipped          input is gzipped (must specify if any of -qfF present) (default: False)
  -b B, --bin-size B     size of each bin in bp for visualization speedup (default: 100)
  --zoomed-in            plot taller traces, cut off pre-anchor regions (default: False)
  --palette ?            custom palette for plotting motifs (default: None)
  -e, --exploded         plot each read separately (default: False)
  --title T              figure title (defaults to input filename) (default: None)
```

Visualizes the density of motifs on each chromosomal arm in the regions covered
by candidate reads, binning the values by windows of *--bin-size*.  
The value of *--palette* can be either none (in which case the maximum of nine
motifs can be plotted with default colors), "paper" or "paper|legend=False" (in
which case motifs known from research can be plotted with custom colors,
matching the colors in the figures in the bioRxiv preprint), or a chained
key-value sequence of "motif=color" and "legend=boolean", for example:
"TTAGGG=green|TGAGGG=#D01000|legend=True".  
Options *--exploded* and *--title* are deprecated.

Option *--zoomed-in* plots taller figures, discards non-telomeric regions, and
visualizes read coverage above each plot. With this option, two custom "debug"
environment variables can be passed to *densityplot* that specify how much
of the surrounding reference coordinates should be included: *PAPER_LEFT_SPAN*
and *PAPER_RIGHT_SPAN*.

Annotates the anchors from the ECX with dashed lines:
* ucsc_mask_anchor == gray,
* fork == blueviolet,
* tract_anchor == red.

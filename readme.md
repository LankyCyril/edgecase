edgeCase
========

edgeCase is a framework for extraction and interpretation of telomeric reads
from long-read single-molecule whole genome sequencing datasets.

![densityplot_sample](assets/densityplot-exploded.png?raw=true "densityplot example")

## Installation

This tool is in active development, so it is not pip-installable yet.
The installation options are:

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

## Overview

edgeCase is a small toolchain that takes a set of aligned reads as the initial
input and describes its telomeric content.

```{sh}
usage: ./edgecase [-h] [-j J] {tailpuller,tailchopper,kmerscanner,densityplot,assembler} ...

positional arguments:
    tailpuller       select overhanging reads
    tailchopper      get overhanging heads/tails of reads
    kmerscanner      perform kmer scan
    densityplot      visualize densities of candidate reads
    assembler        assemble haplotypes at ends of chromosomes

optional arguments:
  -j J, --jobs J     number of jobs to run in parallel (default: 1)
```

Notes:
* The `--jobs` option currently only has effect for `kmerscanner`.
* The `assembler` routine is not available yet (in development).

### ./edgecase tailpuller [options] bam > sam

```{sh}
positional arguments:
  bam                        name of input BAM file

optional arguments:
  -x X, --index X            location of the reference .ecx index
  -F F, --flag-filter F      process only sam entries with none of these flags present (default: 0)
  -m M, --max-read-length M  max read length to consider when selecting lookup regions (default: None)
```

The input files are:
* ECX: a.k.a. the edgeCase indeX, describing anchors of interest in the
reference genome; the format is based on the BED format. Usable 'flag' values
*have* to be among 4096 (hard mask), 8192 (fork), 16384 (telomeric tract). Two
examples of ECX files can be found in the `assets/` subdirectory.
* BAM: reads aligned to the reference genome. Has to have a .bai index.

Outputs a subset SAM file that contains only the reads that overhang anchors
defined in the ECX. If the read overhangs the mask anchor, the 4096 SAM flag is
added; for forks, 8192 is added; for telomeric tracts, 16384.  
For reads on the q arm (i.e., on the 3' end), the 32768 flag is added.  
NB: these flags are unused in the SAM specification and should not clash with
anything. `samtools view` can correctly subset using these flags.

Suggestions:
* use `-F 3844` to skip secondary, supplementary and QC-fail alignments;
* pipe the output through `samtools view -bh -` to compress on the fly;
* supplying `--max-read-length` drastically improves wall time if reads are
significantly shorter than chromosomes.

Bells and whistles:
* **All** edgeCase routines that allow flag filtering (tailpuller, tailchopper,
densityplot) recognize both the numeric flag format (such as 3844) and the
"human-readable" format such as "rev" or "is_q|paired". Combinations are also
understood, for example, "3844|rev".

The full table of flag names (also see
https://broadinstitute.github.io/picard/explain-flags.html):

name               | value | comment
-------------------|-------|---------------------------
paired             | 1     |
mapped_proper_pair | 2     |
unmapped           | 4     |
mate_unmapped      | 8     |
rev                | 16    |
mate_rev           | 32    |
1stmate            | 64    |
2ndmate            | 128   |
secondary          | 256   |
qcfail             | 512   |
pcrdup             | 1024  |
supp               | 2048  |
ucsc_mask_anchor   | 4096  | available after tailpuller
fork               | 8192  | available after tailpuller
tract_anchor       | 16384 | available after tailpuller
is_q               | 32768 | available after tailpuller

### ./edgecase tailchopper [options] bam > fasta

```{sh}
positional arguments:
  bams                        name of input BAM file

optional arguments:
  -x X, --index X             location of the reference .ecx index
  -f f, --flags f             process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g         process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F       process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q       process only entries with MAPQ >= Q (default: 0)
  -t {ucsc_mask_anchor,tract_anchor,cigar,fork}, --target {ucsc_mask_anchor,tract_anchor,cigar,fork}
                              either an ECX flag (cut relative to reference) or cigar (cut clipped ends)
                              (default: tract_anchor)
```

Truncates reads in the tailpuller file to a FASTQ of sequences of
soft/hard-clipped ends or ends overhanging given anchors.

### ./edgecase [-j J] kmerscanner [options] bams > dat

```{sh}
positional arguments:
  bams                   name(s) of input BAM/SAM file(s)

optional arguments:
  --motif M              target motif sequence (default: TTAGGG)
  --head-test H          length of head to use for density filter (if specified) (default: None)
  --tail-test T          length of tail to use for density filter (if specified) (default: None)
  -c C, --cutoff C       use hard cutoff for density (default: None)
  -w W, --window-size W  size of the rolling window (default: 120)
  -n N, --num-reads N    expected number of reads in input (for progress display) (default: None)
```

In a rolling window along each read in a BAM file, calculates densities of given
motifs and outputs a DAT file.  
Optionally filters input by terminal density (outputs data only for reads
exceeding density cutoff).  
By default, outputs data for all input reads.

The value of `--motif` is technically a very simple regex that allows symbols
'A', 'C', 'G', 'T', '.', and '|'. Examples of valid values are: 'TTAGGG',
'TT.GGG|CCC.AA'.

### ./edgecase densityplot [options] dat > pdf

```
positional arguments:
  dat                    input density file

optional arguments:
  -z, --gzipped          input is gzipped (must specify if any of -qfF present) (default: False)
  -x X, --index X        location of the reference .ecx index (default: no default, required)
  -f f, --flags f        process only entries with all these sam flags present (default: 0)
  -g g, --flags-any g    process only entries with any of these sam flags present (default: 65535)
  -F F, --flag-filter F  process only entries with none of these sam flags present (default: 0)
  -q Q, --min-quality Q  process only entries with MAPQ >= Q (default: 0)
  -b B, --bin-size B     size of each bin in bp for visualization speedup (default: 100)
  --title T              figure title (defaults to input filename) (default: None)
  --no-align             plot unaligned (default: False)
```

Visualizes the density of motifs in each read, placing them at their mapping
positions on the reference.  
Annotates the anchors from the ECX with dashed lines:
* hard mask anchor == gray
* fork == red
* telomeric tract == green

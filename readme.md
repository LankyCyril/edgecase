edgeCase
========

*edgeCase* is a framework for extraction and interpretation of telomeric reads
from long-read single-molecule whole genome sequencing datasets. Associated
preprint: https://www.biorxiv.org/content/10.1101/2020.01.31.929307v1

![haplotypes_example](assets/haplotypes-example.png?raw=true "haplotypes example")


## Installation

### Obtaining code

The code can either be downloaded from the
[releases](https://github.com/LankyCyril/edgecase/releases) page
(the [SLIM](https://github.com/LankyCyril/edgecase/releases/download/2021.02.27/edgecase-20210227-SLIM.tar.gz) tarball),  
or cloned with git: `git clone https://github.com/LankyCyril/edgecase`

### Environment setup

#### With Conda (preferred)

```{sh}
$ cd edgecase
$ conda env create --name edgecase --file environment.yaml
$ conda activate edgecase
$ ./edgecase
```

#### By manually installing dependencies

```{sh}
$ cd edgecase
$ pip install numpy scipy scikit-learn statsmodels edlib
$ pip install pandas matplotlib seaborn tqdm regex pysam
$ ./edgecase
```


## Input data and formats

### The extended reference genome

*edgeCase* works with SAM/BAM files aligned to a reference that is annotated
with known subtelomeric regions and uses reads anchored to the outermost ends
of subtelomeres (5' on the *p* arm, 3' on the *q* arm).

```
ref.fa.fai:  a FAI index; create with "samtools faidx ref.fa"
ref.fa.ecx:  an index containing annotations of subtelomere-telomere boundaries
```

*ref.fa.ecx*, a.k.a. the edgeCase indeX, describes anchors of interest in the
reference genome; the format is based on the BED format. Usable "flag" values
*have* to be among 4096 (hard mask), 8192 (fork), 16384 (telomeric tract). Two
examples of ECX files can be found in the "assets" subdirectory.

Specifically, as described in the bioRxiv preprint, the human reference can be
constructed from the hg38/GRCh38 reference genome and subtelomeric assemblies
published by [Stong et al., 2014](https://dx.doi.org/10.1101%2Fgr.166983.113).
To generate this reference, which we call "extended", or *hg38ext*, run
`assets/generate-hg38ext.py --remote > hg38ext.fa`.

### Alignment files

We recommend using *minimap2* to generate BAM files for edgeCase. Another option
is *winnowmap*, but it has not been sufficiently tested yet.

**NB**: currently, it is imperative to supply a BAM file where secondary
alignment entries have read sequences. For example, *minimap2* creates BAMs in
this format with the use of the `-Y` switch.  
We plan to implement a workaround for this requirement in the future.

BAM files must also be indexed (i.e., have a `.bai` file created with
`samtools index`).


## Custom SAM flags

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
mask_anchor        | 4096  | 0x1000    | edgeCase-specific flag; added during pipeline
fork               | 8192  | 0x2000    | edgeCase-specific flag; added during pipeline
tract_anchor       | 16384 | 0x4000    | edgeCase-specific flag; added during pipeline
is_q               | 32768 | 0x8000    | edgeCase-specific flag; added during pipeline

**NB**: these flags are unused in the SAM specification and should not clash with
anything. `samtools view` can correctly subset using these flags.

*Note:* All edgeCase routines that allow flag filtering recognize both the
numeric flag format (such as 3844) and the "human-readable" format such as "rev".
Combinations are also understood, for example, "-F 3844 -F is_q".

*Note:* In the future, custom SAM flags may be superseded with tags.  
A backwards compatibility layer will be provided (i.e., arguments like "-f fork"
or "-F 16384" will still work but interpret and produce appropriate tags).


## The edgeCase pipeline

```
Usage: ./edgecase [-h | --help]
       ./edgecase <command> [<args>...]

Commmands (<command>):
    tailpuller               select overhanging long reads
    tailchopper              get overhanging heads/tails of long reads
    repeatfinder             discover enriched repeats in candidate sequences
    kmerscanner              perform scan of known kmers/motifs
    densityplot              visualize densities of candidate motifs

Development area:
    entropy                  calculate motif entropy among long reads
    levenshtein              calculate pairwise edit distance among long reads
```

All commands output their results to stdout; you must pipe them into other
commands or into the destination file. This applies even to outputs in PDF and
PKL formats.

**NB**: Depending on the aligner used upstream, MAPQ of secondary reads may have
been set to zero regardless of real mapping quality; use this filtering option
with caution. **This warning applies to all edgeCase subroutines that accept
the `-q` filtering flag.**


### tailpuller

Outputs a subset SAM file that contains only the reads that overhang anchors
defined in the ECX. If the read overhangs the mask anchor, the 4096 SAM flag is
added; for forks, 8192 is added; for telomeric tracts, 16384.  
For reads on the *q* arm (i.e., on the 3' end), the 32768 flag is added (see
above for the full list and the explanation of flags).

```
Usage: ./edgecase tailpuller -x filename [-t targetspec]...
                            [-M integer] [--min-map-overlap integer]
                            [-m integer] [--min-telomere-overlap integer]
                            [--output-ambiguous-reads string]
                            [-f flagspec]... [-F flagspec]... [-q integer] <bam>

Required options:
    -x, --index [filename]                   location of the reference .ecx index

Options:
    -t, --target [targetspec]                target reads overlapping these features (ECX flags) [default: tract_anchor]
    -M, --max-read-length [integer]          maximum read length to consider when selecting lookup regions
    --min-map-overlap [integer]              minimum overlap of reference to consider read as mapped [default: 1]
    -m, --min-subtelomere-overlap [integer]  minimum overlap of subtelomere to consider read as candidate [default: 1]
    --min-telomere-overlap [integer]         minimum overlap of telomere to consider read as candidate [default: 1]
    --output-ambiguous-reads [string]        which ambiguously mapping reads to retain (none, all, longest-overlap) [default: none]

Input filtering options:
    -f, --flags [flagspec]                   process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]             process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]              process only entries with this MAPQ or higher [default: 0]
```

Suggestions:
* It is recommended to include secondary and supplementary reads (i.e., leave
  the -F flag as default [0]), because:
    * edgeCase determines unambiguously mapped reads on its own; aligners
      assign the 'supplementary' flag to multi-mapping reads arbitrarily, and
      removing such supplementary reads upstream may lead to loss of information
      in telomeric regions;
    * edgeCase will discard chimeric reads in terminal regions if information
      about supplementary alignments is present.
* Supplying `--max-read-length` drastically improves wall time if reads are
  significantly shorter than chromosomes; for PacBio HiFi (CCS) it is suggested
  to use the value of 30000. If the value is not specified, edgeCase will
  assume *infinity*, and will have to go over the entire content of the BAM file.
* Suggested value of --min-map-overlap for PacBio HiFi: 500.
* Suggested value of --min-(sub)telomere-overlap for PacBio HiFi: 3000.
* Pipe the output through `samtools view -bh -` to compress on the fly.


### tailchopper

Truncates reads in the tailpuller file either to soft/hard-clipped ends (when
--target is "cigar"), or to sequences extending past given anchor (when
--target is "tract_anchor", "fork", or "mask_anchor").

Outputs a SAM file with overhanging tails of candidate reads.

```
Usage: ./edgecase tailchopper -x filename [-t targetspec]
                            [-f flagspec]... [-F flagspec]... [-q integer] <bam>

Required options:
    -x, --index [filename]        location of the reference .ecx index

Options:
    -t, --target [targetspec]     an ECX flag (cut relative to reference) or 'cigar' [default: tract_anchor]

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
```

**NB**: tailchopper outputs a SAM file with unmapped reads (sets the 0x0004 bit
in the flag), but *retains the original mapping position*; do *not* use this
value for downstream analyses unless you know exactly what you are after.

*Suggestion*: pipe the output through `samtools view -bh -` to compress on the
fly.


### repeatfinder

Expects the SAM/BAM file from `tailchopper` as input; however, will also work
on any SAM/BAM file as well as Fasta/Fastq files.

Performs Fisher's exact tests on *k*-mer counts to identify significantly
enriched repeating motifs of lengths from `--min-k` to `--max-k` in the input
file.  
Relies on [jellyfish](http://www.genome.umd.edu/jellyfish.html) to count
*k*-mers. If edgeCase has been installed with the Conda method (by creating
an environment from `environment.yaml`), `jellyfish` is already installed and no
special action is needed. Otherwise, it needs to be installed manually and, if
not in `$PATH`, supplied with the `--jellyfish` option.

Outputs a TSV file with columns:  
`monomer motif length score fraction_explained p p_adjusted`

```
Usage: ./edgecase repeatfinder [-m integer] [-M integer] [-r integer] [-P float]
                               [--jellyfish filename] [--jellyfish-hash-size string]
                               [-n integer] [-j integer] [-q integer]
                               [-f flagspec]... [-F flagspec]... [--fmt string]
                               [--collapse-reverse-complement] <sequencefile>

Options:
    --fmt sam|fastx                     format of input file [default: sam]
    -m, --min-k [integer]               smallest target repeat length [default: 4]
    -M, --max-k [integer]               largest target repeat length [default: 16]
    -r, --min-repeats [integer]         minimum number of consecutive repeats [default: 2]
    -P, --max-p-adjusted [float]        cutoff adjusted p-value [default: .05]
    --jellyfish [filename]              jellyfish binary (unless in $PATH)
    -s, --jellyfish-hash-size [string]  jellyfish initial hash size [default: 2G]
    -n, --max-motifs [integer]          maximum number of motifs to report
    -j, --jobs [integer]                number of jellyfish jobs (parallel threads) [default: 1]
    -C, --collapse-reverse-complement   collapse counts of reverse complement motifs

Input filtering options:
    -f, --flags [flagspec]              process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]        process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]         process only entries with this MAPQ or higher [default: 0]
```


### kmerscanner

Expects the SAM/BAM file from `tailpuller` as input; however, will also work
on any SAM/BAM file as well as Fasta/Fastq files.

Expects the TSV file from `kmerscanner` provided as the `--motif-file` option;
however, one may supply an arbitrary tab-separated file where the first field of
each line is a motif (except for lines starting with "#" which are treated as
comments).

In a rolling window along each read in a BAM file, calculates densities of given
motifs and outputs a tab-separated DAT file with columns:  
`name flag chrom pos mapq motif score clip_5prime clip_3prime b=N`,  
where the last column name contains the value of `--bin-size`, and the column
itself lists all density values along rolling windows for a given motif.

*Note*: it is recommended to pipe the output through `gzip`, as these files
are quite verbose and easily compressible. In the future, we plan to implement
a more space-efficient (and backwards-compatible) format.

```
Usage: ./edgecase kmerscanner [-j integer] --motif-file filename
                              [-b integer] [-n integer]
                              [-f flagspec]... [-F flagspec]... [-q integer]
                              [--fmt string] <sequencefile>

Required options:
    --motif-file [filename]       file with repeated motif sequences (output of `repeatfinder`)

Options:
    --fmt sam|fastx               format of input file [default: sam]
    -b, --bin-size [integer]      size of the rolling window [default: 10]
    -n, --num-reads [integer]     expected number of reads in input (for progress display)
    -j, --jobs [integer]          number of jobs to run in parallel [default: 1]

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
```


### densityplot

Expects the DAT file from `kmerscanner` as input;  
visualizes the density of motifs on each chromosomal arm in the regions covered
by candidate reads.

The value of `--palette` can be either none (in which case the maximum of nine
motifs can be plotted with default colors), "paper", "paper|legend=full",
"paper|legend=density", or "paper|legend=motifs" (in which case motifs known
from research can be plotted with custom colors, matching the colors in the
figures in the paper), or a chained key-value sequence of "motif=color" and
"legend=spec", where "spec" is one of "none", "full", "density", "motifs".
For example: `"TTAGGG=green|TGAGGG=#D01000|legend=full"`.

Annotates the anchors from the ECX with dashed lines:
* mask_anchor == gray,
* fork == blueviolet,
* tract_anchor == red.

Outputs a PDF file (writes it to stdout; you must pipe the output into a file).
Alternatively, can output a Python pickle file (with `--outfmt=pkl`).

```
Usage: ./edgecase densityplot -x filename [-b integer] [--plot-coverage]
                              [--palette palettespec] [--title string]
                              [--n-boot integer] [--chroms-to-plot string]
                              [-f flagspec]... [-F flagspec]... [-q integer]
                              [--figwidth-inches float] [--outfmt string] [-z] <dat>

Required options:
    -x, --index [filename]        location of the reference .ecx index

Options:
    -z, --gzipped                 input is gzipped (must specify if any of -qfF present)
    -b, --bin-size [integer]      size of each bin in bp (overrides bin size in <dat>)
    --n-boot [integer]            number of bootstrap iterations for 95% confidence intervals [default: 1000]
    --palette [palettespec]       custom palette for plotting motifs
    --title [string]              figure title (defaults to input filename)
    --chroms-to-plot [string]     if set, plot chromosomes from this comma-separated list unconditionally
    --plot-coverage               plot coverage by telomeric reads on each arm
    --figwidth-inches [float]     width of figure in inches [default: 13]
    --outfmt [string]             output format (pdf, pkl) [default: pdf]

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
```


### entropy

Expects the DAT file from `kmerscanner` as input; can accept multiple DAT files
at once.

Calculates entropy values of motif assignments per window among reads,
and outputs a TSV file with columns:  
`entropy coverage`.

```
Usage: ./edgecase entropy [-b integer] [-f flagspec]... [-F flagspec]... [-q integer]
                          [-z] <dat>...

Options:
    -z, --gzipped                 input is gzipped (must specify if any of -qfF present)
    -b, --bin-size [integer]      size of each bin in bp (overrides bin size in <dat>)

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
```


### levenshtein

Expects the SAM/BAM file from `tailpuller` as input.

Calculates pairwise relative edit distance (Levenshtein distance) for all pairs
of reads mapping to each chromosomal arm in the input SAM/BAM file.  
Outputs a TSV file with columns:  
`rname qname1 qname2 relative_ld`,  
where `rname` is the name of the chromosome, `qnameN` is the name of a read in
the pair, and `relative_ld` is the distance.

**NB**: this algorithm scales quadratically with the number of input reads and
is computationally infeasible for large datasets.

```
Usage: ./edgecase levenshtein [-f flagspec]... [-F flagspec]... [-q integer]
                              [-j integer] <sequencefile>

Options:
    -j, --jobs [integer]               number of jobs to run in parallel [default: 1]

Input filtering options:
    -f, --flags [flagspec]             process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]       process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]        process only entries with this MAPQ or higher [default: 0]
```

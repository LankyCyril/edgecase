edgeCase
========

From a set of aligned reads (BAM/SAM file) and a reference, edgeCase currently can create:
* 5AC and 3AC SAM files: subsets of reads mapped to ends of main chromosomes and extending into terminal hard-masked regions
* 5OOB and 3OOB FASTQ files: parts of 5AC and 3AC reads completely in the masked region
* density DAT files: densities of given motifs in a rolling window along each read in a set
* density PDF files: densities of motifs, visualized with respect to alignment

```
usage: ./edgecase [-h] [-j J] {kmerscanner,tailpuller,tailchopper,densityplot} ...

positional arguments:
    tailpuller        select overhanging reads
    tailchopper       get clipped heads/tails of reads
    kmerscanner       perform kmer scan
    densityplot       visualize densities of candidate reads

optional arguments:
  -h, --help          show this help message and exit
  -j J, --jobs J      number of jobs to run in parallel (default: 1)
```

### tailpuller

Subsets a sorted BAM/SAM file to either 5AC or 3AC SAM file.

```
usage: ./edgecase tailpuller [options] bams > sam

positional arguments:
  bams                     name(s) of input BAM/SAM file(s)

optional arguments:
  -h, --help               show this help message and exit
  -r R, --reference R      reference FASTA (default: no default, required)
  -p {5,3}, --prime {5,3}  which 'prime' end to output (default: 5)
```

### tailchopper

Truncates reads in a 5AC or 3AC file to a FASTQ of sequences completely overhanging reference hard mask.  
NB! Needs `--prime` to be set, as it has no knowledge of how the SAM was generated.

```
usage: ./edgecase tailchopper [options] bams > fasta

positional arguments:
  bams                     name(s) of input BAM/SAM file(s)

optional arguments:
  -h, --help               show this help message and exit
  -p {5,3}, --prime {5,3}  which 'prime' end to output (default: 5)
```

### kmerscanner

In a rolling window along each read in a SAM file, calculates densities of given motifs and outputs a DAT file.  
Optionally filters input by terminal density (outputs data only for reads exceeding density cutoff).  
By default, outputs data for all input reads.

```
usage: ./edgecase [-j J] kmerscanner [options] bams > dat

positional arguments:
  bams                   name(s) of input BAM/SAM file(s)

optional arguments:
  -h, --help             show this help message and exit
  --motif M              target motif sequence (default: TTAGGG)
  --head-test H          length of head to use for density filter (if specified) (default: None)
  --tail-test T          length of tail to use for density filter (if specified) (default: None)
  -c C, --cutoff C       use hard cutoff for density (default: None)
  -w W, --window-size W  size of the rolling window (default: 120)
  -n N, --num-reads N    expected number of reads in input (for progress display) (default: None)
```

### densityplot

Visualizes density data of reads, placing them to their mapping positions on the reference.

```
usage: ./edgecase densityplot [options] dat > png

positional arguments:
  dat                 input density file

optional arguments:
  -h, --help          show this help message and exit
  -b B, --bin-size B  size of each bin in bp for visualization speedup (default: 100)
  --title TITLE       figure title (defaults to input filename) (default: None)
```

from sys import stdout, stderr
from collections import OrderedDict
from edgecaselib.kmerscanner import get_circular_pattern
from pysam import FastxFile
from numpy import array

__doc__ = """edgeCase shortread: experiments with short reads

Usage: {0} shortread [-m integer] [-M integer] [-r integer]
       {1}           [--motifs string] [--fmt string] <sequencefile>

Output:
    TSV-formatted file with motif incidence per read

Positional arguments:
    <sequencefile>                name of input BAM/SAM/FASTA/FASTQ file

Options:
    --fmt sam|fastx               format of input file [default: fastx]
    -m, --min-k [integer]         smallest target repeat length [default: 4]
    -M, --max-k [integer]         largest target repeat length [default: 16]
    -r, --min-repeats [integer]   minimum number of consecutive repeats [default: 2]
    --motifs [string]             list of target motifs, separated with '|'

Note: if --motifs is specified, options -m and -M have no effect
"""

__docopt_converters__ = [
    lambda min_k: int(min_k),
    lambda max_k: int(max_k),
    lambda min_repeats: int(min_repeats),
]

__docopt_tests__ = {
    lambda min_k, max_k: 0 < min_k < max_k: "not satisfied: 0 < m < M",
    lambda min_repeats: min_repeats > 0: "--min-repeats must be integer > 0",
    lambda fmt: fmt in {"sam", "fastx"}: "unknown value of --fmt",
}


def interpret_arguments(fmt, min_k, max_k, min_repeats, motifs):
    """Parse and check arguments"""
    if fmt == "fastx":
        manager = FastxFile
    else:
        raise NotImplementedError("--fmt={}".format(fmt))
    if motifs:
        print("WARNING: Using --motifs, ignoring -m and -M", file=stderr)
        motif_patterns = OrderedDict([
            [motif, get_circular_pattern(motif)]
            for motif in motifs.split("|")
        ])
    else:
        raise NotImplementedError("Full motif sweep")
    return manager, motif_patterns


def main(sequencefile, fmt, min_k, max_k, min_repeats, motifs, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    manager, motif_patterns = interpret_arguments(
        fmt, min_k, max_k, min_repeats, motifs
    )
    print(*motif_patterns.keys(), sep="\t", file=file)
    with FastxFile(sequencefile) as fastx:
        for entry in fastx:
            motif_counts = array([
                sum(1 for _ in p.finditer(entry.sequence, overlapped=True))
                for p in motif_patterns.values()
            ])
            if (motif_counts != 0).any():
                print(*motif_counts, sep="\t", file=file)

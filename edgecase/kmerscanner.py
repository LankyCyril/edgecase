#!/usr/bin/env python3
from numpy import zeros, array, cumsum
from multiprocessing import Pool
from edgecase.io import ReadFileChain
from pysam import FastxFile
from functools import partial
from tqdm import tqdm
from argparse import ArgumentParser
from sys import stderr


def get_edge_density(read, pattern, head, tail):
    """Calculate density of pattern in head or tail of read"""
    if len(read.sequence) == 0:
        return 0
    if head:
        subsequence = read.sequence[:head]
    elif tail:
        subsequence = read.sequence[-tail:]
    pattern_matches = pattern.findall(subsequence, overlapped=True)
    return len(pattern_matches) / len(subsequence)


def calculate_density(read, pattern, cutoff, window_size, head, tail):
    """Calculate density of pattern hits in a rolling window along given read"""
    if cutoff: # if cutoff specified, filter by hard cutoff
        edge_density = get_edge_density(
            read, pattern, head, tail
        )
        passes_filter = (edge_density > cutoff)
    else: # otherwise, allow all
        passes_filter = True
    if not passes_filter: # nothing to search for
        density_array = zeros(1)
    else: # create array of hits; axis 0 is positions, axis 1 is patterns
        read_length = len(read.sequence)
        canvas = zeros(read_length, dtype=bool)
        pattern_positions = array([
            match.start() for match
            in pattern.finditer(read.sequence, overlapped=True)
        ])
        if len(pattern_positions):
            canvas[pattern_positions] = True
        if read_length <= window_size: # use one window:
            density_array = (canvas.sum(axis=0) / read_length)
        else: # use rolling window:
            roller = cumsum(canvas, axis=0)
            roller[window_size:] = roller[window_size:] - roller[:-window_size]
            density_array = roller[window_size-1:] / window_size
    return read.name, density_array, passes_filter


def pattern_scanner(read_iterator, pattern, cutoff, window_size, head, tail, num_reads, jobs):
    """Calculate density of pattern hits in a rolling window along each read"""
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density, pattern=pattern,
            window_size=window_size, head=head, tail=tail, cutoff=cutoff
        )
        # lazy multiprocess evaluation:
        read_density_iterator = pool.imap_unordered(
            density_calculator, read_iterator
        )
        # iterate pairs (read.name, density_array), same as calculate_density():
        yield from tqdm(
            read_density_iterator, desc="Calculating kmer density",
            unit="read", total=num_reads
        )


USAGE = "python3 {} [options] fastq > txt".format(__file__)

ARG_RULES = {
    ("fastqs",): {
        "help": "name of input FASTQ file",
        "nargs": "+"
    },
    ("--kmer",): {
        "help": "target kmer sequence (default TTAGGG)",
        "default": "TTAGGG", "metavar": "M"
    },
    ("--head-test",): {
        "help": "length of head to use for density filter (if specified)",
        "default": None, "type": int, "metavar": "H"
    },
    ("--tail-test",): {
        "help": "length of tail to use for density filter (if specified)",
        "default": None, "type": int, "metavar": "T"
    },
    ("-c", "--cutoff"): {
        "help": "use hard cutoff for density (default None)",
        "default": None, "type": float, "metavar": "C"
    },
    ("-w", "--window-size"): {
        "help": "size of the rolling window (default 120)",
        "default": 120, "type": int, "metavar": "W"
    },
    ("-n", "--num-reads"): {
        "help": "expected number of reads in input (for progress display)",
        "default": None, "type": int, "metavar": "N"
    },
    ("-j", "--jobs"): {
        "help": "number of jobs to run in parallel (default 1)",
        "default": 1, "type": int, "metavar": "J"
    }
}


if __name__ == "__main__":
    # parse and check arguments:
    parser = ArgumentParser(usage=USAGE)
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    if (args.head_test is not None) and (args.tail_test is not None):
        raise ValueError("Can only specify one of --head, --tail")
    elif (args.cutoff is not None) and (args.head_test is None) and (args.tail_test is None):
        raise ValueError("--cutoff has no effect without a head/tail test")
    elif ((args.head_test is not None) or (args.tail_test is not None)) and (args.cutoff is None):
        print("Warning: --head/--tail has no effect without --cutoff", file=stderr)
    # dispatch data to subroutines:
    pattern = "FIXME"
    # scan fastq for target kmer query, parallelizing on reads:
    with ReadFileChain(args.fastqs, FastxFile) as read_iterator:
        scanner = pattern_scanner(
            read_iterator, pattern,
            window_size=args.window_size, head=args.head, tail=args.tail,
            cutoff=args.cutoff, num_reads=args.num_reads, jobs=args.jobs
        )
        # output densities of reads that pass filter:
        for read_name, density_array, passes_filter in scanner:
            if passes_filter:
                print(read_name, *density_array, sep="\t")

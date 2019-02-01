from numpy import zeros, array, cumsum
from multiprocessing import Pool
from edgecaselib.io import ReadFileChain
from pysam import FastxFile
from functools import partial
from tqdm import tqdm
from sys import stderr


def get_edge_density(read, pattern, head_test, tail_test):
    """Calculate density of pattern in head_test or tail_test of read"""
    if len(read.sequence) == 0:
        return 0
    if head_test:
        subsequence = read.sequence[:head_test]
    elif tail_test:
        subsequence = read.sequence[-tail_test:]
    pattern_matches = pattern.findall(subsequence, overlapped=True)
    return len(pattern_matches) / len(subsequence)


def calculate_density(read, pattern, cutoff, window_size, head_test, tail_test):
    """Calculate density of pattern hits in a rolling window along given read"""
    if cutoff: # if cutoff specified, filter by hard cutoff
        edge_density = get_edge_density(
            read, pattern, head_test, tail_test
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


def pattern_scanner(read_iterator, pattern, cutoff, window_size, head_test, tail_test, num_reads, jobs):
    """Calculate density of pattern hits in a rolling window along each read"""
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density, pattern=pattern,
            window_size=window_size, head_test=head_test, tail_test=tail_test, cutoff=cutoff
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


def main(args):
    # parse and check arguments:
    if (args.head_test is not None) and (args.tail_test is not None):
        raise ValueError("Can only specify one of --head-test, --tail-test")
    elif (args.cutoff is not None) and (args.head_test is None) and (args.tail_test is None):
        raise ValueError("--cutoff has no effect without a head/tail test")
    elif ((args.head_test is not None) or (args.tail_test is not None)) and (args.cutoff is None):
        print("Warning: head/tail test has no effect without --cutoff", file=stderr)
    # dispatch data to subroutines:
    pattern = "FIXME"
    # scan fastq for target kmer query, parallelizing on reads:
    with ReadFileChain(args.fastqs, FastxFile) as read_iterator:
        scanner = pattern_scanner(
            read_iterator, pattern, window_size=args.window_size,
            head_test=args.head_test, tail_test=args.tail_test,
            cutoff=args.cutoff, num_reads=args.num_reads, jobs=args.jobs
        )
        # output densities of reads that pass filter:
        for read_name, density_array, passes_filter in scanner:
            if passes_filter:
                print(read_name, *density_array, sep="\t")

from sys import stdout, stderr
from regex import compile, IGNORECASE
from numpy import zeros, array, cumsum
from multiprocessing import Pool
from edgecaselib.io import ReadFileChain
from pysam import AlignmentFile
from types import SimpleNamespace
from functools import partial
from tqdm import tqdm


def get_circular_pattern(kmer):
    """Convert kmer into circular regex pattern (e.g., r'TCGA|CGAT|GATC|ATCG' for TCGA)"""
    inversions = {kmer[i:]+kmer[:i] for i in range(len(kmer))}
    return compile(r'|'.join(inversions), flags=IGNORECASE)


def get_edge_density(entry, pattern, head_test, tail_test):
    """Calculate density of pattern in head_test or tail_test of read"""
    if len(entry.query_sequence) == 0:
        return 0
    if head_test:
        subsequence = entry.query_sequence[:head_test]
    elif tail_test:
        subsequence = entry.query_sequence[-tail_test:]
    pattern_matches = pattern.findall(subsequence, overlapped=True)
    return len(pattern_matches) / len(subsequence)


def calculate_density(entry, pattern, cutoff, window_size, head_test, tail_test):
    """Calculate density of pattern hits in a rolling window along given read"""
    if cutoff: # if cutoff specified, filter by hard cutoff
        edge_density = get_edge_density(
            entry, pattern, head_test, tail_test
        )
        passes_filter = (edge_density > cutoff)
    else: # otherwise, allow all
        passes_filter = True
    if passes_filter: # calculations will make sense
        read_length = len(entry.query_sequence)
        canvas = zeros(read_length, dtype=bool)
        pattern_positions = array([
            match.start() for match
            in pattern.finditer(entry.query_sequence, overlapped=True)
        ])
        if len(pattern_positions):
            canvas[pattern_positions] = True
        if read_length <= window_size: # use one window:
            density_array = (canvas.sum(axis=0) / read_length)
        else: # use rolling window:
            roller = cumsum(canvas, axis=0)
            roller[window_size:] = roller[window_size:] - roller[:-window_size]
            density_array = roller[window_size-1:] / window_size
        return entry, density_array
    else: # effectively skip
        return None, zeros(1)

def pattern_scanner(entry_iterator, pattern, cutoff, window_size, head_test, tail_test, num_reads, jobs):
    """Calculate density of pattern hits in a rolling window along each read"""
    simple_entry_iterator = (
        SimpleNamespace(
            query_name=entry.query_name,
            flag=entry.flag,
            reference_name=entry.reference_name,
            reference_start=entry.reference_start,
            mapping_quality=entry.mapping_quality,
            query_sequence=entry.query_sequence
        )
        for entry in entry_iterator
    )
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density, pattern=pattern,
            window_size=window_size, head_test=head_test, tail_test=tail_test, cutoff=cutoff
        )
        # lazy multiprocess evaluation:
        read_density_iterator = pool.imap_unordered(
            density_calculator, simple_entry_iterator
        )
        # iterate pairs (entry.query_name, density_array), same as calculate_density():
        yield from tqdm(
            read_density_iterator, desc="Calculating kmer density",
            unit="read", total=num_reads
        )


def main(args, file=stdout):
    # parse and check arguments:
    if (args.head_test is not None) and (args.tail_test is not None):
        raise ValueError("Can only specify one of --head-test, --tail-test")
    elif (args.cutoff is not None) and (args.head_test is None) and (args.tail_test is None):
        raise ValueError("--cutoff has no effect without a head/tail test")
    elif ((args.head_test is not None) or (args.tail_test is not None)) and (args.cutoff is None):
        print("Warning: head/tail test has no effect without --cutoff", file=stderr)
    # dispatch data to subroutines:
    pattern = get_circular_pattern(args.kmer)
    # scan fastq for target kmer query, parallelizing on reads:
    with ReadFileChain(args.bams, AlignmentFile) as entry_iterator:
        scanner = pattern_scanner(
            entry_iterator, pattern, window_size=args.window_size,
            head_test=args.head_test, tail_test=args.tail_test,
            cutoff=args.cutoff, num_reads=args.num_reads, jobs=args.jobs
        )
        # output densities of reads that pass filter:
        for entry, density_array in scanner:
            if entry: # non-null result, entry passed filters
                meta_fields = [
                    entry.query_name, entry.flag, entry.reference_name,
                    entry.reference_start, entry.mapping_quality,
                    args.kmer
                ]
                print(*meta_fields, sep="\t", end="\t", file=file)
                print(*density_array, sep=",", file=file)

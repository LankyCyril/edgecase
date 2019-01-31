#!/usr/bin/env python3
from numpy import zeros, array, cumsum
from multiprocessing import Pool
from edgecase.io import ReadFileChain
from pysam import FastxFile
from itertools import islice
from functools import partial
from tqdm import tqdm


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


def fastq_scanner(fastqs, pattern, cutoff, window_size, head, tail, num_reads, jobs):
    """Thin wrapper over pattern_scanner() providing read_iterator from fastq file"""
    with ReadFileChain(fastqs, FastxFile) as read_iterator:
        # only take the first num_reads entries (`None` takes all):
        capped_read_iterator = islice(read_iterator, num_reads)
        # parallelize scan:
        yield from pattern_scanner(
            capped_read_iterator, pattern, jobs=jobs,
            window_size=window_size, head=head, tail=tail,
            cutoff=cutoff, num_reads=num_reads
        )

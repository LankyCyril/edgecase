from sys import stdout, stderr
from regex import compile, IGNORECASE
from numpy import zeros, array, cumsum
from multiprocessing import Pool
from edgecaselib.util import ReadFileChain
from pysam import AlignmentFile
from types import SimpleNamespace
from functools import partial
from tqdm import tqdm


def get_circular_pattern(motif):
    """Convert motif into circular regex pattern (e.g., r'TCGA|CGAT|GATC|ATCG' for TCGA)"""
    atom_pattern = compile(r'[ACGT.]|\[[ACGT]+\]', flags=IGNORECASE)
    atoms = atom_pattern.findall(motif)
    if "".join(atoms) != motif:
        raise ValueError("Could not parse motif: {}".format(motif))
    inversions = {
        "".join(atoms[i:] + atoms[:i])
        for i in range(len(atoms))
    }
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
        desc = "Calculating density of " + pattern.pattern.split(r'|')[0]
        yield from tqdm(
            read_density_iterator, desc=desc,
            unit="read", total=num_reads
        )


def main(bams, motif="TTAGGG", head_test=None, tail_test=None, cutoff=None, window_size=120, num_reads=None, jobs=1, file=stdout, **kwargs):
    # parse and check arguments:
    if (head_test is not None) and (tail_test is not None):
        raise ValueError("Can only specify one of --head-test, --tail-test")
    elif (cutoff is not None) and (head_test is None) and (tail_test is None):
        raise ValueError("--cutoff has no effect without a head/tail test")
    elif ((head_test is not None) or (tail_test is not None)) and (cutoff is None):
        print("Warning: head/tail test has no effect without --cutoff", file=stderr)
    # dispatch data to subroutines:
    pattern = get_circular_pattern(motif)
    # scan fastq for target motif query, parallelizing on reads:
    with ReadFileChain(bams, AlignmentFile) as entry_iterator:
        scanner = pattern_scanner(
            entry_iterator, pattern, window_size=window_size,
            head_test=head_test, tail_test=tail_test,
            cutoff=cutoff, num_reads=num_reads, jobs=jobs
        )
        # output densities of reads that pass filter:
        for entry, density_array in scanner:
            if entry: # non-null result, entry passed filters
                meta_fields = [
                    entry.query_name, entry.flag, entry.reference_name,
                    entry.reference_start, entry.mapping_quality,
                    motif
                ]
                print(*meta_fields, sep="\t", end="\t", file=file)
                print(*density_array, sep=",", file=file)

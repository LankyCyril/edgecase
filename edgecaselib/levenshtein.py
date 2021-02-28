from sys import stdout
from numba import njit
from pysam import AlignmentFile
from edgecaselib.formats import filter_bam
from numpy import zeros, array, uint32, uint8
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations_with_replacement
from collections import defaultdict
from edgecaselib.util import progressbar


__warning__ = """The `levenshtein` subprogram is in development!
Pairwise distance computation is O(n^2) and is not suited for large scale
experiments."""

__doc__ = """edgeCase levenshtein: pairwise edit distance among telomeric reads

Usage: {0} levenshtein [-f flagspec]... [-F flagspec]... [-q integer]
       {1}             [-j integer] <sequencefile>

Output:
    TSV-formatted file with computed relative pairwise distances (per-arm)

Positional arguments:
    <sequencefile>                     name of input BAM/SAM file

Options:
    -j, --jobs [integer]               number of jobs to run in parallel [default: 1]

Input filtering options:
    -f, --flags [flagspec]             process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]       process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]        process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda min_quality:
        None if (min_quality is None) else int(min_quality),
    lambda jobs: int(jobs),
]


def load_bam_as_dict(alignment, samfilters):
    """Load BAM entries as dictionary: chrom -> qname -> (mappos, int8array)"""
    dna2int = lambda seq: ((array(list(seq.upper())).view(uint32) - 2) >> 1) & 3
    bam_dict = defaultdict(dict)
    for entry in filter_bam(alignment, samfilters, desc="Reading BAM"):
        bam_dict[entry.reference_name][entry.qname] = (
            entry.reference_start,
            dna2int(entry.seq[entry.query_alignment_start:]).astype(uint8),
        )
    return bam_dict


@njit(nogil=True)
def ld(v, w):
    """Calculate levenshtein distance between two int8arrays"""
    m, n = len(v), len(w)
    dp = zeros((m+1, n+1), dtype=uint32)
    for i in range(m+1):
        for j in range(n+1):
            if i == 0:
                dp[i, j] = j
            elif j == 0:
                dp[i, j] = i
            elif v[i-1] == w[j-1]:
                dp[i, j] = dp[i-1, j-1]
            else:
                dp[i, j] = 1 + min(dp[i, j-1], dp[i-1, j], dp[i-1, j-1])
    return dp[m, n]


def get_relative_read_ld(aname, bname, adata, bdata):
    """Calculate relative levenshtein distance between overlapping parts of two reads"""
    sra, A = adata
    srb, B = bdata
    if sra < srb:
        _A, _B = A[srb-sra:], B
    elif sra > srb:
        _A, _B = A, B[sra-srb:]
    else:
        _A, _B = A, B
    overlap_length = min(len(_A), len(_B))
    _A, _B = _A[:overlap_length], _B[:overlap_length]
    if (_A == _B).all():
        return aname, bname, 0
    elif overlap_length > 0:
        return aname, bname, ld(_A, _B) / overlap_length
    else:
        return aname, bname, 1


def calculate_chromosome_lds(chrom, entries, jobs):
    """Calculate pairwise relative levenshtein distances between all reads mapping to one chromosome"""
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        workers = [
            pool.submit(
                get_relative_read_ld,
                aname, bname, entries[aname], entries[bname],
            )
            for aname, bname in combinations_with_replacement(
                sorted(entries.keys()), r=2,
            )
        ]
        iterator = progressbar(
            as_completed(workers), desc=chrom, unit="pair", total=len(workers),
        )
        for worker in iterator:
            yield worker.result()


def main(sequencefile, flags, flag_filter, min_quality, jobs=1, file=stdout, **kwargs):
    print("#rname", "qname1", "qname2", "relative_ld", sep="\t", file=file)
    with AlignmentFile(sequencefile) as alignment:
        bam_dict = load_bam_as_dict(
            alignment, samfilters=[flags, flag_filter, min_quality],
        )
    for chrom, entries in bam_dict.items():
        ld_iterator = calculate_chromosome_lds(chrom, entries, jobs)
        for qname1, qname2, relative_ld in ld_iterator:
            print(chrom, qname1, qname2, relative_ld, sep="\t", file=file)

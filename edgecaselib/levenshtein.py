from sys import stdout
from pysam import AlignmentFile
from edgecaselib.formats import filter_bam
from numpy import zeros, array, uint32, uint8, ndarray
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations_with_replacement
from collections import defaultdict
from edgecaselib.util import progressbar
from pandas import DataFrame
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage


__warning__ = """Pairwise distance computation is O(n^2) and is not suited for
large scale experiments."""

__doc__ = """edgeCase levenshtein: pairwise edit distance among telomeric reads

Usage: {0} levenshtein [-f flagspec]... [-F flagspec]... [-q integer]
       {1}             [-j integer] [-c] <sequencefile>

Output:
    TSV-formatted file with computed relative pairwise distances (per-arm).
    If `-c` (`--cluster`) is passed, converts the distances to a square matrix
    and clusters it with the Ward method over the euclidean metric.

Positional arguments:
    <sequencefile>                     name of input BAM/SAM file

Options:
    -j, --jobs [integer]               number of jobs to run in parallel [default: 1]
    -c, --cluster                      perform clustering after LD calculation

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


try:
    from edlib import align
except (ModuleNotFoundError, ImportError):
    from numba import njit
    dna2int = lambda seq: (
        ((array(list(seq.upper())).view(uint32) - 2) >> 1) & 3
    ).astype(uint8)
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
else:
    dna2int = lambda _:_
    def ld(v, w):
        """Calculate levenshtein distance between two nucleotide sequences"""
        return align(v, w)["editDistance"]


def load_bam_as_dict(alignment, samfilters):
    """Load BAM entries as dictionary: chrom -> qname -> (mappos, int8array)"""
    bam_dict = defaultdict(dict)
    for entry in filter_bam(alignment, samfilters, desc="Reading BAM"):
        bam_dict[entry.reference_name][entry.qname] = (
            entry.reference_start,
            dna2int(entry.seq[entry.query_alignment_start:]),
        )
    return bam_dict


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
    if isinstance(_A, str) and (_A == _B):
        return aname, bname, 0
    elif isinstance(_A, ndarray) and (_A == _B).all():
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


def clustered_rld(ld_iterator, columns, metric="euclidean", method="ward"):
    narrowform = DataFrame(data=ld_iterator, columns=columns)
    pivot_kws = dict(index=columns[0], columns=columns[1], values=columns[2])
    triu_fillna = narrowform.pivot(**pivot_kws).fillna(0)
    squarish = triu_fillna.T + triu_fillna
    Z = linkage(
        squareform(squarish), metric=metric, method=method,
        optimal_ordering=False,
    )
    leaves = dendrogram(Z, no_plot=True)["leaves"]
    data2d = squarish.iloc[leaves, leaves]
    return DataFrame(Z), data2d


def main(sequencefile, flags, flag_filter, min_quality, cluster=False, jobs=1, file=stdout, **kwargs):
    HEADER = ["rname", "qname1", "qname2", "relative_ld"]
    if not cluster:
        print("#", end="", file=file)
        print(*HEADER, sep="\t", file=file)
    with AlignmentFile(sequencefile) as alignment:
        bam_dict = load_bam_as_dict(
            alignment, samfilters=[flags, flag_filter, min_quality],
        )
    for chrom, entries in bam_dict.items():
        ld_iterator = calculate_chromosome_lds(chrom, entries, jobs)
        if cluster:
            Z, data2d = clustered_rld(ld_iterator, columns=HEADER[1:])
            print("##linkage", chrom, sep="\t", file=file)
            print(Z.to_csv(sep="\t", index=False, header=False))
            print("##data2d", chrom, sep="\t", file=file)
            print(data2d.to_csv(sep="\t"))
        else:
            for data in ld_iterator:
                print(chrom, *data, sep="\t", file=file)

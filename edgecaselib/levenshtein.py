from sys import stdout
from numba import njit
from pysam import AlignmentFile
from edgecaselib.formats import filter_bam
from numpy import zeros, array, uint32, uint8, concatenate, nan, isnan
from collections import defaultdict
from tqdm import tqdm
from pandas import DataFrame


def load_bam_as_dict(alignment, samfilters):
    dna2int = lambda seq: ((array(list(seq.upper())).view(uint32) - 2) >> 1) & 3
    bam_dict = defaultdict(dict)
    for entry in filter_bam(alignment, samfilters, desc="Reading BAM"):
        bam_dict[entry.reference_name][entry.qname] = (
            entry.reference_start,
            dna2int(entry.seq[entry.query_alignment_start:]).astype(uint8)
        )
    return bam_dict


@njit(nogil=True)
def ld(v, w):
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


def get_relative_read_ld(sra, A, srb, B, return_bases=False):
    if sra < srb:
        _A, _B = A[srb-sra:], B
    elif sra > srb:
        _A, _B = A, B[sra-srb:]
    else:
        _A, _B = A, B
    overlap_length = min(len(_A), len(_B))
    _A, _B = _A[:overlap_length], _B[:overlap_length]
    if overlap_length > 0:
        distance = ld(_A, _B) / overlap_length
        if return_bases:
            bases = concatenate([_A, _B])
        else:
            bases = zeros(shape=0)
    else:
        distance = 1
        bases = zeros(shape=0)
    return distance, overlap_length, bases


def calculate_chromosome_lds(chrom, entries):
    lds = DataFrame(
        data=nan, columns=sorted(entries.keys()), index=sorted(entries.keys())
    )
    for aname, (sra, A) in tqdm(entries.items(), desc=chrom, unit="read"):
        for bname, (srb, B) in entries.items():
            if aname == bname:
                distance = 0
            elif not isnan(lds.loc[aname, bname]):
                distance = lds.loc[bname, aname]
            else:
                distance, *_ = get_relative_read_ld(sra, A, srb, B)
            lds.loc[aname, bname] = distance
    return lds.fillna(1)


def main(bam, flags, flags_any, flag_filter, min_quality, jobs=1, file=stdout, **kwargs):
    samfilters = [flags, flags_any, flag_filter, min_quality]
    with AlignmentFile(bam) as alignment:
        bam_dict = load_bam_as_dict(alignment, samfilters)
    for chrom, entries in bam_dict.items():
        calculate_chromosome_lds(chrom, entries)

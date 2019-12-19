#!/usr/bin/env python
from sys import argv
from numba import njit
from numpy import array, uint32, zeros, nan, isnan, save
from pysam import AlignmentFile
from collections import defaultdict
from os import path
from tqdm import tqdm
from pandas import DataFrame
from matplotlib.pyplot import switch_backend
from seaborn import clustermap


def load_bam(filename, f, F):
    dna2int = lambda seq: ((array(list(seq.upper())).view(uint32) - 2) >> 1) & 3
    bam_data = defaultdict(dict)
    with AlignmentFile(filename) as bam:
        for entry in bam:
            if (entry.flag & f == f) and (entry.flag & F == 0) and (entry.seq):
                bam_data[entry.reference_name][entry.qname] = (
                    entry.reference_start,
                    dna2int(entry.seq[entry.query_alignment_start:])
                )
    return bam_data


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


def read_ld(sra, A, srb, B):
    if sra < srb:
        _A, _B = A[srb-sra:], B
    elif sra > srb:
        _A, _B = A, B[sra-srb:]
    else:
        _A, _B = A, B
    minlen = min(len(_A), len(_B))
    _A, _B = _A[:minlen], _B[:minlen]
    if minlen > 0:
        return ld(_A, _B) / minlen
    else:
        return 1


def compute_pairwise_lds(entries, chrom="[unspecified chromosome]"):
    pairwise_lds = DataFrame(
        data=nan, columns=sorted(entries.keys()), index=sorted(entries.keys())
    )
    for aname, (sra, A) in tqdm(entries.items(), desc=chrom):
        for bname, (srb, B) in entries.items():
            if aname == bname:
                pairwise_lds.loc[aname, bname] = 0
            elif isnan(pairwise_lds.loc[aname, bname]):
                if not isnan(pairwise_lds.loc[bname, aname]):
                    distance = pairwise_lds.loc[bname, aname]
                else:
                    distance = read_ld(sra, A, srb, B)
                pairwise_lds.loc[aname, bname] = distance
    return pairwise_lds


def main(filename, f=0, F=0, out_dir="."):
    switch_backend("Agg")
    bam_data = load_bam(filename, int(f), int(F))
    for chrom, entries in bam_data.items():
        if len(entries) >= 2:
            pairwise_lds = compute_pairwise_lds(entries, chrom)
            cm = clustermap(data=pairwise_lds, cmap="viridis_r", method="ward")
            cm.ax_heatmap.set(xticks=[], yticks=[])
            pdf = path.join(out_dir, chrom+".pdf")
            cm.fig.savefig(pdf, bbox_inches="tight")
            npy = path.join(out_dir, chrom+".npy")
            save(npy, cm.dendrogram_col.linkage)
            tsv = path.join(out_dir, chrom+".tsv")
            cm.data2d.to_csv(tsv, sep="\t")
    return 0


if __name__ == "__main__":
    returncode = main(*argv[1:])
    exit(returncode)

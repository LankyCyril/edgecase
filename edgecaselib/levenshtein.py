from sys import stdout, stderr
from numba import njit
from pysam import AlignmentFile
from edgecaselib.formats import filter_bam
from numpy import zeros, array, uint32, uint8, concatenate, nan, isnan, unique
from collections import defaultdict
from tqdm import tqdm
from pandas import DataFrame
from matplotlib.pyplot import switch_backend
from seaborn import clustermap
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_score
from scipy.stats import mannwhitneyu


CLUSTERMAP_FIGSIZE = (10, 10)


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


def generate_clustermap(lds, metric="euclidean", method="ward", cmap="viridis_r", vmax=.25):
    cm = clustermap(
        data=lds, metric=metric, method=method,
        cmap=cmap, vmin=0, vmax=.25, figsize=CLUSTERMAP_FIGSIZE
    )
    cm.cax.set(visible=False)
    cm.ax_heatmap.set(xticks=[], yticks=[])
    return cm


def get_mwu_pvals(lds, c1_names, c2_names):
    c1_ingroup = lds.loc[c1_names, c1_names].values.flatten()
    c2_ingroup = lds.loc[c2_names, c2_names].values.flatten()
    outgroup = lds.loc[c1_names, c2_names].values.flatten()
    c1_pval = mannwhitneyu(c1_ingroup, outgroup, alternative="less")[1]
    c2_pval = mannwhitneyu(c2_ingroup, outgroup, alternative="less")[1]
    return c1_pval, c2_pval


def get_major_clusters(lds, linkage):
    """Find two major clusters; so far, only works if there is not more than one outlier read"""
    labels = fcluster(linkage, 2, criterion="maxclust")
    uniq, counts = unique(labels, return_counts=True)
    if 1 in counts:
        labels = fcluster(linkage, 3, criterion="maxclust")
        uniq, counts = unique(labels, return_counts=True)
    good_labels = list(set(uniq[counts!=1]))
    if len(good_labels) != 2:
        return [], [], nan, nan, nan
    bad_labels = uniq[counts==1]
    if len(bad_labels) == 0:
        filtered_labels = labels
    elif len(bad_labels) == 1:
        filtered_labels = labels[labels!=bad_labels[0]]
    else:
        return [], [], nan, nan, nan
    c1c2_names, c1_names, c2_names = (
        lds.index[(labels==good_labels[0]) | (labels==good_labels[1])],
        lds.index[labels==good_labels[0]], lds.index[labels==good_labels[1]]
    )
    filtered_lds = lds.loc[c1c2_names, c1c2_names]
    try:
        c1_pval, c2_pval = get_mwu_pvals(lds, c1_names, c2_names)
        silh_score = silhouette_score(filtered_lds, filtered_labels)
    except ValueError:
        return [], [], nan, nan, nan
    else:
        c1_names, c2_names = (
            list(lds.index[labels==good_labels[0]]),
            list(lds.index[labels==good_labels[1]])
        )
        return c1_names, c2_names, c1_pval, c2_pval, silh_score


def main(bam, output_dir, flags, flags_any, flag_filter, min_quality, jobs=1, file=stdout, **kwargs):
    switch_backend("Agg")
    #samfilters = [flags, flags_any, flag_filter, min_quality]
    #with AlignmentFile(bam) as alignment:
    #    bam_dict = load_bam_as_dict(alignment, samfilters)
    #for chrom, entries in bam_dict.items():
    #    lds = calculate_chromosome_lds(chrom, entries)
    from pickle import load
    with open("chr22.pkl", mode="rb") as pkl:
        lds = load(pkl)
    if True:
        cm = generate_clustermap(lds)
        c1_names, c2_names, c1_pval, c2_pval, silh_score = get_major_clusters(
            lds, cm.dendrogram_row.linkage
        )
        if not isnan(silh_score):
            print(c1_names, c2_names, c1_pval, c2_pval, silh_score)
        else:
            msg = "Warning: not implemented: complex hierarchy in clustering"
            print(msg, file=stderr)
        #pval_text = "N/A" if isnan(pval) else "{}".format(pval)
        #silh_text = "N/A" if isnan(silh_score) else "{:.3f}".format(silh_score)
        #break

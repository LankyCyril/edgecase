from sys import stdout, stderr
from warnings import filterwarnings, resetwarnings
from numba import njit
from pysam import AlignmentFile
from edgecaselib.formats import filter_bam
from numpy import zeros, array, uint32, uint8, concatenate, nan, isnan, unique
from numpy import linspace, vstack
from collections import defaultdict
from edgecaselib.util import progressbar
from pandas import DataFrame, read_csv, concat
from matplotlib.pyplot import switch_backend
from seaborn import clustermap
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_score
from scipy.stats import mannwhitneyu
from os import path
from edgecaselib.densityplot import shorten_chrom_name
from statsmodels.stats.multitest import multipletests


CLUSTERMAP_FIGSIZE = (10, 10)
CLUSTERMAP_CMAP = "viridis_r"
CLUSTERMAP_VMAX = .15


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
    read_iterator = progressbar(entries.items(), desc=chrom, unit="read")
    for aname, (sra, A) in read_iterator:
        for bname, (srb, B) in entries.items():
            if aname == bname:
                distance = 0
            elif not isnan(lds.loc[aname, bname]):
                distance = lds.loc[bname, aname]
            else:
                distance, *_ = get_relative_read_ld(sra, A, srb, B)
            lds.loc[aname, bname] = distance
    return lds.fillna(1)


def generate_clustermap(lds, metric="euclidean", method="ward", cmap=CLUSTERMAP_CMAP, vmax=CLUSTERMAP_VMAX):
    try:
        cm = clustermap(
            data=lds, metric=metric, method=method,
            cmap=cmap, vmin=0, vmax=vmax, figsize=CLUSTERMAP_FIGSIZE
        )
    except (ValueError, ModuleNotFoundError):
        return None
    else:
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


def get_clusters(lds, linkage):
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


def format_pval(pval):
    if isnan(pval):
        return "N/A"
    elif pval >= .00001:
        return "{:.5f}".format(pval)
    else:
        return "{:.1e}".format(pval)


def display_pvals(c1_pval, c2_pval):
    return "$p_1$={}, $p_2$={}".format(
        format_pval(c1_pval), format_pval(c2_pval)
    )


def generate_kmerscanner_file(kmerscanner_file, c1_names, c2_names, output_dir, chrom):
    kmerscanner_dat = read_csv(kmerscanner_file, sep="\t")
    kmerscanner_dat_hap1 = kmerscanner_dat[
        kmerscanner_dat["#name"].isin(c1_names)
    ].copy()
    kmerscanner_dat_hap1["chrom"] = kmerscanner_dat_hap1["chrom"].apply(
        lambda s: s + ":haplotype 1"
    )
    kmerscanner_dat_hap2 = kmerscanner_dat[
        kmerscanner_dat["#name"].isin(c2_names)
    ].copy()
    kmerscanner_dat_hap2["chrom"] = kmerscanner_dat_hap2["chrom"].apply(
        lambda s: s + ":haplotype 2"
    )
    haplo_dat = concat(
        [kmerscanner_dat_hap1, kmerscanner_dat_hap2], axis=0
    )
    haplo_dat.to_csv(
        path.join(output_dir, chrom+".dat.gz"), compression="gzip",
        sep="\t", index=False
    )


def generate_report(report_rows, adj="bonferroni"):
    report = DataFrame(
        data=report_rows,
        columns=[
            "#chrom", "cluster1_size", "cluster2_size", "silhouette_score",
            "cluster1_pvalue", "cluster2_pvalue"
        ]
    )
    pvals = report[["cluster1_pvalue", "cluster2_pvalue"]].values.flatten()
    p_adjusted = multipletests(pvals, method=adj)[1].reshape(-1, 2)
    report["cluster1_p_adjusted"] = p_adjusted[:,0]
    report["cluster2_p_adjusted"] = p_adjusted[:,1]
    return report


def generate_pdf(cm, c1_pval, c2_pval, silh_score, output_dir, chrom, cmap=CLUSTERMAP_CMAP, vmax=CLUSTERMAP_VMAX):
    cm.ax_col_dendrogram.clear()
    cm.ax_col_dendrogram.imshow(
        vstack([linspace(0, 1, 256)]*2), aspect="auto", cmap=cmap
    )
    gp = cm.ax_col_dendrogram.get_position()
    cm.ax_col_dendrogram.set(
        position=[
            gp.x0*.25+gp.x1*.75, gp.y0*.55+gp.y1*.45, gp.x1*.1-gp.x0*.1, .015
        ],
        zorder=float("inf")
    )
    cm.ax_col_dendrogram.text(
        x=-672, y=.8, va="center", ha="left", fontsize=19, s="Distance:  0"
    )
    cm.ax_col_dendrogram.text(
        x=272, y=.8, va="center", ha="left", fontsize=19, s="{}+".format(vmax)
    )
    cm.ax_col_dendrogram.set_axis_off()
    silh_text = "N/A" if isnan(silh_score) else "{:.3f}".format(silh_score)
    cm.ax_col_dendrogram.text(
        x=-672, y=3, va="top", ha="left", fontsize=19,
        s="Silhouette score: "+silh_text
    )
    cm.ax_col_dendrogram.text(
        x=-672, y=6.6, va="top", ha="left", fontsize=19,
        s=display_pvals(c1_pval, c2_pval)
    )
    cm.ax_heatmap.text(
        x=0, y=0, va="center", ha="left", fontsize=20,
        s=shorten_chrom_name(chrom)+"\n\n\n"
    )
    filename = path.join(output_dir, chrom+".pdf")
    cm.fig.savefig(filename, bbox_inches="tight")


def batch_generate_pdfs(cms, report, output_dir):
    if output_dir:
        for _, report_entry in report.iterrows():
            chrom = report_entry["#chrom"]
            if chrom in cms:
                generate_pdf(
                    cms[chrom], report_entry["cluster1_p_adjusted"],
                    report_entry["cluster2_p_adjusted"],
                    report_entry["silhouette_score"],
                    output_dir, chrom
                )


def hide_stats_warnings(state=True):
    if state:
        filterwarnings("ignore", message="invalid value encountered")
        filterwarnings(
            "ignore",
            message="looks suspiciously like an uncondensed distance matrix"
        )
    else:
        resetwarnings()


def warn_about_unsupported_hierarchy(chrom):
    msg = "({}): not implemented: complex clustering hierarchy or too few reads"
    print("Warning", msg.format(chrom), file=stderr)


def main(bam, kmerscanner_file, output_dir, flags, flags_any, flag_filter, min_quality, jobs=1, file=stdout, **kwargs):
    switch_backend("Agg")
    hide_stats_warnings(True)
    msg = "the levenshtein subprogram is in development and very finicky!"
    print("WARNING:", msg, file=stderr)
    samfilters = [flags, flags_any, flag_filter, min_quality]
    with AlignmentFile(bam) as alignment:
        bam_dict = load_bam_as_dict(alignment, samfilters)
    cms, report_rows = {}, []
    for chrom, entries in bam_dict.items():
        lds = calculate_chromosome_lds(chrom, entries)
        cm = generate_clustermap(lds)
        if cm is None:
            c1_names, c2_names = [], []
            c1_pval, c2_pval, silh_score = nan, nan, nan
            warn_about_unsupported_hierarchy(chrom)
        else:
            c1_names, c2_names, c1_pval, c2_pval, silh_score = get_clusters(
                lds, cm.dendrogram_row.linkage
            )
            cms[chrom] = cm
            if not isnan(silh_score):
                if (output_dir is not None) and (kmerscanner_file is not None):
                    generate_kmerscanner_file(
                        kmerscanner_file, set(c1_names), set(c2_names),
                        output_dir, chrom
                    )
            else:
                warn_about_unsupported_hierarchy(chrom)
        report_rows.append([
            chrom, len(c1_names), len(c2_names), silh_score, c1_pval, c2_pval
        ])
    report = generate_report(report_rows)
    batch_generate_pdfs(cms, report, output_dir)
    print(report.to_csv(sep="\t", index=False, na_rep="NA"))
    hide_stats_warnings(False)

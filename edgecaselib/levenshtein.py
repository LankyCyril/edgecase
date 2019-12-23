from sys import stdout, stderr
from warnings import filterwarnings, resetwarnings
from numba import njit
from pysam import AlignmentFile
from edgecaselib.formats import filter_bam
from numpy import zeros, array, uint32, uint8, log, nan, pi, isnan, allclose
from numpy import linspace, vstack, concatenate, unique, tile, triu
from collections import defaultdict
from edgecaselib.util import progressbar
from pandas import Series, DataFrame, read_csv, merge
from matplotlib.pyplot import switch_backend
from matplotlib.patches import Rectangle
from seaborn import clustermap
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_score
from scipy.stats import mannwhitneyu
from os import path
from statsmodels.stats.multitest import multipletests
from glob import glob
from re import search


CLUSTERMAP_FIGSIZE = (10, 10)
CLUSTERMAP_CMAP = "viridis_r"
CLUSTERMAP_VMAX = .15
LOG2PI1 = log(2 * pi) + 1


def load_bam_as_dict(alignment, samfilters):
    """Load BAM entries as dictionary: chrom -> qname -> (mappos, int8array)"""
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


def get_relative_read_ld(sra, A, srb, B, return_bases=False):
    """Calculate relative levenshtein distance between overlapping parts of two reads"""
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
    """Calculate pairwise relative levenshtein distances between all reads mapping to one chromosome"""
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
    """Generate clustermap of pairwise levenshtein distances between reads mapping to one chromosome"""
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


def loglikelihood(*, n, m, f, v, k):
    """Calculate log likelihood for arbitrary cluster"""
    if allclose(v, 0):
        return 0
    else:
        return m * (log(m) - log(n) - 0.5 * (f * log(v) + LOG2PI1)) + 0.5 * k


def cluster_loglikelihood(*, cluster, dataset_size, n_features, k):
    """Calculate log likelihood for subcluster of a cluster"""
    return loglikelihood(
        n=dataset_size, m=cluster.shape[0],
        f=n_features, v=cluster.var(), k=k
    )


def information_criterion(lds, labels, kind):
    """Calculate AIC or BIC for clustering"""
    n, f = lds.shape
    unique_labels = unique(labels)
    k = len(unique_labels)
    if kind == "AIC":
        penalty = - k * (f + 1) / 2 * 2
    elif kind == "BIC":
        penalty = - k * (f + 1) / 2 * log(n)
    else:
        raise ValueError("`kind` must be 'AIC' or 'BIC'")
    return penalty + sum(
        cluster_loglikelihood(
            cluster=lds.loc[labels==label].values,
            dataset_size=n, n_features=f, k=k
        )
        for label in unique_labels
    )


def get_mwu_pval(lds, ingroup_visual_indexer):
    """Calculate Mann-Whitney U p-value between within-cluster and out-of-cluster levenshtein distances"""
    ingroup_indexer = triu(ingroup_visual_indexer, k=1)
    outgroup_indexer = triu(~ingroup_visual_indexer, k=1)
    ingroup = lds.mask(~ingroup_indexer).values.flatten()
    outgroup = lds.mask(~outgroup_indexer).values.flatten()
    u, p = mannwhitneyu(
        ingroup[~isnan(ingroup)], outgroup[~isnan(outgroup)], alternative="less"
    )
    return p


def get_clusters(lds, linkage, min_cluster_size):
    """Find two major clusters; so far, only works if there is not more than one outlier read"""
    bic2k, k2silh, labels = {}, {}, {}
    for k in range(2, len(lds)):
        labels[k] = fcluster(linkage, k, criterion="maxclust")
        label_counts = dict(zip(*unique(labels[k], return_counts=True)))
        if sorted(label_counts.values())[-2] >= min_cluster_size:
            bic2k[information_criterion(lds, labels[k], "BIC")] = k
            k2silh[k] = silhouette_score(lds, labels[k])
            labels[k] = array([
                label if (label_counts[label]>=min_cluster_size) else nan
                for label in labels[k]
            ])
    if len(bic2k) == 0:
        return None, 0, nan, nan
    else:
        best_k = bic2k[max(bic2k)]
        best_silh, best_labels = k2silh[best_k], labels[best_k]
        if isnan(best_labels).any():
            min_label = min(label for label in best_labels if not isnan(label))
            best_labels = array([
                nan if isnan(label) else label-min_label+1
                for label in best_labels
            ])
        ingroup_axis = tile(best_labels, (len(best_labels), 1))
        ingroup_visual_indexer = (ingroup_axis == ingroup_axis.T)
        pval = get_mwu_pval(lds, ingroup_visual_indexer)
        used_k = best_labels[~isnan(best_labels)].max()
        return best_labels, used_k, best_silh, pval


def warn_about_unsupported_hierarchy(chrom):
    """Throw warning if impossible to satisfy clustering conditions"""
    msg = "({}): not implemented: complex clustering hierarchy or too few reads"
    print("\rWarning", msg.format(chrom), file=stderr)


def generate_kmerscanner_file(kmerscanner_file, names, labels, output_dir, chrom):
    """Annotate chromosomes for reads from different haplotypes"""
    kmerscanner_dat = read_csv(kmerscanner_file, sep="\t")
    name_to_label = DataFrame(
        data=[list(names), list(labels)], index=["#name", "label"]
    ).T
    haplo_dat = merge(
        kmerscanner_dat, name_to_label.dropna(), on="#name", how="inner"
    )
    haplo_dat["chrom"] = haplo_dat.apply(
        lambda row: "{}:haplotype {}".format(row["chrom"], int(row["label"])),
        axis=1
    )
    haplo_dat.drop(columns="label").to_csv(
        path.join(output_dir, chrom+".dat.gz"), compression="gzip",
        sep="\t", index=False
    )


def generate_report(report_rows, adj="bonferroni"):
    """Convert raw report to DataFrame and calculate adjusted p-values"""
    report = DataFrame(
        data=report_rows,
        columns=["#chrom", "cluster_count", "silhouette_score", "p"]
    )
    report["cluster_count"] = report["cluster_count"].astype(int)
    report["p_adjusted"] = multipletests(report["p"], method=adj)[1]
    return report


def draw_square(start, end, ax):
    """Draw fancy frame from (start, start) to (end, end)"""
    xy = (start, start)
    width = end - start + 1
    ax.add_patch(
        Rectangle(xy, width, width, fill=False, lw=5, ec="white", clip_on=False)
    )
    ax.add_patch(
        Rectangle(xy, width, width, fill=False, lw=1, ec="black", clip_on=False)
    )


def apply_mask(cm, labels):
    """Draw rectangles around within-cluster pairings"""
    labeled_names = Series(index=cm.data.index, data=labels)
    ordered_labeled_names = labeled_names.loc[cm.data2d.index].to_frame(
        name="label"
    )
    ordered_labeled_names["position"] = range(len(ordered_labeled_names))
    labels_groupby = ordered_labeled_names.groupby("label")
    starts, ends = (
        labels_groupby.min().iloc[:,0],
        labels_groupby.max().iloc[:,0]
    )
    for start, end in zip(starts, ends):
        draw_square(start, end, cm.ax_heatmap)


def generate_pdf(cm, silh_score, labels, output_dir, chrom, cmap=CLUSTERMAP_CMAP, vmax=CLUSTERMAP_VMAX):
    """Annotate clustermap figure and save to file"""
    cm.ax_col_dendrogram.clear()
    cm.ax_col_dendrogram.imshow(
        vstack([linspace(0, 1, 256)]*2), aspect="auto", cmap=cmap
    )
    gp = cm.ax_col_dendrogram.get_position()
    cm.ax_col_dendrogram.set(
        position=[
            gp.x0*.25+gp.x1*.75, gp.y0*.5+gp.y1*.46, gp.x1*.1-gp.x0*.1, .015
        ],
        zorder=float("inf")
    )
    cm.ax_col_dendrogram.text(
        x=-672, y=.8, va="center", ha="left", fontsize=19,
        s="Distance:  0"
    )
    cm.ax_col_dendrogram.text(
        x=272, y=.8, va="center", ha="left", fontsize=19,
        s="{}+".format(vmax)
    )
    cm.ax_col_dendrogram.set_axis_off()
    silh_text = "N/A" if isnan(silh_score) else "{:.3f}".format(silh_score)
    cm.ax_col_dendrogram.text(
        x=-672, y=3, va="top", ha="left", fontsize=19,
        s="Silhouette score: "+silh_text
    )
    cm.ax_col_dendrogram.text(
        x=-672, y=-4.9, va="top", ha="left", fontsize=19, s=chrom
    )
    apply_mask(cm, labels)
    filename = path.join(output_dir, chrom+".pdf")
    cm.fig.savefig(filename, bbox_inches="tight")


def hide_stats_warnings(state=True):
    """Prevent known harmless warnings from being printed to stderr"""
    if state:
        filterwarnings("ignore", message="invalid value encountered")
        filterwarnings(
            "ignore",
            message="looks suspiciously like an uncondensed distance matrix"
        )
        msg1 = "the levenshtein subprogram is in development!"
        msg2 = (
            "pairwise distance computation is O(n^2) " +
            "and not suited for large scale experiments"
        )
        print("WARNING:", msg1, file=stderr)
        print("WARNING:", msg2, file=stderr)
    else:
        resetwarnings()


def process_levenshtein_input(sequencedata, samfilters, output_dir):
    """Iterate over chromosomes and return pairwise levenshtein distances for reads mapped to them"""
    if path.isfile(sequencedata):
        with AlignmentFile(sequencedata) as alignment:
            bam_dict = load_bam_as_dict(alignment, samfilters)
        for chrom, entries in bam_dict.items():
            lds = calculate_chromosome_lds(chrom, entries)
            if output_dir:
                lds.to_csv(path.join(output_dir, chrom+"-matrix.tsv"), sep="\t")
            yield chrom, lds
    elif path.isdir(sequencedata):
        tsv_iterator = progressbar(
            glob(path.join(sequencedata, "*-matrix.tsv")),
            desc="Clustering", unit="chromsome"
        )
        for tsv in tsv_iterator:
            chrom_matcher = search(r'([^/]+)-matrix\.tsv', tsv)
            if chrom_matcher:
                chrom = chrom_matcher.group(1)
                lds = read_csv(tsv, sep="\t", index_col=0)
                is_lds_square = (lds.shape[0] == lds.shape[1])
                if (not is_lds_square) or (lds.index!=lds.columns).any():
                    msg_mask = "({}): malformed matrix? Skipping"
                    print("Warning", msg_mask.format(chrom), file=stderr)
                else:
                    yield chrom, lds
    else:
        raise IOError("Unknown type of input")


def main(sequencedata, min_cluster_size, kmerscanner_file, output_dir, flags, flags_any, flag_filter, min_quality, jobs=1, file=stdout, **kwargs):
    switch_backend("pdf")
    hide_stats_warnings(True)
    report_rows = []
    input_iterator = process_levenshtein_input(
        sequencedata, [flags, flags_any, flag_filter, min_quality], output_dir
    )
    for chrom, lds in input_iterator:
        cm = generate_clustermap(lds)
        if cm is not None:
            labels, k, silh_score, pval = get_clusters(
                lds, cm.dendrogram_row.linkage, min_cluster_size
            )
            if labels is not None:
                if output_dir:
                    generate_pdf(
                        cm, silh_score, labels, output_dir, chrom
                    )
                    if kmerscanner_file:
                        generate_kmerscanner_file(
                            kmerscanner_file, lds.index, labels,
                            output_dir, chrom
                        )
            else:
                warn_about_unsupported_hierarchy(chrom)
        else:
            k, silh_score, pval = 0, nan, nan
            warn_about_unsupported_hierarchy(chrom)
        report_rows.append([chrom, k, silh_score, pval])
    report = generate_report(report_rows)
    print(report.to_csv(sep="\t", index=False, na_rep="NA"))
    hide_stats_warnings(False)

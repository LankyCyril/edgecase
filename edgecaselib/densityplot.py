from sys import stdout, stderr
from edgecaselib.formats import load_kmerscan, interpret_flags
from matplotlib.pyplot import subplots, rc_context
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from seaborn import lineplot
from itertools import count
from tqdm import tqdm
from os import path
from re import split


def motif_subplots(nreads, chrom, max_mapq):
    """Prepare figure with subplots for each read"""
    page, axs = subplots(
        ncols=2, nrows=nreads, squeeze=False,
        figsize=(16, nreads*2/3), sharey=True, frameon=False,
        gridspec_kw={"hspace": 0, "wspace": 0, "width_ratios": (15, 1)}
    )
    # remove subplot borders:
    for ax in axs.flatten():
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.get_xaxis().get_major_formatter().set_scientific(False)
    # remove xticks for all reads except bottom one:
    for ax in axs[:-1, 0]:
        ax.set(xticks=[])
    # annotate chromosome and MAPQ:
    axs[-1, 0].set(xlabel="Chromosome {}".format(chrom))
    meta_twiny = axs[0, 1].twiny()
    meta_twiny.set(title="MAPQ", xlim=(0, max_mapq))
    for spine in meta_twiny.spines.values():
        spine.set_visible(False)
    axs[-1, 1].set(xticks=[])
    return page, axs


def plot_motif_densities(read_data, trace_ax, legend=False):
    """Plot traces for motif densities of one read"""
    trace_data = read_data.iloc[:,8:].T
    trace_data.columns = read_data["motif"]
    try:
        lineplot(data=trace_data, ax=trace_ax, legend=legend, dashes=False)
    except:
        pass
    if legend:
        trace_ax.legend(loc="lower left", bbox_to_anchor=(0, 1.4))
    return trace_data


def highlight_mapped_region(read_data, trace_data, name, trace_ax):
    """Plot a rectangle around mapped region of read"""
    leftmost_map = read_data["pos"].iloc[0]
    map_length = (
        trace_data.dropna().index.max() -
        read_data["clip_3prime"].iloc[0] -
        leftmost_map
    )
    trace_ax.add_patch(
        Rectangle(
            (leftmost_map, -.1), map_length, 1.2,
            edgecolor="gray", facecolor="none"
        )
    )
    trace_ax.text(
        leftmost_map+map_length/2, 1, name,
        verticalalignment="top", horizontalalignment="center"
    )


def plot_read_metadata(read_data, max_mapq, meta_ax):
    """Annotate additional SAM information for read"""
    mapq = read_data["mapq"].iloc[0]
    meta_ax.add_patch(
        Rectangle(
            (0, .2), mapq, .6,
            edgecolor="cornflowerblue", facecolor="cornflowerblue"
        )
    )
    meta_ax.set(xlim=(0, max_mapq))


def chromosome_motif_plot(binned_density_dataframe, chrom, max_mapq, title, no_align, anchors):
    """Render figure with all motif densities of all reads mapping to one chromosome"""
    names = binned_density_dataframe["name"].drop_duplicates()
    page, axs = motif_subplots(len(names), chrom, max_mapq)
    pos_range = binned_density_dataframe.columns[8:]
    for i, name, (trace_ax, meta_ax) in zip(count(), names, axs):
        read_data = binned_density_dataframe[
            binned_density_dataframe["name"]==name
        ]
        legend = "full" if (i==0) else False
        trace_data = plot_motif_densities(read_data, trace_ax, legend)
        if no_align:
            trace_ax.set(ylim=(-.2, 1.2), yticks=[], xticks=[])
        else:
            highlight_mapped_region(read_data, trace_data, name, trace_ax)
            trace_ax.set(
                xlim=(pos_range.min(), pos_range.max()),
                ylim=(-.2, 1.2), yticks=[]
            )
            plot_read_metadata(read_data, max_mapq, meta_ax)
        if anchors is not None:
            trace_ax.axvline(
                anchors.loc[chrom, "anchor"], -.2, 1.2,
                ls="--", lw=2, c="gray", alpha=.7
            )
    axs[0, 0].set(title=title)
    return page


def chromosome_natsort(chrom):
    """Natural order sorting that undestands chr1, chr10, chr14_K*, 7ptel etc"""
    keyoder = []
    for chunk in split(r'(\d+)', chrom): # stackoverflow.com/a/16090640
        if chunk.isdigit():
            keyoder.append(int(chunk))
        elif chunk == "":
            keyoder.append("chr")
        else:
            keyoder.append(chunk.lower())
    return keyoder


def plot_densities(densities, bin_size, title, no_align, anchors, file=stdout.buffer):
    """Plot binned densities as a heatmap"""
    max_mapq = max(d["mapq"].max() for d in densities.values())
    try:
        sorted_chromosomes = sorted(densities.keys(), key=chromosome_natsort)
    except Exception as e:
        msg = "natural sorting failed, pages will be sorted alphanumerically"
        print("Warning: " + msg, file=stderr)
        print("The error was: '{}'".format(e), file=stderr)
        sorted_chromosomes = sorted(densities.keys())
    sorted_densities_iterator = (
        (chrom, densities[chrom]) for chrom in sorted_chromosomes
    )
    decorated_densities_iterator = tqdm(
        sorted_densities_iterator, total=len(densities),
        desc="Plotting", unit="chromosome"
    )
    with rc_context({"figure.max_open_warning": len(densities)+2}):
        with PdfPages(file) as pdf:
            for chrom, binned_density_dataframe in decorated_densities_iterator:
                page = chromosome_motif_plot(
                    binned_density_dataframe, chrom,
                    max_mapq, title, no_align, anchors
                )
                pdf.savefig(page, bbox_inches="tight")


def main(dat, gzipped=None, flags=0, flag_filter=3844, min_quality=0, bin_size=100, title=None, no_align=False, reference=None, names=None, prime=None, file=stdout.buffer, **kwargs):
    """Dispatch data to subroutines"""
    densities = load_kmerscan(
        dat, gzipped, interpret_flags(flags), interpret_flags(flag_filter),
        min_quality, bin_size, no_align
    )
    if title is None:
        title = path.split(dat)[-1]
    if reference:
        if no_align:
            print("`reference` has no effect if no_align is True", file=stderr)
            anchors = None
        elif (names is not None) and (prime is not None):
            raise NotImplementedError
            #anchors, _ = get_anchors(reference, get_mainchroms(names))
            if prime == 5:
                anchors = anchors[["5prime"]]
            elif prime == 3:
                anchors = anchors[["3prime"]]
            else:
                raise ValueError("`prime` can only be 5 or 3")
            anchors.columns = ["anchor"]
        else:
            raise ValueError("For `reference`, must specify `names` & `prime`")
    else:
        anchors = None
    plot_densities(densities, bin_size, title, no_align, anchors, file)

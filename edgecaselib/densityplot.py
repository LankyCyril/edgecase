from sys import stdout
from edgecaselib.formats import load_index, load_kmerscan
from edgecaselib.formats import interpret_flags, FLAG_COLORS, explain_sam_flags
from edgecaselib.util import natsorted_chromosomes
from matplotlib.pyplot import subplots, rc_context
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from seaborn import lineplot
from itertools import count
from tqdm import tqdm
from os import path


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
    read_flag = read_data["flag"].iloc[0]
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
        leftmost_map+map_length/2, 1,
        name + "\n" + explain_sam_flags(read_flag),
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


def chromosome_motif_plot(binned_density_dataframe, ecx, chrom, max_mapq, no_align, title, flags, flag_filter, min_quality):
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
        plottable_flags = ecx.loc[ecx["rname"]==chrom, ["pos", "flag"]]
        for _, pos, flag in plottable_flags.itertuples():
            trace_ax.axvline(
                pos, -.2, 1.2, ls="--", lw=4,
                c=FLAG_COLORS[flag], alpha=.4
            )
    axs[0, 0].set(title="{}\n-f={} -F={} -q={}".format(
        title, flags, flag_filter, min_quality
    ))
    return page


def plot_densities(densities, ecx, bin_size, no_align, title, flags=None, flag_filter=None, min_quality=0, file=stdout.buffer):
    """Plot binned densities as a heatmap"""
    max_mapq = max(d["mapq"].max() for d in densities.values())
    sorted_chromosomes = natsorted_chromosomes(densities.keys())
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
                    binned_density_dataframe, ecx, chrom, max_mapq, no_align,
                    title, flags, flag_filter, min_quality
                )
                pdf.savefig(page, bbox_inches="tight")


def main(dat, gzipped=None, index=None, flags=0, flag_filter=0, min_quality=0, bin_size=100, no_align=False, title=None, file=stdout.buffer, **kwargs):
    """Dispatch data to subroutines"""
    ecx = load_index(index)
    densities = load_kmerscan(
        dat, gzipped, [flags, flag_filter, min_quality], bin_size, no_align
    )
    if title is None:
        title = path.split(dat)[-1]
    plot_densities(
        densities, ecx, bin_size, no_align,
        title, flags, flag_filter, min_quality,
        file
    )

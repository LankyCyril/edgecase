from sys import stdout
from numpy import linspace, array, mean, nan, concatenate, fromstring, full, vstack
from pandas import read_csv, DataFrame, concat, merge
from matplotlib.pyplot import subplots
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from seaborn import lineplot
from itertools import count
from tqdm import tqdm
from os import path


def binned(A, bins, func=mean):
    """Return array data compressed into bins (smoothed by func)"""
    coords = linspace(0, len(A), bins+1).astype(int)
    return array([
        func(A[start:end])
        for start, end in zip(coords, coords[1:])
    ])


def get_binned_density_dataframe(raw_densities, chrom, bin_size, no_align):
    """Subset to densities of reads on `chrom`, bin densities, convert to dataframe"""
    indexer = (raw_densities["chrom"] == chrom)
    densities_subset = raw_densities[indexer].reset_index(drop=True)
    if no_align:
        leftmost_pos = 0
    else:
        offsets = densities_subset["pos"] - densities_subset["clip_5prime"]
        leftmost_pos = offsets.min()
    # align density data to reference and bin it:
    density_arrays = []
    for entry in densities_subset.itertuples():
        unaligned_density = fromstring(entry.density, sep=",", dtype="float32")
        if no_align:
            padder = array([])
        else:
            padder = full(entry.pos - entry.clip_5prime - leftmost_pos, nan)
        aligned_density = concatenate([padder, unaligned_density])
        density_arrays.append(binned(
            aligned_density, bins=aligned_density.shape[0]/bin_size
        ))
    # pad densities on the right so they are all the same length:
    max_density_length = max(bd.shape[0] for bd in density_arrays)
    for i, binned_density in enumerate(density_arrays):
        padder = full(max_density_length - binned_density.shape[0], nan)
        density_arrays[i] = concatenate([binned_density, padder])
    naked_binned_density_dataframe = DataFrame(
        data=vstack(density_arrays),
        columns=[leftmost_pos + j * bin_size for j in range(max_density_length)]
    )
    binned_density_dataframe = concat(
        [densities_subset.iloc[:,:-1], naked_binned_density_dataframe], axis=1
    )
    return binned_density_dataframe.sort_values(
        by=["mapq", "name", "motif"], ascending=[False, True, True]
    )


def load_densities(dat, bin_size, no_align, each_once=True):
    """Load densities from dat file, split into dataframes per chromosome"""
    raw_densities = read_csv(dat, compression="gzip", sep="\t", escapechar="#")
    if each_once:
        raw_densities["length"] = raw_densities["density"].apply(
            lambda d: d.count(",")+1
        )
        groups = raw_densities[["name", "motif", "length"]].groupby(
            ["name", "motif"], as_index=False
        ).max()
        raw_densities = merge(groups, raw_densities).drop(columns="length")
    if no_align:
        raw_densities["chrom"] = "None"
    chromosome_iterator = tqdm(
        raw_densities["chrom"].drop_duplicates(), desc="Reading data",
        unit="chromosome"
    )
    return {
        chrom: get_binned_density_dataframe(
            raw_densities, chrom, bin_size, no_align
        )
        for chrom in chromosome_iterator
    }


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
    lineplot(data=trace_data, ax=trace_ax, legend=legend, dashes=False)
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


def chromosome_motif_plot(binned_density_dataframe, chrom, max_mapq, title, no_align):
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
    axs[0, 0].set(title=title)
    return page


def plot_densities(densities, bin_size, title, no_align, file=stdout.buffer):
    """Plot binned densities as a heatmap"""
    max_mapq = max(d["mapq"].max() for d in densities.values())
    chromosome_iterator = tqdm(
        densities.items(), total=len(densities),
        desc="Plotting", unit="chromosome"
    )
    with PdfPages(file) as pdf:
        for chrom, binned_density_dataframe in chromosome_iterator:
            page = chromosome_motif_plot(
                binned_density_dataframe, chrom, max_mapq, title, no_align
            )
            pdf.savefig(page, bbox_inches="tight")


def main(dat, bin_size=100, title=None, no_align=False, file=stdout.buffer, **kwargs):
    """Dispatch data to subroutines"""
    densities = load_densities(dat, bin_size=bin_size, no_align=no_align)
    if title is None:
        title = path.split(dat)[-1]
    plot_densities(
        densities, bin_size, title, no_align, file
    )

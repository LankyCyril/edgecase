from sys import stdout, stderr
from edgecaselib.formats import load_index, load_kmerscan
from edgecaselib.formats import FLAG_COLORS, explain_sam_flags, interpret_flags
from edgecaselib.formats import DEFAULT_MOTIF_COLORS, DEFAULT_MOTIF_HATCHES
from edgecaselib.formats import split_hatch
from edgecaselib.util import natsorted_chromosomes
from matplotlib.pyplot import subplots, rc_context, switch_backend
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch
from seaborn import lineplot
from itertools import count
from tqdm import tqdm
from os import path
from pandas import concat
from re import search


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
    trace_data = read_data.iloc[:,9:].T
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


def chromosome_motif_plot(binned_density_dataframe, ecx, chrom, max_mapq, no_align, title, samfilters):
    """Render figure with all motif densities of all reads mapping to one chromosome"""
    names = binned_density_dataframe["name"].drop_duplicates()
    page, axs = motif_subplots(len(names), chrom, max_mapq)
    pos_range = binned_density_dataframe.columns[9:]
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
                pos, -.2, 1.2, ls=":", lw=4, c=FLAG_COLORS[flag], alpha=.4
            )
    axs[0, 0].set(title="{}\n-f={} -g={} -F={} -q={}".format(title, *samfilters))
    return page


def make_decorated_densities_iterator(densities):
    """Order chromosomes and wrap with progress bar"""
    sorted_chromosomes = natsorted_chromosomes(densities.keys())
    sorted_densities_iterator = (
        (chrom, densities[chrom]) for chrom in sorted_chromosomes
    )
    return tqdm(
        sorted_densities_iterator, total=len(densities),
        desc="Plotting", unit="chromosome"
    )


def plot_exploded_densities(densities, ecx, no_align, title, samfilters, file=stdout.buffer):
    """Plot binned densities as line plots, one read at a time"""
    max_mapq = max(d["mapq"].max() for d in densities.values())
    decorated_densities_iterator = make_decorated_densities_iterator(densities)
    with rc_context({"figure.max_open_warning": len(densities)+2}):
        with PdfPages(file) as pdf:
            for chrom, binned_density_dataframe in decorated_densities_iterator:
                page = chromosome_motif_plot(
                    binned_density_dataframe, ecx, chrom, max_mapq, no_align,
                    title, samfilters
                )
                pdf.savefig(page, bbox_inches="tight")


def stack_motif_densities(binned_density_dataframe):
    """Prepare densities for plotting the stacked area chart"""
    motif_order = binned_density_dataframe[["motif", "total_count"]] \
        .drop_duplicates().sort_values(by="total_count") \
        ["motif"].drop_duplicates().values
    skinny_bdf = binned_density_dataframe[
        ["name", "motif"] + list(binned_density_dataframe.columns[9:])
    ]
    motif_bdfs = [
        skinny_bdf[skinny_bdf["motif"]==motif] \
            .set_index("name").drop(columns="motif")
        for motif in motif_order
    ]
    for i in range(1, len(motif_bdfs)):
        motif_bdfs[i] += motif_bdfs[i-1]
    for i in range(len(motif_bdfs)):
        motif_bdfs[i] = motif_bdfs[i].reset_index()
        motif_bdfs[i]["motif"] = motif_order[i]
    plottable_df = concat(motif_bdfs, axis=0).melt(
        id_vars=["name", "motif"], var_name="position", value_name="density"
    )
    return plottable_df, motif_order


def get_motif_data(plottable_df, motif):
    """Extract condensed averaged data for one motif (for fill_between)"""
    return plottable_df[plottable_df["motif"]==motif].dropna(how="any") \
        .groupby(["motif", "position"], as_index=False).mean()


def shorten_chrom_name(chrom):
    short_chrom_name_match = search(r'([Cc]hr)?[0-9XYM]+([pq]tel)?', chrom)
    if short_chrom_name_match:
        short_chrom_name = short_chrom_name_match.group()
    else:
        short_chrom_name = chrom[:5]
    if short_chrom_name != chrom:
        short_chrom_name = short_chrom_name + "..."
    return short_chrom_name


def plot_combined_density(binned_density_dataframe, ecx, chrom, title, m_clr, m_hch, target_anchor, is_q, ax):
    """Plot stacked area charts with bootstrapped CIs"""
    plottable_df, motif_order = stack_motif_densities(binned_density_dataframe)
    if len(motif_order) > len(m_clr):
        raise ValueError("Will not plot more motifs than available colors")
    else:
        ordered_m_clr = m_clr[::-1][:len(motif_order)]
        ordered_m_hch = m_hch[::-1][:len(motif_order)]
    for i in range(len(motif_order)+1):
        if i == len(motif_order):
            upper_motif_data = get_motif_data(plottable_df, motif_order[0])
            upper_motif_data["density"], color, hatch = 1, "silver", "*"
        else:
            upper_motif_data = get_motif_data(plottable_df, motif_order[i])
            color, hatch = ordered_m_clr[i], ordered_m_hch[i]
        if i == 0:
            lower_motif_data = upper_motif_data.copy()
            lower_motif_data["density"] = 0
        else:
            lower_motif_data = get_motif_data(plottable_df, motif_order[i-1])
        ax.fill_between(
            x=upper_motif_data["position"],
            y1=lower_motif_data["density"], y2=upper_motif_data["density"],
            color=color, hatch=hatch, alpha=.7
        )
    lineplot(
        data=plottable_df, x="position", y="density", hue="motif", legend=False,
        palette={m: "black" for m in set(plottable_df["motif"])},
        linewidth=.5, alpha=.7, ax=ax
    )
    indexer = (
        (ecx["rname"]==chrom) & (ecx["prime"]==(3 if is_q else 5)) &
        (ecx["flag"]==interpret_flags(target_anchor))
    )
    plottable_flags = ecx.loc[indexer, ["pos", "flag"]]
    for _, pos, flag in plottable_flags.itertuples():
        ax.axvline(pos, -.2, 1.2, ls=":", lw=4, c=FLAG_COLORS[flag], alpha=.4)
    short_chrom_name = shorten_chrom_name(chrom)
    ax.set(xlim=(
        plottable_df["position"].values.min(),
        plottable_df["position"].values.max()
    ))
    ax.set(xlabel="", ylim=(0, 1), ylabel=short_chrom_name, yticks=[])
    return motif_order, ordered_m_clr, ordered_m_hch


def align_subplots(ax2chrom, ecx, target_anchor, is_q):
    """Modify xlim of related axes to make their scales match"""
    prime = 3 if is_q else 5
    anchor_positions, left_spans, right_spans = {}, {}, {}
    for ax, chrom in ax2chrom.items():
        indexer = (
            (ecx["rname"]==chrom) & (ecx["prime"]==prime) &
            (ecx["flag"]==interpret_flags(target_anchor))
        )
        ecx_anchor_positions_set = set(ecx.loc[indexer, "pos"])
        if ecx_anchor_positions_set:
            anchor_position = ecx_anchor_positions_set.pop()
            anchor_positions[ax] = anchor_position
            xlim = ax.get_xlim()
            left_spans[ax] = anchor_position - xlim[0]
            right_spans[ax] = xlim[1] - anchor_position
        else:
            error_msk = "{} not found on {} prime for {} in ECX"
            raise ValueError(error_msk.format(target_anchor, prime, chrom))
    max_left_span = max(left_spans.values())
    max_right_span = max(right_spans.values())
    for ax, chrom in ax2chrom.items():
        ax.set(xlim=(
            anchor_positions[ax] - max_left_span,
            anchor_positions[ax] + max_right_span
        ))


def add_legend(motif_order, m_clr, m_hch, ax, exploded, is_q):
    """Add custom legend"""
    if exploded:
        raise NotImplementedError("Custom legend for exploded density plot")
    else:
        handles = [
            Patch(label=motif, fc=color, ec="black", hatch=hatch, alpha=.7)
            for motif, color, hatch in zip(motif_order, m_clr, m_hch)
        ][::-1]
        background_handle = [
            Patch(label="other", fc="silver", ec="black", hatch="*", alpha=.7)
        ]
        if is_q:
            loc="upper left"
        else:
            loc="upper right"
        ax.legend(handles=handles+background_handle, loc=loc, framealpha=1)
        ax.set(zorder=float("inf"))


def plot_densities(densities, ecx, title, m_clr, m_hch, target_anchor, is_q, file=stdout.buffer):
    """Plot binned densities as bootstrapped line plots, combined per chromosome"""
    decorated_densities_iterator = make_decorated_densities_iterator(densities)
    switch_backend("Agg")
    figure, axs = subplots(
        figsize=(16, len(densities)), gridspec_kw={"hspace": 1},
        nrows=len(densities), squeeze=False, frameon=False
    )
    ax2chrom = {}
    for (chrom, bdf), ax in zip(decorated_densities_iterator, axs[:, 0]):
        for spine in ["top", "left", "right"]:
            ax.spines[spine].set_visible(False)
        motif_order, ordered_m_clr, ordered_m_hch = plot_combined_density(
            bdf, ecx, chrom, title, m_clr, m_hch, target_anchor, is_q, ax=ax
        )
        ax.ticklabel_format(useOffset=False, style="plain")
        ax2chrom[ax] = chrom
    align_subplots(ax2chrom, ecx, target_anchor, is_q)
    add_legend(
        motif_order, ordered_m_clr, ordered_m_hch, axs[0, 0],
        exploded=False, is_q=is_q
    )
    figure.savefig(file, bbox_inches="tight", format="pdf")


def interpret_arguments(exploded, motif_colors, motif_hatches, samfilters, title, dat, no_align):
    """Parse and check arguments"""
    flags2set = lambda f: set(explain_sam_flags(interpret_flags(f)).split("|"))
    potential_target_anchors = {"tract_anchor", "ucsc_mask_anchor", "fork"}
    if exploded:
        for motif_argument in motif_colors, motif_hatches:
            if motif_argument is not None:
                message = "--motif-* with --exploded not fully implemented yet"
                print(message, file=stderr)
        target_anchor, is_q, m_clr, m_hch = None, None, None, None
    else:
        if no_align:
            raise NotImplementedError("--no-align without --exploded")
        if motif_colors:
            m_clr = motif_colors.split("|")
        else:
            m_clr = DEFAULT_MOTIF_COLORS
        if motif_hatches:
            m_hch = split_hatch(motif_hatches)
        else:
            m_hch = DEFAULT_MOTIF_HATCHES
        if len(m_clr) != len(m_hch):
            msg = "The numbers of motif colors and motif hatches do not match"
            raise ValueError(msg)
        flags, _, flag_filter, _ = samfilters
        if "is_q" in (flags2set(flags) - flags2set(flag_filter)):
            is_q = True
        elif "is_q" in (flags2set(flag_filter) - flags2set(flags)):
            is_q = False
        else:
            raise NotImplementedError("non-exploded on both arms")
        target_anchors = (
            potential_target_anchors &
            (flags2set(flags) - flags2set(flag_filter))
        )
        if len(target_anchors) == 1:
            target_anchor = target_anchors.pop()
        else:
            error_message = "non-exploded without a single identifiable anchor"
            raise NotImplementedError(error_message)
    if title is None:
        title = path.split(dat)[-1]
    return target_anchor, is_q, m_clr, m_hch, title


def main(dat, gzipped, index, flags, flags_any, flag_filter, min_quality, bin_size, exploded, motif_colors, motif_hatches, no_align, title, file=stdout.buffer, **kwargs):
    """Dispatch data to subroutines"""
    samfilters = [flags, flags_any, flag_filter, min_quality]
    target_anchor, is_q, m_clr, m_hch, title = interpret_arguments(
        exploded, motif_colors, motif_hatches,
        samfilters, title, dat, no_align
    )
    ecx = load_index(index)
    densities = load_kmerscan(dat, gzipped, samfilters, bin_size, no_align)
    if exploded:
        plot_exploded_densities(
            densities, ecx, no_align, title, samfilters, file
        )
    else:
        plot_densities(
            densities, ecx, title, m_clr, m_hch, target_anchor, is_q, file
        )

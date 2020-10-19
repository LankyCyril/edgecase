from sys import stdout
from edgecaselib.formats import load_index, load_kmerscan
from edgecaselib.formats import FLAG_COLORS, explain_sam_flags, interpret_flags
from edgecaselib.formats import DEFAULT_MOTIF_COLORS, PAPER_PALETTE, BGCOLOR
from collections import OrderedDict
from edgecaselib.util import natsorted_chromosomes, progressbar, motif_revcomp
from matplotlib.pyplot import subplots, rc_context, switch_backend
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch
from numpy import clip
from seaborn import lineplot
from itertools import count
from os import path, getenv
from pandas import concat
from re import search


buffer = getattr(stdout, "buffer", stdout)


__doc__ = """edgeCase densityplot: visualization of motif densities

Usage: {0} densityplot -x filename [-b integer] [-e] [--zoomed-in]
       {1}             [--palette palettespec] [--title string]
       {1}             [--n-boot integer]
       {1}             [-f flagspec] [-F flagspec] [-q integer]
       {1}             [-z] <dat>

Output:
    PDF file with motif densities visualized along chromosomal ends

Positional arguments:
    <dat>                         name of input kmerscanner file

Required options:
    -x, --index [filename]        location of the reference .ecx index

Options:
    -z, --gzipped                 input is gzipped (must specify if any of -qfF present
    -b, --bin-size [integer]      size of each bin in bp for visualization speedup [default: 100]
    --n-boot [integer]            number of bootstrap iterations for 95% confidence intervals [default: 1000]
    -e, --exploded                plot each read separately
    --zoomed-in                   plot taller traces, cut off pre-anchor regions
    --palette [palettespec]       custom palette for plotting motifs
    --title [string]              figure title (defaults to input filename)

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda bin_size: int(bin_size),
    lambda n_boot: int(n_boot),
]


DENSITYPLOT_FIGWIDTH = 17


def simplify_axes(ax, keep=(), keep_scientific=False):
    """Remove frame from axes, set tick formatting to plain"""
    for spine in ["top", "left", "right", "bottom"]:
        if spine not in keep:
            ax.spines[spine].set_visible(False)
    if not keep_scientific:
        ax.ticklabel_format(useOffset=False, style="plain")


def motif_subplots(nreads, chrom, max_mapq):
    """Prepare figure with subplots for each read"""
    page, axs = subplots(
        ncols=2, nrows=nreads, squeeze=False,
        figsize=(DENSITYPLOT_FIGWIDTH, nreads*2/3), sharey=True, frameon=False,
        gridspec_kw={"hspace": 0, "wspace": 0, "width_ratios": (15, 1)},
    )
    # remove subplot borders:
    for ax in axs.flatten():
        simplify_axes(ax)
    # remove xticks for all reads except bottom one:
    for ax in axs[:-1, 0]:
        ax.set(xticks=[])
    # annotate chromosome and MAPQ:
    axs[-1, 0].set(xlabel="Chromosome {}".format(chrom))
    meta_twiny = axs[0, 1].twiny()
    meta_twiny.set(title="MAPQ", xlim=(0, max_mapq))
    simplify_axes(meta_twiny)
    axs[-1, 1].set(xticks=[])
    return page, axs


def chromosome_subplots(nrows, zoomed_in=False):
    """Prepare figure with subplots for each chromosome end"""
    if zoomed_in:
        figsize=(DENSITYPLOT_FIGWIDTH*.6, nrows*3*.65)
        hspace = .35
    else:
        figsize=(DENSITYPLOT_FIGWIDTH*.7, nrows*1.5*.7)
        hspace = .7
    figure, axs = subplots(
        figsize=figsize, gridspec_kw={"hspace": hspace},
        nrows=nrows, squeeze=False, frameon=False,
    )
    for ax in axs[:, 0]:
        simplify_axes(ax, keep={"bottom"})
    return figure, axs


def plot_read_motif_densities(read_data, trace_ax, legend=False):
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
            edgecolor="gray", facecolor="none",
        )
    )
    trace_ax.text(
        leftmost_map+map_length/2, 1,
        name + "\n" + explain_sam_flags(read_flag),
        verticalalignment="top", horizontalalignment="center",
    )


def plot_read_metadata(read_data, max_mapq, meta_ax):
    """Annotate additional SAM information for read"""
    mapq = read_data["mapq"].iloc[0]
    meta_ax.add_patch(
        Rectangle(
            (0, .2), mapq, .6,
            edgecolor="cornflowerblue", facecolor="cornflowerblue",
        )
    )
    meta_ax.set(xlim=(0, max_mapq))


def chromosome_exploded_motif_plot(binned_density_dataframe, ecx, chrom, max_mapq, title, samfilters):
    """Render figure with all motif densities of all reads mapping to one chromosome"""
    names = binned_density_dataframe["name"].drop_duplicates()
    page, axs = motif_subplots(len(names), chrom, max_mapq)
    pos_range = binned_density_dataframe.columns[9:]
    for i, name, (trace_ax, meta_ax) in zip(count(), names, axs):
        read_data = binned_density_dataframe[
            binned_density_dataframe["name"]==name
        ]
        legend = "full" if (i==0) else False
        trace_data = plot_read_motif_densities(read_data, trace_ax, legend)
        highlight_mapped_region(read_data, trace_data, name, trace_ax)
        trace_ax.set(
            xlim=(pos_range.min(), pos_range.max()),
            ylim=(-.2, 1.2), yticks=[],
        )
        plot_read_metadata(read_data, max_mapq, meta_ax)
        plottable_flags = ecx.loc[ecx["rname"]==chrom, ["pos", "flag"]]
        for _, pos, flag in plottable_flags.itertuples():
            trace_ax.axvline(
                pos, -.2, 1.2, ls=":", lw=4, c=FLAG_COLORS[flag], alpha=.4,
            )
    axs[0, 0].set(
        title="{}\n-f={} -F={} -q={}".format(title, *samfilters),
    )
    return page


def make_decorated_densities_iterator(densities):
    """Order chromosomes and wrap with progress bar"""
    sorted_chromosomes = natsorted_chromosomes(densities.keys())
    sorted_densities_iterator = (
        (chrom, densities[chrom]) for chrom in sorted_chromosomes
    )
    return progressbar(
        sorted_densities_iterator, total=len(densities),
        desc="Plotting", unit="chromosome",
    )


def plot_exploded_densities(densities, ecx, title, samfilters, file=buffer):
    """Plot binned densities as line plots, one read at a time"""
    max_mapq = max(d["mapq"].max() for d in densities.values())
    decorated_densities_iterator = make_decorated_densities_iterator(densities)
    with rc_context({"figure.max_open_warning": len(densities)+2}):
        with PdfPages(file) as pdf:
            for chrom, binned_density_dataframe in decorated_densities_iterator:
                page = chromosome_exploded_motif_plot(
                    binned_density_dataframe, ecx, chrom, max_mapq,
                    title, samfilters,
                )
                pdf.savefig(page, bbox_inches="tight")


def stack_motif_densities(binned_density_dataframe):
    """Prepare densities for plotting the stacked area chart; return normal order but stack from bottom"""
    motif_order = (
        binned_density_dataframe[["motif", "score"]]
        .drop_duplicates().sort_values(by="score")["motif"]
        .drop_duplicates().values
    )
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
        id_vars=["name", "motif"], var_name="position", value_name="density",
    )
    return plottable_df, motif_order[::-1]


def get_motif_data(plottable_df, motif):
    """Extract condensed averaged data for one motif (for fill_between)"""
    return plottable_df[plottable_df["motif"]==motif].dropna(how="any") \
        .groupby(["motif", "position"], as_index=False).mean()


def shorten_chrom_name(chrom):
    """Shorten the chromosome name"""
    short_chrom_name_match = search(r'([Cc]hr)?[0-9XYM]+([pq]tel)?', chrom)
    if short_chrom_name_match:
        short_chrom_name = short_chrom_name_match.group()
    else:
        short_chrom_name = chrom[:5]
    if short_chrom_name != chrom:
        short_chrom_name = short_chrom_name + "..."
    return short_chrom_name


def fill_area(plottable_df, i, updated_palette, ax):
    """Fill areas on the area chart, from bottom to top"""
    reversed_motif_order = list(updated_palette.keys())[::-1]
    if i == len(reversed_motif_order):
        upper_motif_data = get_motif_data(plottable_df, reversed_motif_order[0])
        upper_motif_data["density"], color, hatch = 1, BGCOLOR, None
    else:
        motif = reversed_motif_order[i]
        color = updated_palette[motif]
        upper_motif_data = get_motif_data(plottable_df, motif)
        hatch = None
    if i == 0:
        lower_motif_data = upper_motif_data.copy()
        lower_motif_data["density"] = 0
    else:
        lower_motif_data = get_motif_data(
            plottable_df, reversed_motif_order[i-1],
        )
    ax.fill_between(
        x=upper_motif_data["position"],
        y1=lower_motif_data["density"], y2=upper_motif_data["density"],
        color=color, hatch=hatch, alpha=.7,
    )


def coverage_plot(plottable_df, motif_count, is_q, ax, y_offset=.1):
    """Plot read coverage under area chart"""
    covered_positions = plottable_df.loc[
        ~plottable_df["density"].isnull(), "position"
    ]
    coverage_df = covered_positions.value_counts() / motif_count
    coverage_df = coverage_df.astype(int).to_frame().reset_index()
    coverage_df.columns = ["position", "coverage"]
    coverage_df = coverage_df.sort_values(by="position")
    coverage_df["viz_coverage"] = 1 + y_offset + .01 * clip(
        coverage_df["coverage"], a_min=1, a_max=50,
    )
    ax.plot(
        [coverage_df["position"].min(), coverage_df["position"].max()],
        [coverage_df["viz_coverage"].max()]*2, ls="--", color="gray"
    )
    if is_q:
        x, ha = coverage_df["position"].max(), "left"
        ticks_mask, reads_label = " {}", "\n\nreads"
    else:
        x, ha = coverage_df["position"].min(), "right"
        ticks_mask, reads_label = "{} ", "reads\n\n"
    ax.text(
        x=x, y=coverage_df["viz_coverage"].max(), va="center", ha=ha,
        s=ticks_mask.format(coverage_df["coverage"].max())
    )
    ax.text(x=x, y=1+y_offset, va="center", ha=ha, s=ticks_mask.format(1))
    ax.text(
        x=x, y=.5*(1+y_offset+coverage_df["viz_coverage"].max()),
        va="center", ha=ha, rotation=90, s=reads_label,
    )
    ax.fill_between(
        x=coverage_df["position"], y1=1+y_offset,
        y2=coverage_df["viz_coverage"], step="pre", color="gray", alpha=.4,
    )
    return coverage_df["viz_coverage"].max()


def generate_updated_palette(palette, motif_order):
    """Given user palette and motifs in density file, generate palette that satisfies both"""
    updated_palette = OrderedDict()
    if palette is None:
        if len(motif_order) > len(DEFAULT_MOTIF_COLORS):
            error_msk = "Cannot plot over {} motifs with default palette; {}"
            raise ValueError(error_msk.format(
                len(DEFAULT_MOTIF_COLORS), "please provide a custom --palette",
            ))
        else:
            for motif, color in zip(motif_order, DEFAULT_MOTIF_COLORS):
                updated_palette[motif] = color
        updated_motif_order = motif_order
    else:
        for motif in motif_order:
            if motif in palette:
                updated_palette[motif] = palette[motif]
        updated_motif_order = [
            motif for motif in motif_order if motif in palette
        ]
    assert set(updated_palette.keys()) == set(updated_motif_order)
    return updated_palette, updated_motif_order


def plot_combined_density(binned_density_dataframe, n_boot, ecx, title, palette, target_anchor, is_q, display_chrom_name, ecx_chrom_name, zoomed_in, ax):
    """Plot stacked area charts with bootstrapped CIs"""
    plottable_df, motif_order = stack_motif_densities(binned_density_dataframe)
    updated_palette, updated_motif_order = generate_updated_palette(
        palette, motif_order,
    )
    for i in range(len(updated_motif_order)+1):
        fill_area(plottable_df, i, updated_palette, ax)
    lineplot(
        data=plottable_df, x="position", y="density", hue="motif", legend=False,
        palette={m: "black" for m in set(plottable_df["motif"])},
        n_boot=n_boot, linewidth=.5, alpha=.7, ax=ax,
    )
    lineplot(
        x=plottable_df["position"].drop_duplicates().sort_values(), y=1,
        legend=False, color="black", linewidth=1, alpha=1, ax=ax,
    )
    if zoomed_in:
        ymax = coverage_plot(plottable_df, len(motif_order), is_q, ax)
    else:
        ymax = 1
    indexer = (
        (ecx["rname"]==ecx_chrom_name) & (ecx["prime"]==(3 if is_q else 5)) &
        (ecx["flag"]==interpret_flags(target_anchor))
    )
    plottable_flags = ecx.loc[indexer, ["pos", "flag"]]
    for _, pos, flag in plottable_flags.itertuples():
        ax.axvline(pos, 0, 1/ymax, ls=":", lw=4, c=FLAG_COLORS[flag], alpha=.4)
    position_values = plottable_df["position"].values
    ax.set(xlim=(position_values.min(), position_values.max()))
    ax.set(xlabel="", ylim=(0, ymax), ylabel=display_chrom_name, yticks=[])
    return updated_palette


def align_subplots(ax2chrom, ecx, target_anchor, is_q, unit_adjustment=1e6):
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
    max_left_span = str(getenv("PAPER_LEFT_SPAN"))
    if max_left_span.isdigit():
        max_left_span = int(max_left_span)
    else:
        max_left_span = max(left_spans.values())
    max_right_span = str(getenv("PAPER_RIGHT_SPAN"))
    if max_right_span.isdigit():
        max_right_span = int(max_right_span)
    else:
        max_right_span = max(right_spans.values())
    for ax, chrom in ax2chrom.items():
        ax.set(xlim=(
            anchor_positions[ax] - max_left_span,
            anchor_positions[ax] + max_right_span,
        ))
        if unit_adjustment:
            ax.set(
                xticklabels=[
                    "{:.3f}".format(int(xt) / unit_adjustment)
                    for xt in ax.get_xticks().tolist()
                ]
            )


def add_legend(updated_palette, ax, exploded, is_q):
    """Add custom legend"""
    if exploded:
        raise NotImplementedError("Custom legend for exploded density plot")
    else:
        handles = [
            Patch(label=motif, fc=color, ec="black", hatch=None, alpha=.7)
            for motif, color in reversed(updated_palette.items())
        ][::-1]
        background_handle = [Patch(
            label="background", fc=BGCOLOR, ec="black", hatch=None, alpha=.7,
        )]
        if is_q:
            loc="upper left"
        else:
            loc="upper right"
        ax.legend(handles=handles+background_handle, loc=loc, framealpha=.9)
        ax.set(zorder=float("inf"))


def plot_density_scale(ax):
    bar_position = 1 + .002 * DENSITYPLOT_FIGWIDTH
    tick_width = .0003 * DENSITYPLOT_FIGWIDTH
    xtrace = [
        bar_position+tick_width, bar_position,
        bar_position, bar_position+tick_width,
    ]
    ax.plot(
        xtrace, [1, 1, 0, 0], transform=ax.transAxes,
        color="black", linewidth=1, clip_on=False,
    )
    ax.text(
        x=bar_position+tick_width, y=1, transform=ax.transAxes,
        s=" 100%", ha="left", va="center",
    )
    ax.text(
        x=bar_position+tick_width, y=0, transform=ax.transAxes,
        s=" 0%", ha="left", va="center",
    )
    ax.text(
        x=bar_position-tick_width, y=.5, transform=ax.transAxes, rotation=90,
        s="density", ha="right", va="center",
    )


def plot_densities(densities, n_boot, ecx, title, palette, legend, target_anchor, is_q, zoomed_in, file=buffer):
    """Plot binned densities as bootstrapped line plots, combined per chromosome"""
    decorated_densities_iterator = make_decorated_densities_iterator(densities)
    switch_backend("Agg")
    figure, axs = chromosome_subplots(len(densities), zoomed_in)
    ax2chrom = {}
    for (chrom, bdf), ax in zip(decorated_densities_iterator, axs[:, 0]):
        ecx_chrom_name, comment, *_ = chrom.split(":", 1) + [None]
        if zoomed_in:
            short_chrom_name = ecx_chrom_name
        else:
            short_chrom_name = shorten_chrom_name(ecx_chrom_name)
        if comment:
            display_chrom_name = short_chrom_name + "\n" + comment
        else:
            display_chrom_name = short_chrom_name
        updated_palette = plot_combined_density(
            bdf, n_boot, ecx, title, palette, target_anchor, is_q, ax=ax,
            zoomed_in=zoomed_in, display_chrom_name=display_chrom_name,
            ecx_chrom_name=ecx_chrom_name,
        )
        ax2chrom[ax] = ecx_chrom_name
    align_subplots(ax2chrom, ecx, target_anchor, is_q)
    if legend:
        add_legend(
            updated_palette, axs[0, 0], exploded=False, is_q=is_q,
        )
        plot_density_scale(axs[0, 0])
    if title:
        axs[0, 0].set(title=title)
    axs[-1, 0].set(xlabel="Mbp")
    figure.savefig(file, bbox_inches="tight", format="pdf")


def interpret_target(samfilters):
    """For non-exploded densityplots, infer which arm to plot and which anchor to center around"""
    flags2set = lambda f: set(explain_sam_flags(interpret_flags(f)).split("|"))
    potential_target_anchors = {"tract_anchor", "mask_anchor", "fork"}
    flags, flag_filter, _ = samfilters
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
    return is_q, target_anchor


def generate_paper_palette(paper_palette, is_q):
    """Return paper_palette, reverse complement motifs if is_q is False"""
    if is_q is True:
        return paper_palette
    elif is_q is False:
        palette = OrderedDict()
        for motif, color in paper_palette.items():
            palette[motif_revcomp(motif)] = color
        return palette
    else:
        raise NotImplementedError("--palette 'paper' for unknown arm")


def interpret_arguments(palette, exploded, zoomed_in, samfilters, title, dat):
    """Parse and check arguments"""
    if exploded:
        if zoomed_in:
            raise NotImplementedError("--exploded with --zoomed-in")
        target_anchor, is_q = None, None
    else:
        is_q, target_anchor = interpret_target(samfilters)
    if palette is None:
        legend = not zoomed_in
    elif palette in {"paper|legend=False", "paper|legend=True", "paper"}:
        legend = "False" not in palette
        palette = generate_paper_palette(PAPER_PALETTE, is_q)
    else:
        if exploded:
            raise NotImplementedError("--palette with --exploded")
        interpreted_palette = OrderedDict()
        legend = not zoomed_in
        for palette_field in palette.split("|"):
            if palette_field.startswith("legend="):
                if palette_field[7:].lower() == "true":
                    legend = True
                elif palette_field[7:].lower() == "false":
                    legend = False
                else:
                    raise ValueError("Uknown syntax: " + palette_field)
            elif "=" in palette_field:
                motif, color = palette_field.split("=", 1)
                interpreted_palette[motif] = color
            else:
                raise ValueError("Uknown syntax: " + palette_field)
        if interpreted_palette:
            palette = interpreted_palette
        else:
            palette = None
    return target_anchor, is_q, palette, legend, (title or path.split(dat)[-1])


def main(dat, gzipped, index, flags, flag_filter, min_quality, bin_size, n_boot, exploded, zoomed_in, palette, title, file=buffer, **kwargs):
    """Dispatch data to subroutines"""
    samfilters = [flags, flag_filter, min_quality]
    target_anchor, is_q, palette, legend, title = interpret_arguments(
        palette, exploded, zoomed_in, samfilters, title, dat,
    )
    ecx = load_index(index)
    densities = load_kmerscan(dat, gzipped, samfilters, bin_size)
    if exploded:
        plot_exploded_densities(
            densities, ecx, title, samfilters, file,
        )
    else:
        plot_densities(
            densities, n_boot, ecx, title, palette, legend, target_anchor,
            is_q, zoomed_in, file,
        )

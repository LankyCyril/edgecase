from sys import stdout
from edgecaselib.formats import load_index, load_kmerscan
from edgecaselib.formats import FLAG_COLORS, explain_sam_flags, interpret_flags
from edgecaselib.formats import TOL_COLORSCHEME, PAPER_PALETTE, BGCOLOR
from collections import OrderedDict
from edgecaselib.util import natsorted_chromosomes, progressbar, motif_revcomp
from matplotlib.pyplot import subplots, switch_backend, cm
from matplotlib.patches import Patch
from numpy import clip, arange
from seaborn import lineplot
from os import path
from pandas import concat
from re import search
from pickle import dump
buffer = getattr(stdout, "buffer", stdout)


__doc__ = """edgeCase densityplot: visualization of motif densities

Usage: {0} densityplot -x filename [-b integer] [--plot-coverage]
       {1}             [--palette palettespec] [--title string]
       {1}             [--n-boot integer] [--chroms-to-plot string]
       {1}             [-f flagspec]... [-F flagspec]... [-q integer]
       {1}             [--figwidth-inches float] [--outfmt string] [-z] <dat>

Output:
    PDF file with motif densities visualized along chromosomal ends

Positional arguments:
    <dat>                         name of input kmerscanner file

Required options:
    -x, --index [filename]        location of the reference .ecx index

Options:
    -z, --gzipped                 input is gzipped (must specify if any of -qfF present)
    -b, --bin-size [integer]      size of each bin in bp (overrides bin size in <dat>)
    --n-boot [integer]            number of bootstrap iterations for 95% confidence intervals [default: 1000]
    --palette [palettespec]       custom palette for plotting motifs
    --title [string]              figure title (defaults to input filename)
    --chroms-to-plot [string]     if set, plot chromosomes from this comma-separated list unconditionally
    --plot-coverage               plot coverage by telomeric reads on each arm
    --figwidth-inches [float]     width of figure in inches [default: 13]
    --outfmt [string]             output format (pdf, pkl) [default: pdf]

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda bin_size: None if bin_size is None else int(bin_size),
    lambda n_boot: int(n_boot),
    lambda figwidth_inches: float(figwidth_inches),
]

__docopt_tests__ = {
    lambda outfmt:
        outfmt in {"pdf", "pkl"}:
            "--outfmt must be either 'pdf' or 'pkl'",
}


def simplify_axes(ax, keep=(), keep_scientific=False):
    """Remove frame from axes, set tick formatting to plain"""
    for spine in ["top", "left", "right", "bottom"]:
        if spine not in keep:
            ax.spines[spine].set_visible(False)
    if not keep_scientific:
        ax.ticklabel_format(useOffset=False, style="plain")


def chromosome_subplots(nrows, figwidth_inches, plot_coverage=False):
    """Prepare figure with subplots for each chromosome end"""
    if plot_coverage:
        figsize = (figwidth_inches*.7, nrows*2.2*.7*figwidth_inches/15)
        hspace = .3
    else:
        figsize = (figwidth_inches*.7, nrows*1.5*.7*figwidth_inches/15)
        hspace = .7
    figure, axs = subplots(
        figsize=figsize, gridspec_kw={"hspace": hspace},
        nrows=nrows, squeeze=False, frameon=False,
    )
    for ax in axs[:, 0]:
        simplify_axes(ax, keep={"bottom"})
    return figure, axs


def make_decorated_densities_iterator(densities, chroms_to_plot=None):
    """Order chromosomes and wrap with progress bar"""
    if chroms_to_plot:
        sorted_chromosomes = natsorted_chromosomes(
            set(densities.keys()) | set(chroms_to_plot.split(","))
        )
    else:
        sorted_chromosomes = natsorted_chromosomes(densities.keys())
    sorted_densities_iterator = (
        (chrom, densities.get(chrom)) for chrom in sorted_chromosomes
    )
    decorated_densities_iterator = progressbar(
        sorted_densities_iterator, total=len(sorted_chromosomes),
        desc="Plotting", unit="chromosome",
    )
    return decorated_densities_iterator, len(sorted_chromosomes)


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
    return (
        plottable_df[plottable_df["motif"]==motif]
        .dropna(how="any")
        .groupby(["motif", "position"], as_index=False)
        .mean()
    )


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


def coverage_plot(plottable_df, motif_count, is_q, ax, y_offset=.025):
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
    max_viz_coverage = coverage_df["viz_coverage"].max()
    poly, = ax.fill(
        coverage_df["position"], coverage_df["viz_coverage"], color="none",
    )
    gradient = arange(
        coverage_df["coverage"].min(), coverage_df["coverage"].max(), .1,
    )
    img = ax.imshow(
        gradient.reshape(gradient.size, 1), aspect="auto", origin="lower",
        cmap=cm.bone, alpha=.6, vmin=-10, vmax=60, extent=[
            coverage_df["position"].min(),
            coverage_df["position"].max(),
            coverage_df["viz_coverage"].min(),
            coverage_df["viz_coverage"].max(),
        ],
    )
    img.set_clip_path(poly)
    return max_viz_coverage


def generate_updated_palette(palette, motif_order):
    """Given user palette and motifs in density file, generate palette that satisfies both"""
    updated_palette = OrderedDict()
    if palette is None:
        if len(motif_order) > len(TOL_COLORSCHEME):
            error_msk = "Cannot plot over {} motifs with default palette; {}"
            raise ValueError(error_msk.format(
                len(TOL_COLORSCHEME), "please provide a custom --palette",
            ))
        else:
            for motif, color in zip(motif_order, TOL_COLORSCHEME.values()):
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


def plot_anchors(ecx_chrom_name, ecx, target_anchor, is_q, ax):
    """Plot anchor positions over densities"""
    indexer = (
        (ecx["rname"]==ecx_chrom_name) & (ecx["prime"]==(3 if is_q else 5)) &
        (ecx["flag"]==interpret_flags(target_anchor))
    )
    plottable_flags = ecx.loc[indexer, ["pos", "flag"]]
    for _, pos, flag in plottable_flags.itertuples():
        ax.plot(
            [pos, pos], [0, 1], ls="--", lw=3, c=FLAG_COLORS[flag], alpha=.4,
        )
    return ecx.loc[indexer, "pos"].min(), ecx.loc[indexer, "pos"].max()


def plot_combined_density(binned_density_dataframe, n_boot, ecx, title, palette, target_anchor, is_q, ecx_chrom_name, plot_coverage, ax):
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
    if plot_coverage:
        ymax = coverage_plot(plottable_df, len(motif_order), is_q, ax)
    else:
        ymax = 1
    position_values = plottable_df["position"].values
    ax.set(xlim=(position_values.min(), position_values.max()))
    ax.set(xlabel="", ylim=(0, ymax), ylabel="", yticks=[])
    return updated_palette


def align_subplots(ax2chrom, ax2ylabel, ecx, target_anchor, is_q, unit_adjustment=None, unit_fmt=",.0f"):
    """Modify xlim of related axes to make their scales match"""
    prime = 3 if is_q else 5
    anchor_positions, left_spans, right_spans = {}, {}, {}
    ymax = 0
    for ax, chrom in ax2chrom.items():
        ymax = max(ymax, ax.get_ylim()[1])
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
    for ax, chrom in ax2chrom.items():
        ax.set(
            xlim=(
                anchor_positions[ax] - max(left_spans.values()),
                anchor_positions[ax] + max(right_spans.values()),
            ),
            ylim=(0, ymax),
        )
        xmin, xmax = ax.get_xlim()
        ax.text(
            x=xmin-(xmax-xmin)/100, y=.5, transform=ax.transData,
            s=ax2ylabel[ax], ha="right", va="center", rotation=90,
        )
        if unit_adjustment:
            ax.set(
                xticklabels=[
                    format(int(xt) / unit_adjustment, unit_fmt)
                    for xt in ax.get_xticks().tolist()
                ]
            )


def add_motifs_legend(updated_palette, ax, is_q):
    """Add custom legend"""
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


def plot_density_scale(ax, figwidth_inches):
    bar_position = 1 + .002 * figwidth_inches
    tick_width = .0004 * figwidth_inches
    ymax = ax.get_ylim()[1]
    xtrace = [
        bar_position+tick_width, bar_position,
        bar_position, bar_position+tick_width,
    ]
    ax.plot(
        xtrace, [1/ymax, 1/ymax, 0, 0], transform=ax.transAxes,
        color="black", linewidth=1, clip_on=False,
    )
    ax.text(
        x=bar_position+tick_width, y=1/ymax, transform=ax.transAxes,
        s=" 100%", ha="left", va="center",
    )
    ax.text(
        x=bar_position+tick_width, y=0, transform=ax.transAxes,
        s=" 0%", ha="left", va="center",
    )
    ax.text(
        x=bar_position-tick_width, y=.5/ymax, transform=ax.transAxes,
        s="density", ha="right", va="center",
        rotation=90, size=figwidth_inches*2/3,
    )


def format_chrom(chrom):
    """Get printable options based on value of `chrom`"""
    ecx_chrom_name, comment, *_ = chrom.split(":", 1) + [None]
    short_chrom_name = shorten_chrom_name(ecx_chrom_name)
    if comment:
        display_chrom_name = short_chrom_name + "\n" + comment
    else:
        display_chrom_name = short_chrom_name
    return ecx_chrom_name, short_chrom_name, display_chrom_name


def plot_densities(densities, n_boot, ecx, title, palette, motifs_legend, density_legend, target_anchor, is_q, chroms_to_plot, figwidth_inches, plot_coverage, unit="Kbp"):
    """Plot binned densities as bootstrapped line plots, combined per chromosome"""
    decorated_densities_iterator, n_axes = make_decorated_densities_iterator(
        densities, chroms_to_plot,
    )
    figure, axs = chromosome_subplots(n_axes, figwidth_inches, plot_coverage)
    ax2chrom, ax2ylabel = {}, {}
    for (chrom, bdf), ax in zip(decorated_densities_iterator, axs[:, 0]):
        ecx_chrom_name, short_chrom_name, display_chrom_name = format_chrom(
            chrom,
        )
        if bdf is None:
            ax.set(xlabel="", ylim=(0, 1), ylabel=display_chrom_name, yticks=[])
        else:
            updated_palette = plot_combined_density(
                bdf, n_boot, ecx, title, palette, target_anchor, is_q, ax=ax,
                ecx_chrom_name=ecx_chrom_name, plot_coverage=plot_coverage,
            )
        xlim = plot_anchors(ecx_chrom_name, ecx, target_anchor, is_q, ax)
        if bdf is None:
            ax.set(xlim=(xlim[0]-1, xlim[1]+1))
        ax2chrom[ax], ax2ylabel[ax] = ecx_chrom_name, display_chrom_name
    try:
        unit_adjustment = {"Kbp": 1e3, "Mbp": 1e6, "bp": None}[unit]
        unit_fmt = {"Kbp": ",.0f", "Mbp": ",.3f", "bp": ",.0f"}[unit]
    except KeyError:
        raise ValueError("unit", unit)
    align_subplots(
        ax2chrom, ax2ylabel, ecx, target_anchor, is_q,
        unit_adjustment=unit_adjustment, unit_fmt=unit_fmt,
    )
    if motifs_legend:
        add_motifs_legend(updated_palette, axs[0, 0], is_q=is_q)
    if density_legend:
        plot_density_scale(axs[0, 0], figwidth_inches)
    axs[0, 0].set(title=title)
    axs[-1, 0].set(xlabel=unit)
    return figure


def interpret_target(samfilters):
    """Infer which arm to plot and which anchor to center around"""
    flags2set = lambda f: set(explain_sam_flags(interpret_flags(f)))
    flags, flag_filter, _ = samfilters
    if "is_q" in (flags2set(flags) - flags2set(flag_filter)):
        is_q = True
    elif "is_q" in (flags2set(flag_filter) - flags2set(flags)):
        is_q = False
    else:
        raise NotImplementedError("Input contains reads on both arms")
    target_anchors = (
        {"tract_anchor", "mask_anchor", "fork"} &
        (flags2set(flags) - flags2set(flag_filter))
    )
    if len(target_anchors) == 1:
        target_anchor = target_anchors.pop()
    else:
        error_message = f"Multiple anchors to plot: {sorted(target_anchor)}"
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


def interpret_arguments(palette, chroms_to_plot, samfilters, title, outfmt, dat):
    """Parse and check arguments"""
    is_q, target_anchor = interpret_target(samfilters)
    PAPER_PALETTE_AS_PASSED_ARGS = {
        "paper", "paper|legend=none", "paper|legend=full",
        "paper|legend=density", "paper|legend=motifs",
    }
    if palette is None:
        motifs_legend, density_legend = True, True
    elif palette in PAPER_PALETTE_AS_PASSED_ARGS:
        full_legend = "legend=full" in palette
        motifs_legend = ("legend=motifs" in palette) or full_legend
        density_legend = ("legend=density" in palette) or full_legend
        palette = generate_paper_palette(PAPER_PALETTE, is_q)
    else:
        interpreted_palette = OrderedDict()
        motifs_legend, density_legend = True, True
        for palette_field in palette.split("|"):
            if palette_field.startswith("legend="):
                spec = palette_field[7:]
                if spec in {"full", "motifs", "density", "none"}:
                    motifs_legend = spec in {"full", "motifs"}
                    density_legend = spec in {"full", "density"}
                else:
                    raise ValueError("Uknown syntax: " + palette_field)
            elif "=" in palette_field:
                motif, color = palette_field.split("=", 1)
                interpreted_palette[motif] = color
            else:
                raise ValueError("Uknown syntax: " + palette_field)
        palette = interpreted_palette or None
    return (
        target_anchor, is_q, palette, motifs_legend, density_legend,
        (title or path.split(dat)[-1]),
    )


def main(dat, index, gzipped, flags, flag_filter, min_quality, bin_size, n_boot, palette, title, chroms_to_plot, plot_coverage, outfmt, figwidth_inches, file=buffer, **kwargs):
    """Dispatch data to subroutines"""
    samfilters = [flags, flag_filter, min_quality]
    if palette and (palette[0] == "'") and (palette[-1] == "'"):
        palette = palette[1:-1]
    _intrg = interpret_arguments(
        palette, chroms_to_plot, samfilters, title, outfmt, dat,
    )
    target_anchor, is_q, palette, motifs_legend, density_legend, title = _intrg
    ecx = load_index(index)
    densities = load_kmerscan(dat, gzipped, samfilters, bin_size)
    switch_backend("Agg")
    figure = plot_densities(
        densities, n_boot, ecx, title, palette, motifs_legend, density_legend,
        target_anchor, is_q, chroms_to_plot, figwidth_inches, plot_coverage,
    )
    if outfmt == "pdf":
        figure.savefig(file, bbox_inches="tight", format="pdf")
    elif outfmt == "pkl":
        dump(figure, file)
    else:
        raise ValueError(f"--outfmt={outfmt}")
    return 0

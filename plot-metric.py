#!/usr/bin/env python3
from sys import stdout
from re import search
from argparse import ArgumentParser
from numpy import linspace, array, mean, nan, concatenate, fromiter
from pandas import Series, concat
from matplotlib.pyplot import switch_backend, subplots
from seaborn import heatmap

USAGE = "python3 {} [options] txt > png".format(__file__)

ARG_RULES = {
    ("txt",): {
        "help": "name of input (each line is 'name\\tvalue\\tvalue\\tvalue...')"
    },
    ("-b", "--bin-size"): {
        "help": "size of each bin in bp for visualization speedup (100)",
        "default": 100, "type": int, "metavar": "B"
    },
    ("-a", "--align"): {
        "help": "alignment of visualized reads (left)",
        "default": "left", "metavar": "A"
    },
    ("-s", "--figsize"): {
        "help": "figure size (16x9)",
        "default": "16x9", "metavar": "S"
    },
    ("-p", "--palette"): {
        "help": "heatmap palette (viridis)",
        "default": "viridis", "metavar": "P"
    },
    ("--hide-names",): {
        "help": "hide names of reads on heatmap",
        "action": "store_true"
    },
    ("--title",): {
        "help": "figure title (defaults to input filename)"
    },
    ("--xtick-density",): {
        "help": "xtick density compared to heatmap default (.05)",
        "default": .05, "type": float, "metavar": "X"
    }
}


def binned(A, bins, func=mean):
    """Return array data compressed into bins (smoothed by func)"""
    coords = linspace(0, len(A), bins+1).astype(int)
    return array([
        func(A[start:end])
        for start, end in zip(coords, coords[1:])
    ])


def load_metrics(txt, bin_size, align):
    """Load metrics from a text file, bin and convert into dataframe"""
    read_metrics = {}
    # load and bin all metrics first:
    with open(txt, mode="rt") as handle:
        for line in handle:
            name, *metrics = line.strip().split("\t")
            metrics_array = fromiter(map(
                lambda s: float(s) if s!="" else nan,
                metrics
            ), dtype="float32")
            read_metrics[name] = binned(
                metrics_array,
                bins=len(metrics)/bin_size
            )
    # coerce to same lengths and convert into DataFrame:
    maxlen = max(len(d) for d in read_metrics.values())
    rows = []
    for name, metrics in read_metrics.items():
        padder = array([nan] * (maxlen - len(metrics)))
        if align == "left":
            padded_metrics = concatenate([metrics, padder])
        elif align == "right":
            padded_metrics = concatenate([padder, metrics])
        rows.append(Series(padded_metrics, name=name))
    return concat(rows, axis=1).T


def plot_metrics(metrics, figsize, palette, hide_names, bin_size, align, xtick_density, title="", png=stdout.buffer):
    """Plot binned metrics as a heatmap"""
    switch_backend("Agg")
    width, height = tuple(map(int, figsize.split("x")))
    figure, ax = subplots(figsize=(width, height))
    # force vmin and vmax for consistency between plots:
    heatmap(metrics, vmin=0, vmax=1, ax=ax, cmap=palette)
    # adjust xticks according to bin size and alignment:
    xticks = [
        int(tick.get_text()) for tick in ax.get_xticklabels()
    ]
    if align == "right":
        xticks = [
            tick * -1 for tick in reversed(xticks)
        ]
    if xtick_density != 1:
        each = int(len(xticks)*xtick_density)
        from sys import stderr
        print(each, file=stderr)
        xticks = [
            xtick if i%each==0 else ""
            for i, xtick in enumerate(xticks)
        ]
    ax.set_xticklabels(
        [tick * bin_size for tick in xticks],
        rotation=60
    )
    # drop yticks (useful for long names / big sets):
    if hide_names:
        ax.set(yticks=[])
    ax.set(title=title, xlabel="position (bp)")
    figure.savefig(png)


def main(args):
    """Dispatch data to subroutines"""
    metrics = load_metrics(
        args.txt, bin_size=args.bin_size, align=args.align
    )
    if args.title:
        title = args.title
    else:
        title = args.txt.split("/")[-1]
    plot_metrics(
        metrics, args.figsize, args.palette, args.hide_names,
        bin_size=args.bin_size, align=args.align,
        xtick_density=args.xtick_density, title=title, png=stdout.buffer
    )
    return 0


if __name__ == "__main__":
    parser = ArgumentParser(usage=USAGE)
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    if args.align not in {"left", "right"}:
        raise ValueError("--align must be 'left' or 'right'")
    if not search(r'^\d+x\d+$', args.figsize):
        raise ValueError("--figsize must be in format WxH: e.g., 5x7")
    else:
        returncode = main(args)
        exit(returncode)

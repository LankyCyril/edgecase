#!/usr/bin/env python3
from sys import stdout
from re import search
from argparse import ArgumentParser
from numpy import linspace, array, mean, nan, concatenate
from pandas import Series, concat
from matplotlib.pyplot import switch_backend, subplots
from seaborn import heatmap

USAGE = "python3 {} [options] txt > png".format(__file__)

ARG_RULES = {
    ("txt",): {
        "help": "name of input (each line is 'name\\tvalue\\tvalue\\tvalue...')"
    },
    ("-b", "--bin-size"): {
        "help": "size of each bin in visualization (1)",
        "default": 1, "type": int, "metavar": "B"
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
    }
}


def binned(A, bins, func=mean):
    """Return array data compressed into bins (smoothed by func)"""
    if bins > 1:
        coords = linspace(0, len(A), bins+1).astype(int)
        return array([
            func(A[start:end])
            for start, end in zip(coords, coords[1:])
        ])
    else:
        return A


def load_metrics(txt, bin_size=120, align="left"):
    """Load metrics from a text file, bin and convert into dataframe"""
    read_metrics = {}
    # load and bin all metrics first:
    with open(txt, mode="rt") as handle:
        for line in handle:
            name, *metrics = line.strip().split()
            read_metrics[name] = binned(
                array(metrics, dtype="float32"),
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


def plot_metrics(metrics, figsize="16x9", palette="viridis", hide_names=False, title="", xlabel="", png=stdout.buffer):
    """Plot binned metrics as a heatmap"""
    switch_backend("Agg")
    width, height = tuple(map(int, figsize.split("x")))
    figure, ax = subplots(figsize=(width, height))
    heatmap(metrics, vmin=0, vmax=1, ax=ax, cmap=palette)
    ax.set(xlabel=xlabel, title=title)
    if hide_names:
        ax.set(yticks=[])
    figure.savefig(png)


def main(args):
    """Dispatch data to subroutines"""
    metrics = load_metrics(
        args.txt, bin_size=args.bin_size, align=args.align
    )
    if args.bin_size > 1:
        xlabel = "position (x{})".format(args.bin_size)
    else:
        xlabel = "position"
    plot_metrics(
        metrics, args.figsize, args.palette, args.hide_names,
        title=args.txt.split("/")[-1], xlabel=xlabel, png=stdout.buffer
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

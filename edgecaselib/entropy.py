from sys import stdout
from edgecaselib.formats import load_kmerscan
from pandas import DataFrame, concat
from scipy.stats import entropy
from numpy import log, argsort, cumsum, interp
from edgecaselib.util import progressbar
from itertools import chain


__doc__ = """edgeCase entropy: calculation of motif entropy among reads

Usage: {0} entropy [-b integer] [-f flagspec]... [-F flagspec]... [-q integer]
       {1}         [-z] <dat>...

Output:
    TSV file of entropy values and read coverage per bin,
    with coverage-weighted quantiles of entropy in a comment on first line

Positional arguments:
    <dat>                         name of input kmerscanner file(s)

Options:
    -z, --gzipped                 input is gzipped (must specify if any of -qfF present)
    -b, --bin-size [integer]      size of each bin in bp (overrides bin size in <dat>)

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda bin_size: None if bin_size is None else int(bin_size),
]


def calculate_entropies(bdf):
    """Calculate entropies per binned density dataframe"""
    per_read_modes = (
        bdf.groupby("name")
        .apply(lambda block: block.set_index("motif").iloc[:,8:].idxmax(axis=0))
        .dropna(how="all", axis=1)
    )
    N = len(per_read_modes.melt().value.dropna().unique())
    return DataFrame({
        "#entropy": (
            per_read_modes.apply(lambda c: entropy(c.value_counts())) / log(N)
        ),
        "coverage": (~per_read_modes.isnull()).sum(axis=0),
    })


def weighted_quantile(points, weights, q):
    """Calculate quantile `q` of `points`, weighted by `weights`"""
    indsort = argsort(points.values)
    spoints, sweights = points.values[indsort], weights.values[indsort]
    sn = cumsum(sweights)
    pn = (sn - sweights / 2) / sn[-1]
    return interp(q, pn, spoints)


def main(dat, gzipped, flags, flag_filter, min_quality, bin_size, file=stdout, **kwargs):
    """Dispatch data to subroutines"""
    samfilters = [flags, flag_filter, min_quality]
    kmerscans = [load_kmerscan(fn, gzipped, samfilters, bin_size) for fn in dat]
    entropies = concat(
        calculate_entropies(bdf) for bdf in progressbar(
            chain(*(ks.values() for ks in kmerscans)),
            desc="Calculating entropies", unit="arm",
            total=sum(len(ks) for ks in kmerscans),
        )
    )
    quantiles = {
        q: weighted_quantile(
            entropies["#entropy"], entropies["coverage"]-1, q/100,
        )
        for q in progressbar(range(5, 101, 5), desc="Calculating quantiles")
    }
    print("#"+",".join(f"q{k}={v}" for k, v in quantiles.items()), file=file)
    entropies.to_csv(file, sep="\t", index=False)

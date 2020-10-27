from sys import stdout
from edgecaselib.formats import load_index, load_kmerscan, interpret_flags
from pandas import DataFrame, concat
from scipy.stats import entropy
from numpy import log
from edgecaselib.util import progressbar


__doc__ = """edgeCase entropy: calculation of motif entropy among reads

Usage: {0} entropy -x filename [-t targetspec] [-b integer] [-z]
       {1}         [-f flagspec] [-F flagspec] [-q integer] <dat>

Output:
    TSV file of entropy values and read coverage per bin,
    with coverage-weighted quantiles of entropy in a comment on first line

Positional arguments:
    <dat>                         name of input kmerscanner file

Required options:
    -x, --index [filename]        location of the reference .ecx index

Options:
    -z, --gzipped                 input is gzipped (must specify if any of -qfF present)
    -b, --bin-size [integer]      size of each bin in bp (overrides bin size in <dat>)
    -t, --target [targetspec]     an ECX flag past which to evaluate motifs [default: tract_anchor]

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda bin_size: None if bin_size is None else int(bin_size),
]


def calculate_entropies(bdf, chrom, ecx, integer_target):
    """Calculate entropies per binned density dataframe"""
    flags = bdf["flag"].drop_duplicates()
    is_q = set(flag & 0x8000 == 0x8000 for flag in flags)
    if len(is_q) != 1:
        raise ValueError(f"p- and q-arm reads mixed on {chrom} in DAT")
    else:
        prime = 3 if is_q.pop() else 5
    anchors = ecx.loc[
        ((ecx["rname"]==chrom) & (ecx["prime"]==prime) &
        (ecx["flag"]==integer_target)), "pos"
    ]
    if len(anchors) == 0:
        raise ValueError(f"No relevant {prime}-prime info on {chrom} in ECX")
    elif len(anchors) > 1:
        raise ValueError(f"Conflicting {prime}-prime info on {chrom} in ECX")
    else:
        anchor = anchors.iloc[0]
    if prime == 3:
        positions = [c for c in bdf.columns[9:] if c >= anchor]
    else:
        positions = [c for c in bdf.columns[9:] if c <= anchor]
    bdf_sliced = bdf[list(bdf.columns[:9])+positions]
    per_read_modes = (
        bdf_sliced.groupby("name")
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


def main(dat, index, gzipped, flags, flag_filter, min_quality, bin_size, target, file=stdout, **kwargs):
    """Dispatch data to subroutines"""
    samfilters = [flags, flag_filter, min_quality]
    ecx = load_index(index)
    integer_target = interpret_flags(target)
    densities = load_kmerscan(dat, gzipped, samfilters, bin_size)
    entropies = concat(
        calculate_entropies(bdf, chrom, ecx, integer_target)
        for chrom, bdf in progressbar(
            densities.items(), total=len(densities),
            desc="Calculating entropies", unit="arm",
        )
    )
    entropies.to_csv(file, sep="\t", index=False)

from sys import stderr
from numpy import linspace, array, mean, nan, concatenate, fromstring, full, vstack
from pandas import read_csv, merge, concat, DataFrame
from gzip import open as gzopen
from tempfile import TemporaryDirectory
from os import path
from tqdm import tqdm
from functools import reduce
from operator import __or__


ECX_FLAGS = {
    "ucsc_mask_anchor": 0x1000,
    "fork": 0x2000,
    "tract_anchor": 0x4000,
    "is_q": 0x8000
}


def passes_filter(entry_flag, entry_mapq, flags, flag_filter, min_quality):
    """Check if entry flags pass filters"""
    passes_quality = (entry_mapq is None) or (entry_mapq >= min_quality)
    if not passes_quality:
        return False
    else:
        return (entry_flag is None) or (
            (entry_flag & flags == flags) and
            (entry_flag & flag_filter == 0)
        )


def filter_and_read_csv(dat, gzipped, flags, flag_filter, min_quality):
    """If filters supplied, subset CSV first, then read with pandas"""
    number_retained = 0
    if gzipped:
        opener = gzopen
    else:
        opener = open
    with opener(dat, mode="rt") as dat_handle:
        with TemporaryDirectory() as tempdir:
            datflt_name = path.join(tempdir, "dat.gz")
            with gzopen(datflt_name, mode="wt") as datflt:
                decorated_line_iterator = tqdm(
                    dat_handle, desc="Filtering", unit=" records"
                )
                for line in decorated_line_iterator:
                    if line[0] == "#":
                        print(line, end="", file=datflt)
                    else:
                        fields = line.split("\t")
                        line_passes_filter = passes_filter(
                            int(fields[1]), int(fields[4]),
                            flags, flag_filter, min_quality
                        )
                        if line_passes_filter:
                            number_retained += 1
                            print(line, end="", file=datflt)
                print("Kept {} records".format(number_retained), file=stderr)
            return read_csv(datflt_name, sep="\t", escapechar="#")


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


def load_kmerscan(dat, gzipped, flags, flag_filter, min_quality, bin_size, no_align, each_once=True):
    """Load densities from dat file, split into dataframes per chromosome"""
    if (flags == 0) and (flag_filter == 0) and (min_quality == 0):
        raw_densities = read_csv(dat, sep="\t", escapechar="#")
    else:
        raw_densities = filter_and_read_csv(
            dat, gzipped, flags, flag_filter, min_quality
        )
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


def interpret_flags(flags):
    """If flags are not a decimal number, assume strings and convert to number"""
    if isinstance(flags, int) or flags.isdigit():
        return int(flags)
    elif not isinstance(flags, str):
        raise ValueError("Unknown flags: {}".format(repr(flags)))
    elif flags[:2] == "0x":
        return int(flags, 16)
    elif flags[:2] == "0b":
        return int(flags, 2)
    elif "|" in flags:
        flag_set = set(map(interpret_flags, flags.split("|")))
        return reduce(__or__, flag_set | {0})
    elif flags in ECX_FLAGS:
        return ECX_FLAGS[flags]
    else:
        raise ValueError("Unknown flags: {}".format(repr(flags)))

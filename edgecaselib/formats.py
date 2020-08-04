from sys import stderr
from collections import OrderedDict
from contextlib import contextmanager, ExitStack
from itertools import chain
from numpy import linspace, array, mean, concatenate, fromstring, full, vstack
from numpy import nan, unique
from pandas import read_csv, merge, concat, DataFrame
from gzip import open as gzopen
from tempfile import TemporaryDirectory
from os import path
from edgecaselib.util import progressbar
from functools import reduce
from operator import __or__


ALL_SAM_FLAGS = [
    "paired", "mapped_proper_pair", "unmapped", "mate_unmapped", "rev",
    "mate_rev", "1stmate", "2ndmate", "secondary", "qcfail", "pcrdup", "supp",
    "ucsc_mask_anchor", "fork", "tract_anchor", "is_q"
]

FLAG_COLORS = {0x1000: "gray", 0x2000: "blueviolet", 0x4000: "red"}

DEFAULT_MOTIF_COLORS = [
    "#117733", "#88CCEE", "#AA4499", "#DDCC77", "#332288", "#882255",
    "#44AA99", "#CC6677", "#EEEEEE",
]

BGCOLOR = "#BBBBCA"

PAPER_PALETTE = OrderedDict([
    ("TTAGGG", "#117733"), ("TTGGGG", "#AA4499"), ("TTAGGGG", "#332288"),
    ("TGAGGG", "#DDCC77"), ("TCAGGG", "#44AA99"), ("TTAGGGTTAGGGG", "#EEEEEE"),
    ("CGCGG", "#88CCEE"),
])

PAPER_PALETTE_RC = OrderedDict([
    ("CCCTAA", "#117733"), ("CCCCAA", "#AA4499"), ("CCCCTAA", "#332288"),
    ("CCCTCA", "#DDCC77"), ("CCCTGA", "#44AA99"), ("CCCCTAACCCTAA", "#EEEEEE"),
    ("CCGCG", "#88CCEE"),
])


class EmptyKmerscanError(ValueError):
    """Raised when supplied kmerscanner file is empty"""
    pass


def explain_sam_flags(flag, sep="|"):
    """Convert an integer flag into string"""
    return sep.join(ALL_SAM_FLAGS[i] for i in range(16) if flag & 2**i != 0)


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
    elif flags in ALL_SAM_FLAGS:
        return 2**ALL_SAM_FLAGS.index(flags)
    else:
        raise ValueError("Unknown flags: {}".format(repr(flags)))


def entry_filters_ok(entry_flag, entry_mapq, integer_samfilters):
    """Check if entry flags and mapq pass filters"""
    if (entry_mapq is None) or (entry_flag is None):
        return True
    if entry_mapq < integer_samfilters[-1]:
        return False
    else:
        return (
            # -f, flags (all must be present):
            (entry_flag & integer_samfilters[0] == integer_samfilters[0]) and
            # -F, flag_filter (all must be absent):
            (entry_flag & integer_samfilters[1] == 0)
        )


def filter_and_read_tsv(dat, gzipped, samfilters):
    """If filters supplied, subset DAT first, then read with pandas"""
    number_retained = 0
    if gzipped:
        opener = gzopen
    else:
        opener = open
    integer_samfilters = list(map(interpret_flags, samfilters))
    with opener(dat, mode="rt") as dat_handle:
        with TemporaryDirectory() as tempdir:
            datflt_name = path.join(tempdir, "dat.gz")
            with gzopen(datflt_name, mode="wt") as datflt:
                decorated_line_iterator = progressbar(
                    dat_handle, desc="Filtering", unit=" lines",
                )
                for line in decorated_line_iterator:
                    if line[0] == "#":
                        print(line, end="", file=datflt)
                    else:
                        fields = line.split("\t")
                        line_passes_filter = entry_filters_ok(
                            int(fields[1]), int(fields[4]), integer_samfilters,
                        )
                        if line_passes_filter:
                            number_retained += 1
                            print(line, end="", file=datflt)
                print("Kept {} records".format(number_retained), file=stderr)
            print("Loading DAT...", file=stderr, flush=True)
            return read_csv(datflt_name, sep="\t", escapechar="#")


def binned(A, bins, func=mean):
    """Return array data compressed into bins (smoothed by func)"""
    coords = linspace(0, len(A), bins+1).astype(int)
    return array([
        func(A[start:end])
        for start, end in zip(coords, coords[1:])
    ])


def get_binned_density_dataframe(raw_densities, chrom, bin_size, no_align=False):
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
            aligned_density, bins=aligned_density.shape[0]//bin_size,
        ))
    # pad densities on the right so they are all the same length:
    max_density_length = max(bd.shape[0] for bd in density_arrays)
    for i, binned_density in enumerate(density_arrays):
        padder = full(max_density_length - binned_density.shape[0], nan)
        density_arrays[i] = concatenate([binned_density, padder])
    naked_binned_density_dataframe = DataFrame(
        data=vstack(density_arrays),
        columns=[leftmost_pos+j*bin_size for j in range(max_density_length)],
    )
    binned_density_dataframe = concat(
        [densities_subset.iloc[:,:-1], naked_binned_density_dataframe], axis=1,
    )
    return binned_density_dataframe.sort_values(
        by=["mapq", "name", "motif"], ascending=[False, True, True],
    )


def are_motifs_consistent(raw_densities):
    checker = raw_densities[["chrom", "motif"]].drop_duplicates()
    if len(unique(checker["chrom"].value_counts().values)) != 1:
        return False
    elif len(unique(checker["motif"].value_counts().values)) != 1:
        return False
    else:
        return True


def load_kmerscan(dat, gzipped, samfilters, bin_size, no_align=False, each_once=True):
    """Load densities from dat file, split into dataframes per chromosome"""
    if not any(samfilters): # all zero / None
        print("Loading DAT...", file=stderr, flush=True)
        raw_densities = read_csv(dat, sep="\t", escapechar="#")
    else:
        raw_densities = filter_and_read_tsv(dat, gzipped, samfilters)
    if len(raw_densities) == 0:
        raise EmptyKmerscanError
    if not are_motifs_consistent(raw_densities):
        raise NotImplementedError(
            "Inconsistent number of motifs in DAT; plotting of reads " +
            "identified de novo with kmerscanner is not implemented"
        )
    if each_once:
        raw_densities["length"] = raw_densities["density"].apply(
            lambda d: d.count(",")+1
        )
        groups = raw_densities[["name", "motif", "length"]].groupby(
            ["name", "motif"], as_index=False,
        ).max()
        raw_densities = merge(groups, raw_densities).drop(columns="length")
    if no_align:
        raw_densities["chrom"] = "None"
    chromosome_iterator = progressbar(
        raw_densities["chrom"].drop_duplicates(), desc="Interpreting data",
        unit="chromosome",
    )
    return {
        chrom: get_binned_density_dataframe(
            raw_densities, chrom, bin_size, no_align,
        )
        for chrom in chromosome_iterator
    }


def load_index(index_filename, as_filter_dict=False):
    """Load ECX index; convert to simpler dictionary for filtering if requested"""
    if not path.isfile(index_filename):
        raise FileNotFoundError(index_filename)
    else:
        ecx = read_csv(
            index_filename, sep="\t", skiprows=1,
            escapechar="#", na_values="-",
        )
        ecx = ecx[ecx["blacklist"].isnull()]
        if as_filter_dict:
            slim_ecx = ecx[["rname", "pos", "flag", "prime"]]
            rnames = set(slim_ecx["rname"])
            return {
                rname: {
                    5: slim_ecx[
                        (slim_ecx["rname"]==rname) & (slim_ecx["prime"]==5)
                    ].drop(columns=["rname", "prime"]),
                    3: slim_ecx[
                        (slim_ecx["rname"]==rname) & (slim_ecx["prime"]==3)
                    ].drop(columns=["rname", "prime"])
                }
                for rname in rnames
            }
        else:
            return ecx


@contextmanager
def ReadFileChain(filenames, manager):
    """Chain records from all filenames in list bams, replicating behavior of pysam context managers"""
    with ExitStack() as stack:
        yield chain(*(
            stack.enter_context(manager(filename))
            for filename in filenames
        ))


def filter_bam(alignment, samfilters, desc=None):
    """Wrap alignment iterator with a flag and quality filter"""
    integer_samfilters = list(map(interpret_flags, samfilters))
    filtered_iterator = (
        entry for entry in alignment
        if entry_filters_ok(entry.flag, entry.mapq, integer_samfilters)
    )
    if desc is None:
        return filtered_iterator
    else:
        return progressbar(filtered_iterator, desc=desc, unit="read")

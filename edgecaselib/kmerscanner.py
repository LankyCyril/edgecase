from sys import stdout
from numpy import zeros, array, cumsum, nan
from multiprocessing import Pool
from edgecaselib.util import get_circular_pattern
from edgecaselib.formats import filter_bam
from edgecaselib.tailchopper import get_cigar_clip_length
from pysam import AlignmentFile, FastxFile
from types import SimpleNamespace
from functools import partial
from edgecaselib.util import progressbar
from collections import OrderedDict
from pandas import read_csv


__doc__ = """edgeCase kmerscanner: calculation of motif densities

Usage: {0} kmerscanner [-j integer] --motif-file filename
       {1}             [-b integer] [-n integer]
       {1}             [-f flagspec]... [-F flagspec]... [-q integer]
       {1}             [--fmt string] <sequencefile>

Output:
    DAT file with calculated motif densities along rolling windows

Positional arguments:
    <sequencefile>                name of input BAM/SAM/FASTA/FASTQ file

Required options:
    --motif-file [filename]       file with repeated motif sequences (output of `repeatfinder`)

Options:
    --fmt sam|fastx               format of input file [default: sam]
    -b, --bin-size [integer]      size of the rolling window [default: 10]
    -n, --num-reads [integer]     expected number of reads in input (for progress display)
    -j, --jobs [integer]          number of jobs to run in parallel [default: 1]

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda bin_size: int(bin_size),
    lambda num_reads: None if (num_reads is None) else int(num_reads),
    lambda jobs: int(jobs),
    lambda min_quality: None if (min_quality is None) else int(min_quality),
]


DAT_HEADER = [
    "#name", "flag", "chrom", "pos", "mapq", "motif", "score",
    "clip_5prime", "clip_3prime",
]


def calculate_density(entry, pattern, bin_size, positions_accounted_for):
    """Calculate density of pattern hits in a rolling window along given read"""
    if entry.query_sequence is not None: # calculations will make sense
        read_length = len(entry.query_sequence)
        canvas = zeros(read_length, dtype=bool)
        pattern_positions = array([
            match.start() for match
            in pattern.finditer(entry.query_sequence, overlapped=True)
            if match.start() not in positions_accounted_for
        ])
        if len(pattern_positions):
            canvas[pattern_positions] = True
        if read_length <= bin_size: # use one window:
            density_array = (canvas.sum(axis=0) / read_length).reshape(1)
        else: # use rolling window:
            roller = cumsum(canvas, axis=0)
            roller[bin_size:] = roller[bin_size:] - roller[:-bin_size]
            density_array = roller[bin_size-1:] / bin_size
        return True, density_array, pattern_positions
    else: # effectively skip
        return False, zeros(1), array([])


def calculate_density_of_patterns(entry, motif_patterns, bin_size):
    """Calculate density of hits of each pattern in a rolling window along given read"""
    entry_set = []
    positions_accounted_for = set()
    for motif, pattern in motif_patterns.items():
        passes_filter, density_array, pattern_positions = calculate_density(
            entry, pattern, bin_size, positions_accounted_for,
        )
        if passes_filter:
            entry_set.append([entry, motif, density_array])
        if len(pattern_positions):
            positions_accounted_for |= set(pattern_positions)
    return entry_set


def pattern_scanner(entry_iterator, fmt, samfilters, motif_patterns, bin_size, num_reads, jobs):
    """Calculate density of pattern hits in a rolling window along each read"""
    if fmt == "sam":
        filtered_iterator = filter_bam(entry_iterator, samfilters)
    else:
        filtered_iterator = entry_iterator
    simple_entry_iterator = (
        SimpleNamespace(
            query_name=getattr(
                entry, "query_name", getattr(entry, "name", None)
            ),
            flag=getattr(entry, "flag", None),
            reference_name=getattr(entry, "reference_name", None),
            reference_start=getattr(entry, "reference_start", None),
            mapping_quality=getattr(entry, "mapping_quality", None),
            query_sequence=getattr(
                entry, "query_sequence", getattr(entry, "sequence", None)
            ),
            cigarstring=getattr(entry, "cigarstring", ""),
        )
        for entry in filtered_iterator
    )
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density_of_patterns, motif_patterns=motif_patterns,
            bin_size=bin_size,
        )
        # lazy multiprocess evaluation:
        read_density_iterator = pool.imap_unordered(
            density_calculator, simple_entry_iterator,
        )
        # iterate pairs (entry.query_name, density_array), same as calculate_density_of_patterns():
        desc = "Calculating density"
        yield from progressbar(
            read_density_iterator, desc=desc, unit="read", total=num_reads,
        )


def interpret_arguments(fmt, motif_file):
    """Parse and check arguments"""
    if fmt == "sam":
        manager = AlignmentFile
    elif fmt == "fastx":
        manager = FastxFile
    else:
        raise ValueError("--fmt can only be 'sam' or 'fastx'")
    motif_data = read_csv(motif_file, sep="\t", escapechar="#")
    if "motif" not in motif_data.columns:
        if "monomer" in motif_data.columns:
            motif_data = motif_data.rename(columns={"monomer": "motif"})
        else:
            raise KeyError("No motif column found in motif file")
    if "length" not in motif_data.columns:
        motif_data["length"] = motif_data["motif"].apply(lambda m: len(m))
    motif_data = motif_data.sort_values(by="length", ascending=False)
    motif_patterns = OrderedDict([
        [motif, get_circular_pattern(motif)] for motif in motif_data["motif"]
    ])
    if "score" in motif_data.columns:
        scores = dict(zip(
            motif_data["motif"], motif_data["score"]
        ))
    else:
        scores = {m: nan for m in motif_data["motif"]}
    return manager, motif_patterns, scores


def main(sequencefile, fmt, flags, flag_filter, min_quality, motif_file, bin_size, num_reads, jobs=1, file=stdout, **kwargs):
    # parse and check arguments:
    manager, motif_patterns, scores = interpret_arguments(fmt, motif_file)
    print(*DAT_HEADER, f"b={bin_size}", sep="\t", file=file)
    # scan fastq for target motif queries, parallelizing on reads:
    with manager(sequencefile) as entry_iterator:
        scanner = pattern_scanner(
            entry_iterator, motif_patterns=motif_patterns,
            samfilters=[flags, flag_filter, min_quality],
            fmt=fmt, bin_size=bin_size, num_reads=num_reads, jobs=jobs,
        )
        # output densities of reads that pass filter:
        for entry_set in scanner:
            for entry, motif, density_array in entry_set:
                if entry: # non-null result, entry passed filters
                    meta_fields = [
                        entry.query_name, entry.flag, entry.reference_name,
                        entry.reference_start, entry.mapping_quality,
                        motif, scores[motif],
                        get_cigar_clip_length(entry, 5),
                        get_cigar_clip_length(entry, 3),
                    ]
                    print(*meta_fields, sep="\t", end="\t", file=file)
                    print(*density_array, sep=",", file=file)

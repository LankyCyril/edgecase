from sys import stdout
from os import path
from edgecaselib.formats import load_index, filter_bam
from pysam import AlignmentFile
from functools import reduce
from operator import __or__
from copy import deepcopy
from edgecaselib.util import progressbar
from itertools import chain
from numpy import isnan, inf
from pandas import DataFrame


__doc__ = """edgeCase tailpuller: selection of candidate telomeric long reads

Usage: {0} tailpuller -x filename [-f flagspec] [-F flagspec] [-q integer]
       {1}            [-M integer] [-m integer] [-a] <bam>

Output:
    SAM-formatted file with reads overhanging anchors defined in index

Positional arguments:
    <bam>                             name of input BAM/SAM file; must have a .bai index

Required options:
    -x, --index [filename]            location of the reference .ecx index

Options:
    -M, --max-read-length [integer]   maximum read length to consider when selecting lookup regions
    -m, --min-overlap [integer]       minimum overlap of subtelomere to consider read [default: 0]
    -a, --allow-ambiguous             include reads that map to multiple terminal positions / other arms

Input filtering options:
    -f, --flags [flagspec]            process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]      process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]       process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda min_quality:
        None if (min_quality is None) else int(min_quality),
    lambda min_overlap:
        0 if (min_overlap is None) else int(min_overlap),
    lambda max_read_length:
        inf if (max_read_length is None) else int(max_read_length),
]

__docopt_tests__ = {
    lambda bam:
        path.isfile(bam + ".bai"): "BAM index (.bai) not found",
    lambda max_read_length:
        max_read_length > 0: "--max-read-length below 0",
}


def get_terminal_pos(entry, cigarpos):
    """Calculate the position of clipped start/end of read relative to the reference"""
    if not entry.cigartuples: # no CIGAR available
        return None
    # measure the clipped stretch on the left (0) or right (-1):
    cigartype, clip = entry.cigartuples[cigarpos]
    if (cigartype != 4) and (cigartype != 5): # not a soft/hard clip
        clip = 0
    # determine location of start/end of the read relative to reference:
    if cigarpos == 0: # start
        return entry.reference_start - clip
    elif cigarpos == -1: # end
        return entry.reference_end + clip
    else:
        raise ValueError("get_terminal_pos(): cigarpos can only be 0 or -1")


def updated_entry(entry, flags, is_q=False):
    """Add ECX flags to entry"""
    new_entry = deepcopy(entry)
    new_entry.flag |= reduce(__or__, flags)
    if is_q:
        new_entry.flag |= 0x8000
    return new_entry


def filter_entries(bam_data, ecxfd, min_overlap, samfilters):
    """Only pass reads extending past regions specified in the ECX"""
    for entry in filter_bam(bam_data, samfilters):
        if entry.reference_name in ecxfd:
            if (entry.reference_length or 0) >= min_overlap:
                # find positions of start and end of read relative to reference:
                p_pos = get_terminal_pos(entry, cigarpos=0)
                q_pos = get_terminal_pos(entry, cigarpos=-1)
                # collect ECX flags where anchor is to right of read start:
                ecx_t_p5 = ecxfd[entry.reference_name][5]
                p_flags = set(ecx_t_p5.loc[ecx_t_p5["pos"]>=p_pos, "flag"])
                if p_flags:
                    yield updated_entry(entry, p_flags)
                # collect ECX flags where anchor is to left of read end
                ecx_t_p3 = ecxfd[entry.reference_name][3]
                q_flags = set(ecx_t_p3.loc[ecx_t_p3["pos"]<q_pos, "flag"])
                if q_flags:
                    yield updated_entry(entry, q_flags, is_q=True)


def get_bam_chunk(bam_data, chrom, ecxfd, reflens, max_read_length):
    """Subset bam_data to a region where reads of interest can occur"""
    if chrom not in reflens:
        return []
    elif (max_read_length is None) or (max_read_length == inf):
        return bam_data.fetch(chrom, None, None)
    else:
        p_innermost_pos = ecxfd[chrom][5]["pos"].max() + max_read_length
        if p_innermost_pos > reflens[chrom] - 1:
            p_innermost_pos = reflens[chrom] - 1
        elif p_innermost_pos < 0:
            p_innermost_pos = 0
        q_innermost_pos = ecxfd[chrom][3]["pos"].min() - max_read_length
        if q_innermost_pos > reflens[chrom] - 1:
            q_innermost_pos = reflens[chrom] - 1
        elif q_innermost_pos < 0:
            q_innermost_pos = 0
        if isnan(p_innermost_pos) and isnan(q_innermost_pos):
            return None
        elif (not isnan(p_innermost_pos)) and isnan(q_innermost_pos):
            return bam_data.fetch(chrom, 0, p_innermost_pos)
        elif isnan(p_innermost_pos) and (not isnan(q_innermost_pos)):
            return bam_data.fetch(chrom, q_innermost_pos, None)
        elif p_innermost_pos >= q_innermost_pos:
            return bam_data.fetch(chrom, None, None)
        else:
            return chain(
                bam_data.fetch(chrom, 0, p_innermost_pos),
                bam_data.fetch(chrom, q_innermost_pos, None)
            )


def parse_bam_with_ambiguity(bam, ecxfd, max_read_length, min_overlap, samfilters):
    """Parse BAM file, select overhanging reads, possibly mapping to more than one arm"""
    with AlignmentFile(bam) as bam_data:
        reflens = dict(zip(bam_data.references, bam_data.lengths))
        bam_header_string = str(bam_data.header).rstrip("\n")
        decorated_bam_iterator = progressbar(
            ecxfd, total=len(ecxfd), desc="Pulling", unit="chromosome",
        )
        candidates = []
        for chrom in decorated_bam_iterator:
            bam_chunk = get_bam_chunk(
                bam_data, chrom, ecxfd, reflens, max_read_length,
            )
            candidates.extend(
                filter_entries(bam_chunk, ecxfd, min_overlap, samfilters),
            )
    return bam_header_string, candidates


def get_unambiguous_candidates(candidates):
    """Subset candidate entries to those that map unambiguously"""
    dispatcher = DataFrame(
        columns=["qname", "flag"],
        data=[[entry.qname, entry.flag] for entry in candidates],
    )
    primaries = set(dispatcher.loc[
        dispatcher["flag"] & 3844 == 0, "qname",
    ])
    primaries_with_duplicates = dispatcher.loc[
        dispatcher["qname"].isin(primaries), "qname",
    ]
    name_counts = primaries_with_duplicates.value_counts()
    unique_candidate_names = set(name_counts[name_counts==1].index)
    for entry in candidates:
        if entry.qname in unique_candidate_names:
            yield entry


def main(bam, index, flags, flag_filter, min_quality, max_read_length, min_overlap, allow_ambiguous, file=stdout, **kwargs):
    # dispatch data to subroutines:
    ecxfd = load_index(index, as_filter_dict=True)
    samfilters = [flags, flag_filter, min_quality]
    bam_header_string, candidates = parse_bam_with_ambiguity(
        bam, ecxfd, max_read_length, min_overlap, samfilters,
    )
    print(bam_header_string, file=file)
    if allow_ambiguous:
        for entry in candidates:
            print(entry.to_string(), file=file)
    else:
        for entry in get_unambiguous_candidates(candidates):
            print(entry.to_string(), file=file)

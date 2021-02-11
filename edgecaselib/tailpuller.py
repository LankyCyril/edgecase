from sys import stdout, stderr
from os import path
from edgecaselib.formats import interpret_flags, load_index, filter_bam
from pysam import AlignmentFile
from functools import reduce
from operator import __or__
from copy import deepcopy
from edgecaselib.util import progressbar
from itertools import chain
from numpy import isnan, inf
from pandas import DataFrame, merge


__doc__ = """edgeCase tailpuller: selection of candidate telomeric long reads

Usage: {0} tailpuller -x filename [-t targetspec]...
       {1}            [-M integer] [--min-map-overlap integer]
       {1}            [-m integer] [--min-telomere-overlap integer]
       {1}            [-f flagspec]... [-F flagspec]... [-q integer] <bam>

Output:
    SAM-formatted file with reads overhanging anchors defined in index

Positional arguments:
    <bam>                                    name of input BAM/SAM file; must have a .bai index

Required options:
    -x, --index [filename]                   location of the reference .ecx index

Options:
    -t, --target [targetspec]                target reads overlapping these features (ECX flags) [default: tract_anchor]
    -M, --max-read-length [integer]          maximum read length to consider when selecting lookup regions *
    --min-map-overlap [integer]              minimum overlap of reference to consider read as mapped [default: 1] **
    -m, --min-subtelomere-overlap [integer]  minimum overlap of subtelomere to consider read as candidate [default: 1] ***
    --min-telomere-overlap [integer]         minimum overlap of telomere to consider read as candidate [default: 1] ***

Input filtering options:
    -f, --flags [flagspec]                   process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]             process only entries with none of these sam flags present [default: 0] ****
    -q, --min-quality [integer]              process only entries with this MAPQ or higher [default: 0] *****

Notes:
    * Suggested value of --max-read-length for PacBio HiFi: 30000;
      if not specified, will assume +infinity (will be slow).
   ** Suggested value of --min-map-overlap for PacBio HiFi: 500;
  *** Suggested value of --min-(sub)telomere-overlap for PacBio HiFi: 3000;
 **** It is recommended to include secondary and supplementary reads (i.e.,
      leave the -F flag as default [0]), because:
        **** edgeCase determines unambiguously mapped reads on its own; aligners
             assign the 'supplementary' flag to multi-mapping reads arbitrarily,
             and removing such supplementary reads upstream may lead to loss of
             information in telomeric regions;
        **** edgeCase will discard chimeric reads in terminal regions if
             information about supplementary alignments is present.
***** Depending on the aligner used, MAPQ of secondary reads may have been set
      to zero regardless of real mapping quality; use this filtering option with
      caution.
"""

__docopt_converters__ = [
    lambda min_quality:
        None if (min_quality is None) else int(min_quality),
    lambda max_read_length:
        inf if (max_read_length is None) else int(max_read_length),
    lambda min_map_overlap:
        1 if (min_map_overlap is None) else int(min_map_overlap),
    lambda min_subtelomere_overlap:
        1 if (min_subtelomere_overlap is None) else int(min_subtelomere_overlap),
    lambda min_telomere_overlap:
        None if (min_telomere_overlap is None) else int(min_telomere_overlap),
]

__docopt_tests__ = {
    lambda bam:
        path.isfile(bam + ".bai"): "BAM index (.bai) not found",
    lambda max_read_length:
        max_read_length > 0: "--max-read-length below 0",
    lambda target:
        set(target) <= {"mask_anchor", "fork", "tract_anchor"}:
            "unknown value of --target",
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


def filter_entries(bam_data, ecxfd, targets, samfilters, min_map_overlap):
    """Only pass reads extending past regions specified in the ECX"""
    for entry in filter_bam(bam_data, samfilters):
        if entry.reference_name in ecxfd:
            if entry.reference_length >= min_map_overlap:
                # find positions of start and end of read relative to reference:
                p_pos = get_terminal_pos(entry, cigarpos=0)
                q_pos = get_terminal_pos(entry, cigarpos=-1)
                # collect ECX flags where anchor is to right of read start:
                ecx_t_p5 = ecxfd[entry.reference_name][5]
                p_flags = set(ecx_t_p5.loc[ecx_t_p5["pos"]>=p_pos, "flag"])
                if targets & p_flags:
                    yield updated_entry(entry, p_flags)
                # collect ECX flags where anchor is to left of read end
                ecx_t_p3 = ecxfd[entry.reference_name][3]
                q_flags = set(ecx_t_p3.loc[ecx_t_p3["pos"]<q_pos, "flag"])
                if targets & q_flags:
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


def parse_bam_with_ambiguity(bam, ecxfd, max_read_length, min_map_overlap, targets, samfilters):
    """Parse BAM file, select overhanging reads, possibly mapping to more than one arm"""
    with AlignmentFile(bam) as bam_data:
        reflens = dict(zip(bam_data.references, bam_data.lengths))
        bam_header_string = str(bam_data.header).rstrip("\n")
        decorated_bam_iterator = progressbar(
            ecxfd, total=len(ecxfd), desc="Pulling", unit="chromosome",
        )
        entries = []
        for chrom in decorated_bam_iterator:
            bam_chunk = get_bam_chunk(
                bam_data, chrom, ecxfd, reflens, max_read_length,
            )
            entries.extend(
                filter_entries(
                    bam_chunk, ecxfd, targets, samfilters, min_map_overlap,
                ),
            )
    return bam_header_string, entries


def make_entry_dispatchers(entries, ecx):
    """Prepare entry descriptions for ambiguity resolution"""
    entry_dispatcher = DataFrame(
        columns=["qname", "rname", "mappos", "flag", "entry"],
        data=[[e.qname, e.reference_name, e.pos, e.flag, e] for e in entries],
    )
    entry_dispatcher["prime"] = (entry_dispatcher["flag"] & 0x8000 == 0) * 2 + 3
    rname_dispatcher = merge(
        entry_dispatcher[["qname", "rname", "prime"]],
        ecx[["rname", "prime", "chromosome"]].drop_duplicates(), how="outer",
    )
    arm_counts = (
        rname_dispatcher.groupby("qname")[["chromosome", "prime"]]
        .nunique().reset_index()
    )
    valid_qnames = arm_counts.loc[
        (arm_counts["chromosome"]==1) & (arm_counts["prime"]==1), "qname",
    ]
    assert valid_qnames.value_counts().max() == 1
    return entry_dispatcher, valid_qnames


def get_unambiguous_entries(entries, ecx):
    """Subset candidate entries to those that map unambiguously"""
    entry_dispatcher, valid_qnames = make_entry_dispatchers(entries, ecx)
    chromosomes = set(ecx["chromosome"].drop_duplicates())
    for qname in progressbar(valid_qnames, desc="Filtering", unit="read"):
        entry_mappings = entry_dispatcher[entry_dispatcher["qname"]==qname]
        entry_mapped_to_main = entry_mappings["rname"].isin(chromosomes)
        if entry_mapped_to_main.any(): # prefer canonical chromosomes
            target_entry_mappings = entry_mappings[entry_mapped_to_main]
        else: # fall back to forks / subtelomeres
            target_entry_mappings = entry_mappings
        rnames = target_entry_mappings["rname"].drop_duplicates()
        primes = target_entry_mappings["prime"].drop_duplicates()
        if (len(rnames) == 1) and (len(primes) == 1):
            if primes.iloc[0] == 3: # find innermost mappos on same q arm
                target_mappos = target_entry_mappings["mappos"].min()
            else: # find innermost mappos on same p arm
                target_mappos = target_entry_mappings["mappos"].max()
            entry_candidates = target_entry_mappings.loc[
                target_entry_mappings["mappos"]==target_mappos, "entry",
            ]
            if len(entry_candidates) == 1:
                entry = entry_candidates.iloc[0]
                if entry.flag & 0x800 == 0: # non-supplementary
                    yield entry


def main(bam, index, flags, flag_filter, min_quality, max_read_length, min_map_overlap, min_subtelomere_overlap, min_telomere_overlap, target, file=stdout, **kwargs):
    # dispatch data to subroutines:
    ecxfd = load_index(index, as_filter_dict=True)
    ecx = load_index(index, as_filter_dict=False)
    bam_header_string, entries = parse_bam_with_ambiguity(
        bam, ecxfd, max_read_length, min_map_overlap,
        targets={interpret_flags(t) for t in target},
        samfilters=[flags, flag_filter, 0], # allow any quality on first pass
    )
    print(bam_header_string, file=file)
    n_orphaned_entries = 0
    for entry in get_unambiguous_entries(entries, ecx):
        if entry.query_length - entry.reference_length >= min_telomere_overlap:
            if entry.reference_length >= min_subtelomere_overlap:
                if entry.mapq >= min_quality: # enforce quality on second pass
                    print(entry.to_string(), file=file)
                    if entry.seq is None:
                        n_orphaned_entries += 1
    if n_orphaned_entries == 0:
        return 0
    else:
        warning = [
            f"CRITICAL: {n_orphaned_entries} entries have no sequence data;",
            "          this will cause information loss downstream.",
            "          Please re-submit a SAM/BAM with sequences reported for ",
            "          all alignments.",
        ]
        print("\n".join(warning), file=stderr)
        return 1

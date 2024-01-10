from sys import stdout, stderr
from pysam import AlignmentFile
from numpy import array, where, errstate
from re import search, split
from edgecaselib.formats import load_index, filter_bam, interpret_flags
from edgecaselib.util import progressbar


__doc__ = """edgeCase tailchopper: selection of overhanging heads/tails of reads

Usage: {0} tailchopper -x filename [-t targetspec]
       {1}             [-f flagspec]... [-F flagspec]... [-q integer] <bam>

Output:
    SAM-formatted file with tails of candidate reads overhanging anchors defined
    in index

Positional arguments:
    <bam>                         name of input BAM/SAM file

Required options:
    -x, --index [filename]        location of the reference .ecx index

Options:
    -t, --target [targetspec]     an ECX flag (cut relative to reference) or 'cigar' [default: tract_anchor]

Input filtering options:
    -f, --flags [flagspec]        process only entries with all these sam flags present [default: 0]
    -F, --flag-filter [flagspec]  process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]   process only entries with this MAPQ or higher [default: 0]

Notes:
  * Depending on the aligner used, MAPQ of secondary reads may have been set to
    zero regardless of real mapping quality; use this filtering option with
    caution.
"""

__docopt_converters__ = [
    lambda min_quality:
        None if (min_quality is None) else int(min_quality),
]

__docopt_tests__ = {
    lambda target:
        target in {"mask_anchor", "fork", "tract_anchor", "cigar"}:
            "unknown value of --target",
}


def get_cigar_clip_length(entry, prime):
    """Return only clipped part of sequence: deprecated method still used by kmerscanner"""
    if prime not in {5, 3}:
        raise ValueError("`prime` can only be 5 or 3")
    if prime == 5:
        cigar_clip = search(r'^(\d+[SH])+', entry.cigarstring)
    elif prime == 3:
        cigar_clip = search(r'(\d+[SH])+$', entry.cigarstring)
    if not cigar_clip:
        return 0
    else:
        return sum(
            int(clip) for clip in split(r'[SH]', cigar_clip.group())
            if clip != ""
        )


def update_aligned_segment(entry, ref_map_pos, bounds):
    """Update sequence, cigar, quality string in-place"""
    start, end = bounds
    if (end is not None) and (start is not None) and (end < start):
        start, end = end, start
    qualities_substring = entry.query_qualities[start:end]
    entry.query_sequence = entry.query_sequence[start:end]
    if entry.query_sequence:
        entry.cigarstring = str(len(entry.query_sequence)) + "S"
    else:
        entry.cigarstring = None
    if ref_map_pos is not None:
        entry.reference_start = ref_map_pos
    entry.query_qualities = qualities_substring


def cigar_chopper(entry, ecx, integer_target):
    """Return only clipped part of sequence"""
    is_q = (entry.flag & 0x8000 != 0)
    error = None
    if is_q:
        cigar_clip = search(r'(\d+[SH])+$', entry.cigarstring)
    else:
        cigar_clip = search(r'^(\d+[SH])+', entry.cigarstring)
    if not cigar_clip:
        update_aligned_segment(entry, None, (0, 0))
        error = "No clipped sequence"
    else:
        clip_length = sum(
            int(clip) for clip in split(r'[SH]', cigar_clip.group())
            if clip != ""
        )
        if clip_length > 0:
            if is_q:
                update_aligned_segment(
                    entry, entry.reference_end-1, (-clip_length, None),
                )
            else:
                update_aligned_segment(
                    entry, entry.reference_start, (None, clip_length),
                )
        else:
            update_aligned_segment(entry, None, (0, 0))
            error = "No clipped sequence"
    return entry, error


def find_map_and_cut_positions(entry, anchor_pos, is_q):
    """Find closest mapping position to anchor"""
    ref_positions = array(
        entry.get_reference_positions(full_length=True), # [read_pos->ref_pos]
        dtype=float,
    )
    try:
        if is_q:
            read_cut_pos = where(ref_positions>=anchor_pos)[0][0]
        else:
            read_cut_pos = where(ref_positions<=anchor_pos)[0][-1]
    except IndexError:
        ref_map_pos, read_cut_pos = None, None
    else:
        ref_map_pos = ref_positions[read_cut_pos].astype(int)
    return ref_map_pos, read_cut_pos


def relative_chopper(entry, ecx, integer_target):
    """Return only part of sequence past the anchor"""
    is_q = (entry.flag & 0x8000 != 0)
    error = None
    prime = 3 if is_q else 5
    indexer = (
        (ecx["rname"]==entry.reference_name) &
        (ecx["prime"]==prime) & (ecx["flag"]==integer_target)
    )
    anchor_positions = ecx.loc[indexer, "pos"]
    if len(anchor_positions) > 1:
        raise ValueError("Ambiguous index entry: {}".format(anchor_positions))
    elif len(anchor_positions) == 0:
        update_aligned_segment(entry, None, (0, 0))
        error = "No anchor data in index"
    else:
        ref_map_pos, read_cut_pos = find_map_and_cut_positions(
            entry, anchor_positions.iloc[0], is_q,
        )
        if is_q:
            update_aligned_segment(entry, ref_map_pos, (read_cut_pos, None))
        else:
            update_aligned_segment(entry, ref_map_pos, (None, read_cut_pos))
    return entry, error


def main(bam, index, flags, flag_filter, min_quality, target, file=stdout, **kwargs):
    """Interpret arguments and dispatch data to subroutines"""
    if target == "cigar":
        chopper, integer_target = cigar_chopper, None
    else:
        chopper, integer_target = relative_chopper, interpret_flags(target)
    ecx = load_index(index)
    with AlignmentFile(bam) as alignment:
        print(str(alignment.header).rstrip("\n"), file=file)
        n_skipped = 0
        bam_iterator = progressbar(
            filter_bam(alignment, [flags, flag_filter, min_quality]),
            desc="Chopping", unit="read",
        )
        with errstate(invalid="ignore"):
            for entry in bam_iterator:
                if entry.query_sequence:
                    chopped_entry, error = chopper(
                        entry, ecx, integer_target,
                    )
                    if chopped_entry.query_sequence:
                        print(chopped_entry.to_string(), file=file)
                    else:
                        n_skipped += 1
    if n_skipped:
        msg_mask = "Skipped {} reads to be safe (unsure where to chop)"
        print(msg_mask.format(n_skipped), file=stderr)
    warning = [
        "WARNING: Read mapping positions were adjusted and retained;",
        "         this is needed to comply with the SAM spec.",
        "         Do not use these positions for analyses outside of edgeCase!",
    ]
    print("\n".join(warning), file=stderr)
    return 0

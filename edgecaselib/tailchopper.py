from sys import stdout, stderr
from pysam import AlignmentFile
from re import search, split
from edgecaselib.formats import load_index, filter_bam, interpret_flags


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


def update_aligned_segment(entry, map_pos, start=None, end=None):
    """Update sequence, cigar, quality string in-place"""
    if (end is not None) and (start is not None) and (end < start):
        start, end = end, start
    qualities_substring = entry.query_qualities[start:end]
    entry.query_sequence = entry.query_sequence[start:end]
    if entry.query_sequence:
        entry.cigarstring = str(len(entry.query_sequence)) + "S"
    else:
        entry.cigarstring = None
    entry.flag |= 4
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
        update_aligned_segment(entry, None, 0, 0)
        error = "No clipped sequence"
    else:
        clip_length = sum(
            int(clip) for clip in split(r'[SH]', cigar_clip.group())
            if clip != ""
        )
        if clip_length > 0:
            if is_q:
                map_pos = entry.reference_end - 1
                update_aligned_segment(entry, map_pos, -clip_length, None)
            else:
                map_pos = entry.reference_start
                update_aligned_segment(entry, map_pos, None, clip_length)
        else:
            update_aligned_segment(entry, None, 0, 0)
            error = "No clipped sequence"
    return entry, error


def find_map_and_cut_positions(entry, anchor_pos, is_q):
    """Find closest mapping position within `relax_radius` of anchor"""
    positions = entry.get_reference_positions(full_length=True)
    try:
        cut_pos = positions.index(anchor_pos)
        if is_q:
            return cut_pos, cut_pos
        else:
            return entry.reference_start, cut_pos
    except ValueError:
        if entry.reference_start <= anchor_pos < entry.reference_end:
            ref_map_length = entry.reference_end - 1 - entry.reference_start
            try:
                read_map_length = (
                    positions.index(entry.reference_end-1) -
                    positions.index(entry.reference_start)
                )
            except ValueError:
                raise IndexError("SAM error? Last mapped position not on read")
            ref_anchor_distance = anchor_pos - entry.reference_start
            cut_pos = int(round(
                ref_anchor_distance * read_map_length / ref_map_length
            ))
            if is_q:
                return cut_pos, cut_pos
            else:
                return entry.reference_start, cut_pos
        elif is_q:
            if entry.reference_start >= anchor_pos:
                return entry.reference_start, entry.reference_start
            elif entry.reference_end < anchor_pos:
                raise ValueError("Anchor position beyond mapped portion")
        else:
            if entry.reference_end < anchor_pos:
                return entry.reference_start, entry.reference_end - 1
            elif entry.reference_start > anchor_pos:
                raise ValueError("Anchor position beyond mapped portion")
    raise ValueError("Satisfactory cutting position not found")


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
        update_aligned_segment(entry, 0, 0)
        error = "No anchor data in index"
    else:
        anchor_pos = anchor_positions.iloc[0]
        try:
            map_pos, cut_pos = find_map_and_cut_positions(
                entry, anchor_pos, is_q
            )
        except ValueError as e:
            update_aligned_segment(entry, None, 0, 0)
            error = str(e)
        else:
            if is_q:
                update_aligned_segment(entry, map_pos, cut_pos, None)
            else:
                update_aligned_segment(entry, map_pos, None, -cut_pos)
    return entry, error


def main(bam, index, flags, flags_any, flag_filter, min_quality, target, jobs=1, file=stdout, **kwargs):
    """Interpret arguments and dispatch data to subroutines"""
    if target == "cigar":
        chopper, integer_target = cigar_chopper, None
    else:
        chopper, integer_target = relative_chopper, interpret_flags(target)
    ecx = load_index(index)
    with AlignmentFile(bam) as alignment:
        print(str(alignment.header).rstrip("\n"), file=file)
        samfilters = [flags, flags_any, flag_filter, min_quality]
        for entry in filter_bam(alignment, samfilters):
            if entry.query_sequence:
                chopped_entry, error = chopper(
                    entry, ecx, integer_target
                )
                if chopped_entry.query_sequence:
                    print(chopped_entry.to_string(), file=file)
                else:
                    error_message = "Skipping {} ({})".format(
                        entry.query_name, error or "Zero length chop"
                    )
                    print(error_message, file=stderr)

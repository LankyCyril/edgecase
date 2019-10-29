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


def update_aligned_segment(entry, start=None, end=None):
    """Update sequence, cigar, quality string in-place"""
    qualities_substring = entry.query_qualities[start:end]
    entry.query_sequence = entry.query_sequence[start:end]
    entry.cigarstring = None
    entry.flag |= 4
    entry.query_qualities = qualities_substring


def cigar_chopper(entry, ecx, integer_target, relax_radius=0):
    """Return only clipped part of sequence"""
    is_q = (entry.flag & 0x8000 != 0)
    if is_q:
        cigar_clip = search(r'(\d+[SH])+$', entry.cigarstring)
    else:
        cigar_clip = search(r'^(\d+[SH])+', entry.cigarstring)
    if not cigar_clip:
        update_aligned_segment(entry, 0, 0)
    else:
        clip_length = sum(
            int(clip) for clip in split(r'[SH]', cigar_clip.group())
            if clip != ""
        )
        if clip_length > 0:
            if is_q:
                update_aligned_segment(entry, -clip_length, None)
            else:
                update_aligned_segment(entry, None, clip_length)
        else:
            update_aligned_segment(entry, 0, 0)
    return entry


def find_closest_position(positions, anchor_pos, relax_radius):
    """Find closest mapping position within `relax_radius` of anchor"""
    for step in range(0, relax_radius+1):
        for adjusted_anchor_pos in anchor_pos-step, anchor_pos, anchor_pos+step:
            try:
                read_pos = positions.index(adjusted_anchor_pos)
            except ValueError:
                pass
            else:
                return read_pos
    else:
        raise ValueError("Mapped position not found within radius")


def relative_chopper(entry, ecx, integer_target, relax_radius):
    """Return only part of sequence past the anchor"""
    is_q = (entry.flag & 0x8000 != 0)
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
    else:
        anchor_pos = anchor_positions.iloc[0]
        positions = entry.get_reference_positions(full_length=True)
        try:
            read_pos = find_closest_position(
                positions, anchor_pos, relax_radius=relax_radius
            )
        except ValueError:
            update_aligned_segment(entry, 0, 0)
        else:
            if is_q:
                update_aligned_segment(entry, read_pos, None)
            else:
                update_aligned_segment(entry, None, -read_pos)
    return entry


def main(bam, index, flags, flags_any, flag_filter, min_quality, target, relax_radius, jobs=1, file=stdout, **kwargs):
    """Interpret arguments and dispatch data to subroutines"""
    if target == "cigar":
        if relax_radius != 0:
            raise ValueError("Cannot relax radius with --target=cigar")
        chopper, integer_target = cigar_chopper, None
    else:
        chopper, integer_target = relative_chopper, interpret_flags(target)
    ecx = load_index(index)
    with AlignmentFile(bam) as alignment:
        print(str(alignment.header).rstrip("\n"), file=file)
        samfilters = [flags, flags_any, flag_filter, min_quality]
        for entry in filter_bam(alignment, samfilters):
            if entry.query_sequence:
                chopped_entry = chopper(
                    entry, ecx, integer_target, relax_radius
                )
                if chopped_entry.query_sequence:
                    print(chopped_entry.to_string(), file=file)
                else:
                    print("Unchoppable:", entry.query_name, file=stderr)

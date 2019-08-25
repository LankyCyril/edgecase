from sys import stdout, stderr
from types import SimpleNamespace
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


def cigar_chopper(entry, ecx, target, integer_target):
    """Return only clipped part of sequence"""
    name_mask = "{}/chop=cigar/flag={}"
    is_q = (entry.flag & 0x8000 != 0)
    if is_q:
        cigar_clip = search(r'(\d+[SH])+$', entry.cigarstring)
    else:
        cigar_clip = search(r'^(\d+[SH])+', entry.cigarstring)
    if not cigar_clip:
        return SimpleNamespace(
            name=name_mask.format(entry.query_name, entry.flag),
            sequence=""
        )
    else:
        clip_length = sum(
            int(clip) for clip in split(r'[SH]', cigar_clip.group())
            if clip != ""
        )
        if clip_length > 0:
            if is_q:
                return SimpleNamespace(
                    name=name_mask.format(entry.query_name, entry.flag),
                    sequence=entry.query_sequence[-clip_length:]
                )
            else:
                return SimpleNamespace(
                    name=name_mask.format(entry.query_name, entry.flag),
                    sequence=entry.query_sequence[:clip_length]
                )
        else:
            return SimpleNamespace(
                name=name_mask.format(entry.query_name, entry.flag),
                sequence=""
            )


def relative_chopper(entry, ecx, target, integer_target):
    """Return only part of sequence past the anchor"""
    name_mask = "{}/chop=" + target + "/flag={}"
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
        return SimpleNamespace(
            name=name_mask.format(entry.query_name, entry.flag), sequence=""
        )
    else:
        anchor_pos = anchor_positions.iloc[0]
        positions = entry.get_reference_positions(full_length=True)
        try:
            read_pos = positions.index(anchor_pos)
            if is_q:
                sequence = entry.query_sequence[read_pos:]
            else:
                sequence = entry.query_sequence[:-read_pos]
        except ValueError:
            sequence = ""
        return SimpleNamespace(
            name=name_mask.format(entry.query_name, entry.flag),
            sequence=sequence
        )


def main(bam, index, flags, flags_any, flag_filter, min_quality, target, file=stdout, **kwargs):
    """Interpret arguments and dispatch data to subroutines"""
    if target == "cigar":
        chopper, integer_target = cigar_chopper, None
    else:
        chopper, integer_target = relative_chopper, interpret_flags(target)
    ecx = load_index(index)
    with AlignmentFile(bam) as alignment:
        samfilters = flags, flags_any, flag_filter, min_quality
        for entry in filter_bam(alignment, samfilters):
            chopped_entry = chopper(entry, ecx, target, integer_target)
            if len(chopped_entry.sequence):
                print(">" + chopped_entry.name, file=file)
                print(chopped_entry.sequence, file=file)
            else:
                warn_mask = "WARNING: omitting {} chopped to zero length"
                print(warn_mask.format(chopped_entry.name), file=stderr)

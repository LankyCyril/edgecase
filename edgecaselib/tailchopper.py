from sys import stdout
from pysam import AlignmentFile
from re import search, split
from edgecaselib.formats import load_index, filter_bam, interpret_flags
from edgecaselib.util import get_executable
from tempfile import TemporaryDirectory


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


def cigar_chopper(entry, ecx, integer_target):
    """Return only clipped part of sequence"""
    is_chopped = True
    is_q = (entry.flag & 0x8000 != 0)
    if is_q:
        cigar_clip = search(r'(\d+[SH])+$', entry.cigarstring)
    else:
        cigar_clip = search(r'^(\d+[SH])+', entry.cigarstring)
    if not cigar_clip:
        is_chopped = False
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
            is_chopped = False
    return entry, is_chopped


def relative_chopper(entry, ecx, integer_target):
    """Return only part of sequence past the anchor"""
    is_chopped = True
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
        is_chopped = False
    else:
        anchor_pos = anchor_positions.iloc[0]
        positions = entry.get_reference_positions(full_length=True)
        try:
            read_pos = positions.index(anchor_pos)
        except ValueError:
            is_chopped = False
        else:
            if is_q:
                update_aligned_segment(entry, read_pos, None)
            else:
                update_aligned_segment(entry, None, -read_pos)
    return entry, is_chopped


def prechop(alignment, samfilters, ecx, chopper, integer_target):
    """First round of chopping: chop only those sequences that have exact alignment to the reference anchor"""
    chopped_entries, unchopped_entries = [], []
    for entry in filter_bam(alignment, samfilters):
        if entry.query_sequence:
            maybe_chopped_entry, is_chopped = chopper(
                entry, ecx, integer_target
            )
            if len(maybe_chopped_entry.query_sequence) >= 6:
                if is_chopped:
                    chopped_entries.append(maybe_chopped_entry)
                elif maybe_chopped_entry.query_sequence:
                    unchopped_entries.append(maybe_chopped_entry)
    return chopped_entries, unchopped_entries


def rechop_mafft(chopped_entries, unchopped_entries, mafft, mafft_options, jobs):
    """Realign entries with MAFFT and chop previously unchopped entries"""
    with open("chopped.fa", mode="wt") as testfile:
        for entry in chopped_entries:
            if entry.reference_name == "chr4":
                print(">{}\n{}".format(entry.qname, entry.seq), file=testfile)
    with open("unchopped.fa", mode="wt") as testfile:
        for entry in unchopped_entries:
            if entry.reference_name == "chr4":
                print(">{}\n{}".format(entry.qname, entry.seq), file=testfile)
    exit(0)


def main(bam, index, flags, flags_any, flag_filter, min_quality, target, min_length, mafft, mafft_options, jobs=1, file=stdout, **kwargs):
    """Interpret arguments and dispatch data to subroutines"""
    if target == "cigar":
        chopper, integer_target = cigar_chopper, None
    else:
        chopper, integer_target = relative_chopper, interpret_flags(target)
    mafft = get_executable("mafft", mafft)
    ecx = load_index(index)
    with AlignmentFile(bam) as alignment:
        #print(str(alignment.header).rstrip("\n"), file=file)
        samfilters = [flags, flags_any, flag_filter, min_quality]
        chopped_entries, unchopped_entries = prechop(
            alignment, samfilters, ecx, chopper, integer_target
        )
    if len(chopped_entries) and len(unchopped_entries):
        chopped_entries = rechop_mafft(
            chopped_entries, unchopped_entries, mafft, mafft_options, jobs
        )
    for chopped_entry in chopped_entries:
        print(chopped_entry.to_string(), file=file)

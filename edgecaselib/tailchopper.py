from sys import stdout, stderr
from types import SimpleNamespace
from edgecaselib.util import ReadFileChain
from edgecaselib.tailpuller import is_good_entry
from pysam import AlignmentFile
from re import search, split


def chop(entry, prime):
    """Return only clipped part of sequence"""
    if prime not in {5, 3}:
        raise ValueError("`prime` can only be 5 or 3")
    if prime == 5:
        cigar_clip = search(r'^(\d+[SH])+', entry.cigarstring)
    elif prime == 3:
        cigar_clip = search(r'(\d+[SH])+$', entry.cigarstring)
    if not cigar_clip:
        return SimpleNamespace(name=entry.query_name, sequence="")
    else:
        clip_length = sum(
            int(clip) for clip in split(r'[SH]', cigar_clip.group())
            if clip != ""
        )
        assert clip_length > 0, "Clip length <=0, malformed CIGAR?"
        if prime == 5:
            return SimpleNamespace(
                name=entry.query_name,
                sequence=entry.query_sequence[:clip_length]
            )
        elif prime == 3:
            return SimpleNamespace(
                name=entry.query_name,
                sequence=entry.query_sequence[-clip_length:]
            )


def main(bams, prime, file=stdout, **kwargs):
    with ReadFileChain(bams, AlignmentFile) as bam_data:
        for entry in bam_data:
            if is_good_entry(entry):
                chopped_entry = chop(entry, prime)
                if len(chopped_entry.sequence):
                    print(">" + chopped_entry.name, file=file)
                    print(chopped_entry.sequence, file=file)
                else:
                    warn_mask = "WARNING: omitting {} chopped to zero length"
                    print(warn_mask.format(chopped_entry.name), file=stderr)

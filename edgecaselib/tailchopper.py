from types import SimpleNamespace
from edgecaselib.io import ReadFileChain
from edgecaselib.tailpuller import get_anchors, is_good_entry
from pysam import AlignmentFile


def chop(entry, anchors, prime):
    """Return only part of sequence extending past anchor"""
    if prime == 5:
        chop_pos = anchors.loc[entry.reference_name, "5prime"]
        return SimpleNamespace(
            name=entry.query_name,
            sequence=entry.query_sequence[:chop_pos]
        )
    elif prime == 3:
        chop_pos = anchors.loc[entry.reference_name, "3prime"]
        return SimpleNamespace(
            name=entry.query_name,
            sequence=entry.query_sequence[chop_pos:]
        )


def main(args):
    anchors = get_anchors(args.reference)
    with ReadFileChain(args.bams, AlignmentFile) as bam_data:
        for entry in bam_data:
            if is_good_entry(entry):
                chopped_entry = chop(entry, anchors, args.prime)
                print(">" + chopped_entry.name)
                print(chopped_entry.sequence)

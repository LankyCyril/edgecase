#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager, ExitStack
from itertools import chain, takewhile, filterfalse
from tqdm import tqdm
from pysam import FastxFile, AlignmentFile
from re import compile, IGNORECASE
from pandas import DataFrame
from types import SimpleNamespace

USAGE = "python3 {} [options] bams > sam/fasta".format(__file__)

ARG_RULES = {
    ("bams",): {
        "help": "name of input BAM/SAM files",
        "nargs": "+"
    },
    ("-r", "--reference"): {
        "help": "reference FASTA (required)",
        "required": True, "metavar": "R"
    },
    ("-t", "--target"): {
        "help": "what to output: 5AC, 3AC, 5OOB, 3OOB (required)",
        "metavar": "T"
    },
    ("-f", "--format"): {
        "help": "output format (FASTA, SAM) (required)",
        "required": True, "metavar": "F"
    },
    ("-j", "--jobs"): {
        "help": "number of jobs to run in parallel (default 1)",
        "default": 1, "type": int, "metavar": "J"
    }
}

MAINCHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


@contextmanager
def AlignmentFileChain(bams):
    """Chain BAM records from all filenames in list bams, replicating behavior of pysam.AlignmentFile"""
    with ExitStack() as stack:
        yield chain(*(
            stack.enter_context(AlignmentFile(bam))
            for bam in bams
        ))


def get_anchors(reference):
    """Get coordinates of hard-masked bounds at each end of each main chromosome"""
    pattern = compile(r'[^n]', flags=IGNORECASE)
    bar = tqdm(desc="Finding anchors", total=len(MAINCHROMS), unit="chromosome")
    anchor_data = {}
    with FastxFile(reference) as genome:
        for entry in genome:
            if entry.name in MAINCHROMS:
                bar.update()
                bound_5prime = pattern.search(entry.sequence).span()[0]
                bound_3prime = (
                    len(entry.sequence) -
                    pattern.search(entry.sequence[::-1]).span()[0]
                )
                anchor_data[entry.name] = bound_5prime, bound_3prime
    return DataFrame(data=anchor_data, index=["5prime", "3prime"]).T


def is_good_entry(entry):
    """Simple filter"""
    if entry.is_unmapped or entry.is_secondary or entry.is_supplementary:
        return False
    elif entry.reference_name in MAINCHROMS:
        return True
    else:
        return False


def filter_bams(bams, anchors, target):
    """Only pass reads extending past anchors"""
    isnone = lambda p: p is None
    with AlignmentFileChain(args.bams) as alignment:
        for entry in alignment:
            if is_good_entry(entry):
                positions = entry.get_reference_positions(full_length=True)
                left_clip = sum(
                    True for _ in takewhile(isnone, positions)
                )
                left_mappos = next(filterfalse(isnone, positions))
                right_clip = sum(
                    True for _ in takewhile(isnone, reversed(positions))
                )
                right_mappos = next(filterfalse(isnone, reversed(positions)))
                if target[0] == "5":
                    anchor = anchors.loc[entry.reference_name, "5prime"]
                    if left_mappos - left_clip < anchor:
                        yield entry
                if target[0] == "3":
                    anchor = anchors.loc[entry.reference_name, "3prime"]
                    if right_mappos + right_clip > anchor:
                        yield entry


def chop(entry, anchors, target):
    """Return only part of sequence extending past anchor"""
    if target[0] == "5":
        chop_pos = anchors.loc[entry.reference_name, "5prime"]
        return SimpleNamespace(
            name=entry.query_name,
            sequence=entry.query_sequence[:chop_pos]
        )
    if target[0] == "3":
        chop_pos = anchors.loc[entry.reference_name, "3prime"]
        return SimpleNamespace(
            name=entry.query_name,
            sequence=entry.query_sequence[chop_pos:]
        )


def main(args):
    """Dispatch data to subroutines"""
    anchors = get_anchors(args.reference)
    if args.format == "SAM":
        with AlignmentFile(args.bams[0]) as bam:
            print(str(bam.header).rstrip("\n"))
    for entry in filter_bams(args.bams, anchors, args.target):
        if args.target.endswith("AC"):
            print(entry.to_string())
        elif args.target.endswith("OOB"):
            chopped_entry = chop(entry, anchors, args.target)
            print(">{}\n{}".format(chopped_entry.name, chopped_entry.sequence))


if __name__ == "__main__":
    parser = ArgumentParser(usage=USAGE)
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    if args.target.endswith("AC") and (args.format == "FASTA"):
        raise NotImplementedError("{} in FASTA format".format(args.target))
    elif args.target.endswith("OOB") and (args.format == "SAM"):
        raise NotImplementedError("{} in SAM format".format(args.target))
    else:
        main(args)

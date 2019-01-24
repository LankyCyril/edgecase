#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager, ExitStack
from itertools import chain
from tqdm import tqdm
from pysam import FastxFile, AlignmentFile
from re import compile, IGNORECASE
from pandas import DataFrame

USAGE = "python3 {} --ref [reffile] [other options] bams > bam".format(__file__)

ARG_RULES = {
    ("bams",): {
        "help": "name of input BAM files",
        "nargs": "+"
    },
    ("-r", "--reference"): {
        "help": "reference FASTA (required)",
        "required": True, "metavar": "R"
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
    anchor_data = {}
    bar = tqdm(desc="Finding anchors", total=len(MAINCHROMS), unit="chromosome")
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


def main(args):
    """Dispatch data to subroutines"""
    anchors = get_anchors(args.reference)


if __name__ == "__main__":
    parser = ArgumentParser(usage=USAGE)
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    main(args)

#!/usr/bin/env python3
from argparse import ArgumentParser

ARG_RULES = {
    ("input-fastq",): {
        "help": "name of input FASTQ file",
    },
    ("--kmer",): {
        "help": "target kmer sequence (default TTAGGG)",
        "default": "TTAGGG", "metavar": "M"
    },
    ("--background-kmers",): {
        "help": "number of random kmers for significance threshold estimation (default 10)",
        "default": 10, "type": int, "metavar": "B"
    },
    ("--density-kde-plot",): {
        "help": "if filename provided, will plot densities of checked kmers",
        "default": None, "metavar": "P"
    },
    ("-n", "--num-reads"): {
        "help": "process only N reads (default all)",
        "default": None, "type": int, "metavar": "N"
    },
    ("-j", "--jobs"): {
        "help": "number of jobs to run in parallel (default 1)",
        "default": 1, "type": int, "metavar": "J"
    }
}

def main(args):
    args = get_args()
    return 0

if __name__ == "__main__":
    parser = ArgumentParser()
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    returncode = main(args)
    exit(returncode)

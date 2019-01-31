#!/usr/bin/env python3
from argparse import ArgumentParser
from edgecase.kmerscanner import fastq_scanner
from sys import stderr

USAGE = "python3 {} [options] fastq > txt".format(__file__)

ARG_RULES = {
    ("fastqs",): {
        "help": "name of input FASTQ file",
        "nargs": "+"
    },
    ("--kmer",): {
        "help": "target kmer sequence (default TTAGGG)",
        "default": "TTAGGG", "metavar": "M"
    },
    ("--head-test",): {
        "help": "length of head to use for density filter (if specified)",
        "default": None, "type": int, "metavar": "H"
    },
    ("--tail-test",): {
        "help": "length of tail to use for density filter (if specified)",
        "default": None, "type": int, "metavar": "T"
    },
    ("-c", "--cutoff"): {
        "help": "use hard cutoff for density (default None)",
        "default": None, "type": float, "metavar": "C"
    },
    ("-w", "--window-size"): {
        "help": "size of the rolling window (default 120)",
        "default": 120, "type": int, "metavar": "W"
    },
    ("-n", "--num-reads"): {
        "help": "process at most N reads (default all)",
        "default": None, "type": int, "metavar": "N"
    },
    ("-j", "--jobs"): {
        "help": "number of jobs to run in parallel (default 1)",
        "default": 1, "type": int, "metavar": "J"
    }
}


def main(args):
    """Dispatch data to subroutines"""
    pattern = "FIXME"
    # decide on density cutoff:
    if (args.head_test is not None) or (args.tail_test is not None):
        cutoff = args.cutoff
    else: # no cutoffs, process all reads
        cutoff = None
    # scan fastq for target kmer query, parallelizing on reads:
    scanner = fastq_scanner(
        args.fastqs, pattern,
        window_size=args.window_size, head=args.head_test, tail=args.tail_test,
        cutoff=cutoff, num_reads=args.num_reads, jobs=args.jobs,
    )
    # output densities of reads that pass filter:
    for read_name, density_array, passes_filter in scanner:
        if passes_filter:
            print(read_name, *density_array, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser(usage=USAGE)
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    if (args.head_test is not None) and (args.tail_test is not None):
        raise ValueError("Can only specify one of --head, --tail")
    elif (args.cutoff is not None) and (args.head_test is None) and (args.tail_test is None):
        raise ValueError("--cutoff has no effect without a head/tail test")
    elif (args.cutoff is None) and ((args.head_test is not None) or (args.tail_test is not None)):
        print("Warning: --head/--tail has no effect without --cutoff", file=stderr)
    main(args)

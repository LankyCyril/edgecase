#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import product
from regex import compile, IGNORECASE

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

class KmerIdentity:
    """Hold all possible kmers of length k and their groupings by circular shift"""
    _one2many = {}
    _many2one = {}
 
    def __init__(self, k=6, alphabet=list("ACGT")):
        """Generate all possible kmers and group into identical by circular shift"""
        for kmer in map("".join, product(alphabet, repeat=k)):
            shifted_kmers = {kmer[i:]+kmer[:i] for i in range(k)}
            for shifted_kmer in shifted_kmers:
                if shifted_kmer in self._one2many:
                    anchor_kmer = shifted_kmer
                    self._one2many[anchor_kmer] |= shifted_kmers
                    break
            else:
                anchor_kmer = shifted_kmer
                self._one2many[anchor_kmer] = shifted_kmers
            for shifted_kmer in shifted_kmers:
                self._many2one[shifted_kmer] = anchor_kmer
 
    def pattern(self, kmer, flags=IGNORECASE):
        """Return a pattern matching kmer and its shifts"""
        anchor_kmer = self._many2one[kmer]
        shifted_kmers = self._one2many[anchor_kmer]
        return compile(r'|'.join(shifted_kmers), flags=flags)

def main(args):
    kmer_identity = KmerIdentity(k=len(args.kmer))
    return 0

if __name__ == "__main__":
    parser = ArgumentParser()
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    returncode = main(args)
    exit(returncode)

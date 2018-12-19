#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import product
from regex import compile, IGNORECASE
from random import sample
from multiprocessing import Pool
from tqdm import tqdm
from pysam import FastxFile
from functools import partial

ARG_RULES = {
    ("fastq",): {
        "help": "name of input FASTQ file",
    },
    ("--kmer",): {
        "help": "target kmer sequence (default TTAGGG)",
        "default": "TTAGGG", "metavar": "M"
    },
    ("-w", "--window-size"): {
        "help": "size of the rolling window (default 120)",
        "default": 120, "type": int, "metavar": "W"
    },
    ("-b", "--n-background-kmers"): {
        "help": "number of random kmers for significance threshold estimation (default 9)",
        "default": 9, "type": int, "metavar": "B"
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
 
    def anchors(self, exclude=()):
        """Return all anchor kmers except for the ones matching `exclude`"""
        anchor_kmers = set(self._one2many.keys())
        excluded_kmers = {self._many2one[kmer] for kmer in exclude}
        return anchor_kmers - excluded_kmers


def choose_background_kmers(kmer_identity, kmer, n_background_kmers):
    """Choose background kmers at random, throw sensible errors"""
    population = kmer_identity.anchors(exclude={kmer})
    if n_background_kmers > len(population):
        raise ValueError(
            "Requested number of background kmers " +
            "bigger than number of possible kmers"
        )
    else:
        return sample(population, n_background_kmers)


def calculate_density(read, overlapped, window_size):
    return read.name, [len(read.sequence)]


def pattern_scanner(read_iterator, pattern, overlapped=True, window_size=120, jobs=1, num_reads=None, desc=None):
    """Calculate density of pattern hits in a rolling window along each read"""
    with Pool(jobs) as pool:
        density_calculator = partial(
            calculate_density, overlapped=overlapped, window_size=window_size
        )
        read_density_iterator = pool.imap_unordered(
            density_calculator, read_iterator
        )
        if desc is None:
            desc = pattern.pattern
        yield from tqdm(
            read_density_iterator, desc=desc, unit="read", total=num_reads
        )


def main(args):
    kmer_identity = KmerIdentity(k=len(args.kmer))
    background_kmers = choose_background_kmers(
        kmer_identity, args.kmer, args.n_background_kmers
    )
    for kmer in [args.kmer]+background_kmers:
        with FastxFile(args.fastq) as read_iterator:
            scanner = pattern_scanner(
                read_iterator, kmer_identity.pattern(kmer),
                window_size=args.window_size, jobs=args.jobs,
                num_reads=args.num_reads, desc=kmer
            )
            for read_name, kmer_density in scanner:
                print(kmer, read_name, *kmer_density, sep="\t")
    return 0


if __name__ == "__main__":
    parser = ArgumentParser()
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    returncode = main(args)
    exit(returncode)

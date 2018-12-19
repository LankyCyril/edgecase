#!/usr/bin/env python3
from sys import stderr
from argparse import ArgumentParser
from itertools import product, islice
from regex import compile, IGNORECASE
from random import sample
from multiprocessing import Pool
from tqdm import tqdm
from pysam import FastxFile
from functools import partial
from numpy import array, zeros, cumsum
from pandas import DataFrame, concat
from matplotlib.pyplot import subplots, switch_backend
from seaborn import kdeplot
from collections import OrderedDict

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
    ("-n", "--num-reads"): {
        "help": "process at most N reads (default all)",
        "default": None, "type": int, "metavar": "N"
    },
    ("-j", "--jobs"): {
        "help": "number of jobs to run in parallel (default 1)",
        "default": 1, "type": int, "metavar": "J"
    },
    ("--dump",): {
        "help": "(temporary) dump all output to text file",
        "action": "store_true"
    },
    ("--density-kde-plot",): {
        "help": "if filename provided, will plot densities of checked kmers",
        "default": None, "metavar": "P"
    },
}

ALPHABET = list("ACGT")
COMPLEMENT = dict(zip(ALPHABET, ALPHABET[::-1]))


def revcomp(sequence):
    """Reverse complement a nucleotide sequence"""
    return "".join(COMPLEMENT[c] for c in reversed(sequence))


class KmerIdentity:
    """Holds all possible kmers of length k and their groupings by circular shift"""
    _one2many = {}
    _many2one = {}
 
    def __init__(self, k, alphabet=ALPHABET):
        """Generate all possible kmers and group into identical by circular shift"""
        kmer_iterator = map("".join, product(alphabet, repeat=k))
        total = len(alphabet)**k
        desc = "Generating kmer identity map"
        for kmer in tqdm(kmer_iterator, desc=desc, total=total):
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
        excluded_anchors = {self._many2one[kmer] for kmer in exclude}
        return anchor_kmers - excluded_anchors


def choose_background_kmers(kmer_identity, n_background_kmers, exclude):
    """Choose background kmers at random, throw sensible errors"""
    population = kmer_identity.anchors(exclude=exclude)
    if n_background_kmers > len(population):
        raise ValueError(
            "Requested number of background kmers " +
            "bigger than number of possible kmers"
        )
    else:
        return sample(population, n_background_kmers)


def calculate_density(read, patterns, overlapped, window_size):
    """Calculate density of pattern hits in a rolling window along given read"""
    read_length = len(read.sequence)
    if read_length == 0: # nothing to search for
        density_array = zeros((1, len(patterns)))
    else: # create array of hits; axis 0 is positions, axis 1 is patterns:
        canvas = zeros((read_length, len(patterns)), dtype=bool)
        for i, pattern in enumerate(patterns.values()):
            pattern_positions = array([
                match.start() for match
                in pattern.finditer(read.sequence, overlapped=overlapped)
            ])
            if len(pattern_positions):
                canvas[pattern_positions, i] = True
        if read_length <= window_size: # use one window:
            # [None, :] ensures two dimensions in output:
            density_array = (canvas.sum(axis=0) / read_length)[None, :]
        else: # use rolling window:
            roller = cumsum(canvas, axis=0)
            roller[window_size:] = roller[window_size:] - roller[:-window_size]
            density_array = roller[window_size-1:] / window_size
    return read.name, list(patterns.keys()), density_array


def pattern_scanner(read_iterator, patterns, overlapped=True, window_size=120, num_reads=None, jobs=1):
    """Calculate density of pattern hits in a rolling window along each read"""
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density, patterns=patterns,
            overlapped=overlapped, window_size=window_size
        )
        # lazy multiprocess evaluation:
        read_density_iterator = pool.imap_unordered(
            density_calculator, read_iterator
        )
        # iterate pairs (read.name, density_array), same as calculate_density():
        yield from tqdm(
            read_density_iterator, desc="Calculating kmer density",
            unit="read", total=num_reads
        )


def fastq_scanner(fastq, patterns, overlapped=True, window_size=120, num_reads=None, jobs=1):
    """Thin wrapper over pattern_scanner() providing read_iterator from fastq file"""
    with FastxFile(args.fastq) as read_iterator:
        # only take the first num_reads entries (`None` takes all):
        capped_read_iterator = islice(read_iterator, args.num_reads)
        # parallelize scan (all patterns at once):
        yield from pattern_scanner(
            capped_read_iterator, patterns, jobs=args.jobs,
            window_size=args.window_size, num_reads=args.num_reads
        )


def to_narrow_dataframe(scanner):
    """Convert multilevel data reads->kmers->densities to narrow-form DataFrame"""
    sections = []
    for read_name, kmers, density_array in scanner:
        section = DataFrame(data=density_array, columns=kmers)
        section.index.name = "position"
        section = section.reset_index().melt(
            id_vars="position", var_name="kmer", value_name="density"
        )
        section["read_name"] = read_name
        sections.append(section)
    print("Merging DataFrame entries... ", end="", file=stderr, flush=True)
    matrix = concat(sections)[["kmer", "read_name", "position", "density"]]
    print("done", file=stderr, flush=True)
    return matrix


def dump_scan_results(scanner):
    """Print scanner results to stdout"""
    for read_name, kmers, density_array in scanner:
        for kmer, kmer_density in zip(kmers, density_array.T):
            print(read_name, kmer, *kmer_density, sep="\t")


def plot_density_kdes(densities_matrix, target_kmers, imgfile, figsize=(10, 12)):
    """Plot combined density kde to imgfile"""
    print("Plotting density KDEs... ", end="", file=stderr, flush=True)
    switch_backend("Agg")
    figure, [target_ax, background_ax] = subplots(
        nrows=2, figsize=figsize, sharex=True
    )
    target_indexer = densities_matrix["kmer"].isin(target_kmers)
    kdeplot(
        densities_matrix.loc[target_indexer, "density"],
        color="green", ax=target_ax, legend=False
    )
    kdeplot(
        densities_matrix.loc[~target_indexer, "density"],
        color="gray", ax=background_ax, legend=False
    )
    title_mask = "Density distribution of {}"
    target_ax.set(
        title=title_mask.format(target_kmers)
    )
    background_ax.set(
        title=title_mask.format("background kmers")
    )
    figure.savefig(imgfile)
    print("done", file=stderr, flush=True)


def main(args):
    """Dispatch data to subroutines"""
    target_kmer, target_rc_kmer = args.kmer, revcomp(args.kmer)
    # precompute circular shifts of all possible kmers of same length as target:
    kmer_identity = KmerIdentity(k=len(target_kmer))
    # randomly select n_background_kmers kmers different from target:
    background_kmers = choose_background_kmers(
        kmer_identity, args.n_background_kmers,
        exclude={target_kmer, target_rc_kmer}
    )
    # prepare scanner patterns:
    patterns = OrderedDict((
        (kmer, kmer_identity.pattern(kmer))
        for kmer in [target_kmer, target_rc_kmer] + background_kmers
    ))
    # scan fastq for each unique kmer query, parallelizing on reads:
    scanner = fastq_scanner(
        args.fastq, patterns, jobs=args.jobs,
        window_size=args.window_size, num_reads=args.num_reads
    )
    if args.dump: # print everything to stdout:
        dump_scan_results(scanner)
    elif args.density_kde_plot: # make narrow-form DataFrame and plot:
        densities_matrix = to_narrow_dataframe(scanner)
        plot_density_kdes(
            densities_matrix, target_kmers={target_kmer, target_rc_kmer},
            imgfile=args.density_kde_plot
        )
    else: # perform searches without further action:
        print("Dry run", file=stderr, flush=True)
        for _ in scanner:
            pass


if __name__ == "__main__":
    parser = ArgumentParser()
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    if args.dump and args.density_kde_plot:
        raise ValueError("Cannot both dump and plot at the same time")
    else:
        main(args)

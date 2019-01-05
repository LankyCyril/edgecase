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
from pandas import DataFrame, Series, concat
from matplotlib.pyplot import subplots, switch_backend
from seaborn import kdeplot

ARG_RULES = {
    ("fastq",): {
        "help": "name of input FASTQ file",
    },
    ("--kmer",): {
        "help": "target kmer sequence (default TTAGGG)",
        "default": "TTAGGG", "metavar": "M"
    },
    ("--head",): {
        "help": "length of head to use for density filter (default 768)",
        "default": None, "type": int, "metavar": "H",
    },
    ("--tail",): {
        "help": "length of tail to use for density filter (default 768)",
        "default": None, "type": int, "metavar": "T",
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
    }
}

ALPHABET = list("ACGT")


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


def get_edge_density(sequence, pattern, overlapped, head=None, tail=None):
    """Calculate density of pattern in head or tail of sequence"""
    if len(sequence) == 0:
        return 0
    if head:
        subsequence = sequence[:head]
    elif tail:
        subsequence = sequence[-tail:]
    pattern_matches = pattern.findall(subsequence, overlapped=overlapped)
    return len(pattern_matches) / len(subsequence)


def get_edge_densities(read, patterns, overlapped, head=None, tail=None):
    """Thin wrapper for get_edge_density(): calculate edge densities of multiple patterns in read"""
    return {
        kmer: get_edge_density(read.sequence, pattern, overlapped, head, tail)
        for kmer, pattern in patterns.items()
    }


def determine_cutoff(fastq, patterns, target_kmer, overlapped=True, head=None, tail=None, num_reads=None, jobs=1):
    """Determine significant density cutoff"""
    with FastxFile(fastq) as read_iterator:
        # only take the first num_reads entries (`None` takes all):
        capped_read_iterator = islice(read_iterator, num_reads)
        # parallelize scan (all patterns at once):
        with Pool(jobs) as pool:
            # imap_unordered() only accepts single-argument functions:
            edge_densities_calculator = partial(
                get_edge_densities, patterns=patterns,
                overlapped=overlapped, head=head, tail=tail
            )
            # lazy multiprocess evaluation:
            edge_densities_iterator = pool.imap_unordered(
                edge_densities_calculator, read_iterator
            )
            # combine results into DataFrame:
            decorated_iterator = tqdm(
                edge_densities_iterator, desc="Calculating significance cutoff",
                unit="read", total=num_reads
            )
            edge_densities_dataframe = concat(
                (Series(ed) for ed in decorated_iterator),
                axis=1
            )
    print(edge_densities_dataframe.T.to_csv(sep="\t", index=None))


def calculate_density(read, pattern, overlapped, window_size, head=None, tail=None, cutoff=0):
    """Calculate density of pattern hits in a rolling window along given read"""
    edge_density = get_edge_density(
        read.sequence, pattern, overlapped, head, tail
    )
    passes_filter = (edge_density >= cutoff)
    if not passes_filter: # nothing to search for
        density_array = zeros(1)
    else: # create array of hits; axis 0 is positions, axis 1 is patterns
        read_length = len(read.sequence)
        canvas = zeros(read_length, dtype=bool)
        pattern_positions = array([
            match.start() for match
            in pattern.finditer(read.sequence, overlapped=overlapped)
        ])
        if len(pattern_positions):
            canvas[pattern_positions] = True
        if read_length <= window_size: # use one window:
            density_array = (canvas.sum(axis=0) / read_length)
        else: # use rolling window:
            roller = cumsum(canvas, axis=0)
            roller[window_size:] = roller[window_size:] - roller[:-window_size]
            density_array = roller[window_size-1:] / window_size
    return read.name, density_array, passes_filter


def pattern_scanner(read_iterator, pattern, overlapped=True, window_size=120, head=None, tail=None, cutoff=0, num_reads=None, jobs=1):
    """Calculate density of pattern hits in a rolling window along each read"""
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density, pattern=pattern,
            overlapped=overlapped, window_size=window_size,
            head=head, tail=tail, cutoff=cutoff
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


def fastq_scanner(fastq, pattern, overlapped=True, window_size=120, head=None, tail=None, cutoff=0, num_reads=None, jobs=1):
    """Thin wrapper over pattern_scanner() providing read_iterator from fastq file"""
    with FastxFile(fastq) as read_iterator:
        # only take the first num_reads entries (`None` takes all):
        capped_read_iterator = islice(read_iterator, num_reads)
        # parallelize scan:
        yield from pattern_scanner(
            capped_read_iterator, pattern, jobs=jobs,
            window_size=window_size, head=head, tail=tail, cutoff=cutoff,
            num_reads=num_reads
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
    kmer_identity = KmerIdentity(k=len(args.kmer))
    # randomly select n_background_kmers kmers different from target:
    background_kmers = choose_background_kmers(
        kmer_identity, args.n_background_kmers, exclude={args.kmer}
    )
    # prepare scanner patterns for target and background kmers:
    patterns = {
        kmer: kmer_identity.pattern(kmer)
        for kmer in [args.kmer] + background_kmers
    }
    # decide on density cutoff:
    cutoff = determine_cutoff(
        args.fastq, patterns, target_kmer=args.kmer, jobs=args.jobs,
        head=args.head, tail=args.tail,
        num_reads=args.num_reads
    )
    raise NotImplementedError("Determining cutoff")
    # scan fastq for target kmer query, parallelizing on reads:
    scanner = fastq_scanner(
        args.fastq, kmer_identity.pattern(args.kmer), jobs=args.jobs,
        window_size=args.window_size, head=args.head, tail=args.tail,
        cutoff=cutoff, num_reads=args.num_reads
    )
    # output densities of reads that pass filter:
    for read_name, density_array, passes_filter in scanner:
        if passes_filter:
            print(read_name, *density_array, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    for rule_args, rule_kwargs in ARG_RULES.items():
        parser.add_argument(*rule_args, **rule_kwargs)
    args = parser.parse_args()
    if (args.head is None) and (args.tail is None):
        raise ValueError("Must specify --head or --tail for filtering")
    elif (args.head is not None) and (args.tail is not None):
        raise ValueError("Can only specify one of --head, --tail")
    else:
        main(args)

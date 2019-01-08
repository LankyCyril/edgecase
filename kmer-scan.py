#!/usr/bin/env python3
from sys import stderr
from argparse import ArgumentParser
from itertools import product, islice
from regex import compile, IGNORECASE
from multiprocessing import Pool
from tqdm import tqdm
from pysam import FastxFile
from functools import partial
from numpy import array, zeros, cumsum, fromiter
from sklearn.mixture import GaussianMixture

ARG_RULES = {
    ("fastq",): {
        "help": "name of input FASTQ file",
    },
    ("--kmer",): {
        "help": "target kmer sequence (default TTAGGG)",
        "default": "TTAGGG", "metavar": "M"
    },
    ("--head-test",): {
        "help": "length of head to use for density filter (default 768)",
        "default": None, "type": int, "metavar": "H"
    },
    ("--tail-test",): {
        "help": "length of tail to use for density filter (default 768)",
        "default": None, "type": int, "metavar": "T"
    },
    ("--output-gmm",): {
        "help": "if filename provided, store GMM components in it",
        "default": None, "type": str, "metavar": "G"
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


def get_edge_density(read, pattern, overlapped, head=None, tail=None):
    """Calculate density of pattern in head or tail of read"""
    if len(read.sequence) == 0:
        return 0
    if head:
        subsequence = read.sequence[:head]
    elif tail:
        subsequence = read.sequence[-tail:]
    pattern_matches = pattern.findall(subsequence, overlapped=overlapped)
    return len(pattern_matches) / len(subsequence)


def output_gmm_components(gmm, X, edge_densities, tsv):
    """Save discovered components to file"""
    labels = gmm.predict(X)
    with open(tsv, mode="wt") as handle:
        for row in zip(edge_densities, labels):
            print(*row, sep="\t", file=handle)


def train_gmm(fastq, pattern, overlapped=True, head=None, tail=None, num_reads=None, output_gmm=None, jobs=1):
    """Train Gaussian Mixture to determine distribution components"""
    with FastxFile(fastq) as read_iterator:
        # only take the first num_reads entries (`None` takes all):
        capped_read_iterator = islice(read_iterator, num_reads)
        # parallelize scan (all patterns at once):
        with Pool(jobs) as pool:
            # imap_unordered() only accepts single-argument functions:
            edge_densities_calculator = partial(
                get_edge_density, pattern=pattern,
                overlapped=overlapped, head=head, tail=tail
            )
            # lazy multiprocess evaluation:
            edge_densities_iterator = pool.imap_unordered(
                edge_densities_calculator, read_iterator
            )
            # collect results into an array:
            decorated_iterator = tqdm(
                edge_densities_iterator,
                desc="Collecting edge densities for training",
                unit="read", total=num_reads
            )
            edge_densities = fromiter(decorated_iterator, dtype="float32")
    # find best n_components based on AIC:
    X = edge_densities.reshape(-1, 1)
    aics = []
    for n_components in tqdm(range(2, 6), desc="Training GMM"):
        gmm = GaussianMixture(n_components=n_components).fit(X)
        aics.append(gmm.aic(X))
    n_components = array(aics).argmin() + 2
    # train final Gaussian Mixture model and determine significant component:
    gmm = GaussianMixture(n_components=n_components).fit(X)
    target_component = gmm.predict([[edge_densities.max()]])
    if output_gmm is not None: # visualize components if requested
        output_gmm_components(gmm, X, edge_densities, tsv=output_gmm)
    return gmm, target_component


def calculate_density(read, pattern, gmm, target_component, overlapped, window_size, head=None, tail=None, cutoff=0):
    """Calculate density of pattern hits in a rolling window along given read"""
    edge_density = get_edge_density(
        read, pattern, overlapped, head, tail
    )
    passes_filter = (gmm.predict([[edge_density]])[0] == target_component)
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


def pattern_scanner(read_iterator, pattern, gmm, target_component, overlapped=True, window_size=120, head=None, tail=None, cutoff=0, num_reads=None, jobs=1):
    """Calculate density of pattern hits in a rolling window along each read"""
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density, pattern=pattern,
            overlapped=overlapped, window_size=window_size,
            head=head, tail=tail,
            gmm=gmm, target_component=target_component
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


def fastq_scanner(fastq, pattern, gmm, target_component, overlapped=True, window_size=120, head=None, tail=None, num_reads=None, jobs=1):
    """Thin wrapper over pattern_scanner() providing read_iterator from fastq file"""
    with FastxFile(fastq) as read_iterator:
        # only take the first num_reads entries (`None` takes all):
        capped_read_iterator = islice(read_iterator, num_reads)
        # parallelize scan:
        yield from pattern_scanner(
            capped_read_iterator, pattern, jobs=jobs,
            window_size=window_size, head=head, tail=tail,
            gmm=gmm, target_component=target_component,
            num_reads=num_reads
        )


def main(args):
    """Dispatch data to subroutines"""
    kmer_identity = KmerIdentity(k=len(args.kmer))
    # decide on density cutoff:
    gmm, target_component = train_gmm(
        args.fastq, kmer_identity.pattern(args.kmer), jobs=args.jobs,
        head=args.head_test, tail=args.tail_test,
        output_gmm=args.output_gmm, num_reads=args.num_reads
    )
    # scan fastq for target kmer query, parallelizing on reads:
    scanner = fastq_scanner(
        args.fastq, kmer_identity.pattern(args.kmer), jobs=args.jobs,
        window_size=args.window_size, head=args.head_test, tail=args.tail_test,
        gmm=gmm, target_component=target_component,
        num_reads=args.num_reads
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
    if (args.head_test is None) and (args.tail_test is None):
        raise ValueError("Must specify --head or --tail for filtering")
    elif (args.head_test is not None) and (args.tail_test is not None):
        raise ValueError("Can only specify one of --head, --tail")
    else:
        main(args)

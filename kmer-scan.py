#!/usr/bin/env python3
from sys import stderr
from argparse import ArgumentParser
from itertools import product, islice, chain
from contextlib import contextmanager, ExitStack
from regex import compile, IGNORECASE
from multiprocessing import Pool
from tqdm import tqdm
from pysam import FastxFile
from functools import partial
from numpy import array, zeros, cumsum, fromiter
from sklearn.mixture import GaussianMixture

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
    ("-p", "--pmax"): {
        "help": "p-value cutoff (1 - proba) for target GMM component (1e-5)",
        "default": 1e-5, "type": float, "metavar": "P"
    },
    ("-c", "--cutoff"): {
        "help": "use hard cutoff for density instead of GMM training",
        "default": None, "type": float, "metavar": "C"
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

GMM_NCOMPONENTS_RANGE = range(2, 6)
GMM_TRAIN_ROUNDS = 5


@contextmanager
def FastxChain(fastxs):
    """Chain fastx records from all filenames in list fastxs, replicating behavior of pysam.FastxFile"""
    with ExitStack() as stack:
        yield chain(*(
            stack.enter_context(FastxFile(fastx))
            for fastx in fastxs
        ))


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


def output_gmm_components(gmm, edge_densities, tsv):
    """Save discovered components to file"""
    X = edge_densities.reshape(-1, 1)
    labels = gmm.predict(X)
    probas = gmm.predict_proba(X)
    with open(tsv, mode="wt") as handle:
        for ed, label, ps in zip(edge_densities, labels, probas):
            print(ed, label, *ps, sep="\t", file=handle)


def train_gmm(mixed_distribution):
    """Train Gaussian Mixture to determine distribution components; choose model with lowest AIC"""
    aic2gmm = {}
    decorated_iterator = tqdm(
        chain(*[GMM_NCOMPONENTS_RANGE]*GMM_TRAIN_ROUNDS),
        desc="Training GMM", total=len(GMM_NCOMPONENTS_RANGE)*GMM_TRAIN_ROUNDS
    )
    for n_components in decorated_iterator:
        gmm = GaussianMixture(n_components=n_components).fit(mixed_distribution)
        aic2gmm[gmm.aic(mixed_distribution)] = gmm
    min_aic = min(aic2gmm.keys())
    return aic2gmm[min_aic]


def train_density_gmm(fastqs, pattern, overlapped=True, head=None, tail=None, num_reads=None, output_gmm=None, jobs=1):
    """Train Gaussian Mixture to determine component containing significant edge densities"""
    with FastxChain(fastqs) as read_iterator:
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
                edge_densities_calculator, capped_read_iterator
            )
            # collect results into an array:
            decorated_iterator = tqdm(
                edge_densities_iterator,
                desc="Collecting edge densities for training",
                unit="read", total=num_reads
            )
            edge_densities = fromiter(decorated_iterator, dtype="float32")
    # train Gaussian Mixture model and determine significant component:
    gmm = train_gmm(edge_densities.reshape(-1, 1))
    target_component = gmm.predict([[edge_densities.max()]])
    # visualize components if requested:
    if output_gmm is not None:
        output_gmm_components(gmm, edge_densities, tsv=output_gmm)
    return gmm, target_component


def calculate_density(read, pattern, gmm, target_component, pmax, cutoff, overlapped, window_size, head=None, tail=None):
    """Calculate density of pattern hits in a rolling window along given read"""
    if gmm: # if GMM trained, filter by predict_proba
        edge_density = get_edge_density(
            read, pattern, overlapped, head, tail
        )
        probas = gmm.predict_proba([[edge_density]])
        passes_filter = (1 - probas[0][target_component] < pmax)
    elif cutoff: # if cutoff specified, filter by hard cutoff
        edge_density = get_edge_density(
            read, pattern, overlapped, head, tail
        )
        passes_filter = (edge_density > cutoff)
    else: # otherwise, allow all
        passes_filter = True
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


def pattern_scanner(read_iterator, pattern, gmm, target_component, pmax, cutoff=None, overlapped=True, window_size=120, head=None, tail=None, num_reads=None, jobs=1):
    """Calculate density of pattern hits in a rolling window along each read"""
    with Pool(jobs) as pool:
        # imap_unordered() only accepts single-argument functions:
        density_calculator = partial(
            calculate_density, pattern=pattern,
            overlapped=overlapped, window_size=window_size,
            head=head, tail=tail, cutoff=cutoff,
            gmm=gmm, target_component=target_component, pmax=args.pmax
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


def fastq_scanner(fastqs, pattern, gmm, target_component, pmax, cutoff=None, overlapped=True, window_size=120, head=None, tail=None, num_reads=None, jobs=1):
    """Thin wrapper over pattern_scanner() providing read_iterator from fastq file"""
    with FastxChain(fastqs) as read_iterator:
        # only take the first num_reads entries (`None` takes all):
        capped_read_iterator = islice(read_iterator, num_reads)
        # parallelize scan:
        yield from pattern_scanner(
            capped_read_iterator, pattern, jobs=jobs,
            window_size=window_size, head=head, tail=tail,
            gmm=gmm, target_component=target_component, pmax=args.pmax,
            cutoff=cutoff, num_reads=num_reads
        )


def main(args):
    """Dispatch data to subroutines"""
    kmer_identity = KmerIdentity(k=len(args.kmer))
    # decide on density cutoff:
    if (args.head_test is not None) or (args.tail_test is not None):
        if args.cutoff is None: # use GMM instead of hard cutoff
            gmm, target_component = train_density_gmm(
                args.fastqs, kmer_identity.pattern(args.kmer), jobs=args.jobs,
                head=args.head_test, tail=args.tail_test,
                output_gmm=args.output_gmm, num_reads=args.num_reads
            )
            cutoff = None
        else: # use hard cutoff
            gmm, target_component, cutoff = None, None, args.cutoff
    else: # no cutoffs, process all reads
        gmm, target_component, cutoff = None, None, None
    # scan fastq for target kmer query, parallelizing on reads:
    scanner = fastq_scanner(
        args.fastqs, kmer_identity.pattern(args.kmer), jobs=args.jobs,
        window_size=args.window_size, head=args.head_test, tail=args.tail_test,
        gmm=gmm, target_component=target_component, pmax=args.pmax,
        cutoff=cutoff, num_reads=args.num_reads
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
    elif (args.cutoff is not None) and (args.output_gmm is not None):
        raise ValueError("Cannot output GMM if hard cutoff specified")
    else:
        main(args)

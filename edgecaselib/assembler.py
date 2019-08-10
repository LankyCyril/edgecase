from sys import stdout
from edgecaselib.densityplot import chromosome_natsort
from pysam import FastaFile
from collections import defaultdict
from tqdm import tqdm


def get_reference_seq(reference, chromosome):
    """Pick out only the sequence of targeted chromosome"""
    with FastaFile(reference) as reference_fasta:
        if chromosome in reference_fasta:
            return reference_fasta[chromosome]
        else:
            raise KeyError("No '{}' in reference".format(chromosome))


def count_kmers(sequence, k, desc=None):
    """Count kmers in sequence"""
    counts = defaultdict(int)
    if desc:
        sequence_iterator = tqdm(
            sequence[k:], desc="Counting kmers in {}".format(desc)
        )
    else:
        sequence_iterator = iter(sequence[k:])
    kmer = sequence[:k].upper()
    counts[kmer] = 1
    for base in sequence_iterator:
        kmer = kmer[1:] + base.upper()
        counts[kmer] += 1
    return dict(counts)


def main(bam, reference, kmer_size, chromosomes, output_prefix, jobs=1, file=stdout, **kwargs):
    # parse and check arguments:
    target_chromosomes = sorted(
        set(chromosomes.split("|")), key=chromosome_natsort
    )
    for chromosome in target_chromosomes:
        reference_seq = get_reference_seq(reference, chromosome)
        reference_kmers = count_kmers(reference_seq, kmer_size, desc=chromosome)

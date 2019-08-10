from sys import stdout
from edgecaselib.util import natsorted_chromosomes
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


def count_kmers(sequence, k, with_ordered=False, desc=None):
    """Count kmers in sequence; optionally return them in order (useful for reference sequences)"""
    counts, ordered_kmers = defaultdict(int), []
    if desc:
        sequence_iterator = tqdm(sequence[k:], desc=desc)
    else:
        sequence_iterator = iter(sequence[k:])
    kmer = sequence[:k].upper()
    counts[kmer] = 1
    if with_ordered:
        ordered_kmers.append(kmer)
    for base in sequence_iterator:
        kmer = kmer[1:] + base.upper()
        counts[kmer] += 1
        if with_ordered:
            ordered_kmers.append(kmer)
    if with_ordered:
        return dict(counts), ordered_kmers
    else:
        return dict(counts)


def main(bam, reference, kmer_size, chromosomes, output_prefix, jobs=1, file=stdout, **kwargs):
    # parse and check arguments:
    if chromosomes:
        target_chromosomes = natsorted_chromosomes(chromosomes.split("|"))
    else:
        with FastaFile(reference) as fa:
            target_chromosomes = natsorted_chromosomes(fa.references)
    for chromosome in target_chromosomes:
        reference_seq = get_reference_seq(reference, chromosome)
        reference_kmers, ordered_reference_kmers = count_kmers(
            reference_seq, kmer_size, with_ordered=True,
            desc="Counting kmers in {}".format(chromosome)
        )

from sys import stdout
from edgecaselib.util import natsorted_chromosomes
from edgecaselib.formats import interpret_flags, filter_bam
from pysam import FastaFile, AlignmentFile
from collections import defaultdict
from tqdm import tqdm


def get_reference_seq(reference, chromosome, minpos=None, maxpos=None):
    """Pick out only the sequence of targeted chromosome"""
    with FastaFile(reference) as reference_fasta:
        if chromosome in reference_fasta:
            return reference_fasta[chromosome][minpos:maxpos]
        else:
            return None


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


def find_minmax_pos(bam, chromosome, flags, flag_filter, min_quality, desc="find_minmax_pos"):
    """Determine leftmost and rightmost mapping positions in given chunk of the tailpuller file"""
    minpos, maxpos = float("inf"), 0
    with AlignmentFile(bam) as alignment:
        entry_iterator = filter_bam(
            alignment.fetch(chromosome), flags, flag_filter, min_quality,
            desc=desc
        )
        for entry in entry_iterator:
            minpos = min(minpos, entry.reference_start)
            maxpos = max(maxpos, entry.reference_end)
    if minpos < maxpos:
        return minpos, maxpos
    else:
        return None, None


def assemble_around_chromosome(bam, chromosome, reference, kmer_size, flags, flag_filter, min_quality, output_prefix):
    """Perform local reference-guided assembly on one chromosome"""
    minpos, maxpos = find_minmax_pos(
        bam, chromosome, flags, flag_filter, min_quality,
        desc="determining min/max positions on {}".format(chromosome)
    )
    if minpos and maxpos:
        reference_seq = get_reference_seq(reference, chromosome, minpos, maxpos)
        if reference_seq:
            reference_kmers, ordered_reference_kmers = count_kmers(
                reference_seq, kmer_size, with_ordered=True,
                desc="Counting kmers in {}".format(chromosome)
            )


def main(bam, reference, flags, flag_filter, min_quality, kmer_size, chromosomes, output_prefix, jobs=1, file=stdout, **kwargs):
    # parse and check arguments:
    if chromosomes:
        target_chromosomes = natsorted_chromosomes(chromosomes.split("|"))
    else:
        with FastaFile(reference) as fa:
            target_chromosomes = natsorted_chromosomes(fa.references)
    for chromosome in target_chromosomes:
        assemble_around_chromosome(
            bam, chromosome, reference, kmer_size,
            interpret_flags(flags), interpret_flags(flag_filter), min_quality,
            output_prefix
        )

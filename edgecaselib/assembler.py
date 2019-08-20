from sys import stdout
from edgecaselib.util import natsorted_chromosomes
from edgecaselib.formats import filter_bam
from pysam import FastaFile, AlignmentFile
from tqdm import tqdm
from pickle import dump


def get_reference_seq(reference, chromosome, minpos=None, maxpos=None):
    """Pick out only the sequence of targeted chromosome"""
    with FastaFile(reference) as reference_fasta:
        if chromosome in reference_fasta:
            return reference_fasta[chromosome][minpos:maxpos]
        else:
            return None


def collect_unique_kmers(sequence, k, desc=None, dinvert=False):
    """Collect identities and positions of unique k-mers in a sequence"""
    kmer_table, repeated_kmers = {}, set()
    if desc:
        sequence_iterator = tqdm(sequence[k:], desc=desc)
    else:
        sequence_iterator = iter(sequence[k:])
    kmer = sequence[:k].upper()
    kmer_table[kmer] = 0
    for pos, base in enumerate(sequence_iterator, start=1):
        kmer = kmer[1:] + base.upper()
        if kmer not in repeated_kmers:
            if kmer in kmer_table:
                repeated_kmers.add(kmer)
                del kmer_table[kmer]
            else:
                kmer_table[kmer] = pos
    if dinvert:
        return {p: m for m, p in kmer_table.items()}
    else:
        return kmer_table


def find_minmax_pos(bam, chromosome, samfilters, desc="find_minmax_pos"):
    """Determine leftmost and rightmost mapping positions in given chunk of the tailpuller file"""
    minpos, maxpos = float("inf"), 0
    with AlignmentFile(bam) as alignment:
        entry_iterator = filter_bam(
            alignment.fetch(chromosome), samfilters, desc=desc
        )
        for entry in entry_iterator:
            minpos = min(minpos, entry.reference_start)
            maxpos = max(maxpos, entry.reference_end)
    if minpos < maxpos:
        return minpos, maxpos
    else:
        return None, None


def assemble_around_chromosome(bam, chromosome, reference, reference_guided, kmer_size, samfilters, output_prefix):
    """Perform local reference-guided assembly on one chromosome"""
    pickle_prefix = "{}-{}-{:02d}k".format(output_prefix, chromosome, kmer_size)
    minpos, maxpos = find_minmax_pos(
        bam, chromosome, samfilters,
        desc="determining min/max positions on {}".format(chromosome)
    )
    if minpos and maxpos:
        reference_seq = get_reference_seq(reference, chromosome, minpos, maxpos)
        if reference_guided and reference_seq:
            reference_kmers = collect_unique_kmers(
                reference_seq, kmer_size,
                desc="Counting kmers in {}".format(chromosome)
            )
            with open(pickle_prefix + "-ref_kmers.pkl", mode="wb") as pkl:
                dump(reference_kmers, pkl)
        read_kmers = {}
        with AlignmentFile(bam) as alignment:
            decorated_entry_iterator = filter_bam(
                alignment.fetch(chromosome), samfilters,
                desc="Counting kmers in reads mapped to {}".format(chromosome)
            )
            for entry in decorated_entry_iterator:
                read_kmers[entry.qname] = collect_unique_kmers(
                    entry.seq, kmer_size
                )
        with open(pickle_prefix + "-read_kmers.pkl", mode="wb") as pkl:
            dump(read_kmers, pkl)


def main(bam, reference, reference_guided, flags, flags_any, flag_filter, min_quality, kmer_size, chromosomes, output_prefix, jobs=1, file=stdout, **kwargs):
    # parse and check arguments:
    if chromosomes:
        target_chromosomes = natsorted_chromosomes(chromosomes.split("|"))
    else:
        with FastaFile(reference) as fa:
            target_chromosomes = natsorted_chromosomes(fa.references)
    for chromosome in target_chromosomes:
        assemble_around_chromosome(
            bam, chromosome, reference, reference_guided, kmer_size,
            [flags, flags_any, flag_filter, min_quality],
            output_prefix
        )

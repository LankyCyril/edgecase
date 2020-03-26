from sys import stdout, stderr
from collections import OrderedDict
from edgecaselib.util import get_executable, progressbar
from edgecaselib.kmerscanner import get_circular_pattern
from pysam import FastxFile
from numpy import array
from edgecaselib.repeatfinder import lowest_alpha_inversion
from functools import lru_cache
from tempfile import TemporaryDirectory
from os import path
from subprocess import check_output
from pandas import DataFrame

__doc__ = """edgeCase shortread: experiments with short reads

Usage: {0} shortread [-m integer] [-M integer] [-r integer]
       {1}           [--kmer-counter string] [--motifs string]
       {1}           [--fmt string] [-n integer] <sequencefile>

Output:
    TSV-formatted file with motif incidence per read

Positional arguments:
    <sequencefile>                name of input BAM/SAM/FASTA/FASTQ file

Options:
    --fmt sam|fastx               format of input file [default: fastx]
    -m, --min-k [integer]         smallest target repeat length [default: 4]
    -M, --max-k [integer]         largest target repeat length [default: 16]
    -r, --min-repeats [integer]   minimum number of consecutive repeats [default: 2]
    --kmer-counter [string]       kmer-counter binary (unless in $PATH)
    --motifs [string]             list of target motifs, separated with '|'
    -n, --num-reads [integer]     expected number of reads in input (for progress display)

Notes:
* if --motifs is specified, options -m and -M have no effect
* kmer-counter is available at https://github.com/alexpreynolds/kmer-counter
"""

__docopt_converters__ = [
    lambda min_k: int(min_k),
    lambda max_k: int(max_k),
    lambda min_repeats: int(min_repeats),
    lambda num_reads: None if (num_reads is None) else int(num_reads),
]

__docopt_tests__ = {
    lambda min_k, max_k: 0 < min_k < max_k: "not satisfied: 0 < m < M",
    lambda min_repeats: min_repeats > 0: "--min-repeats must be integer > 0",
    lambda fmt: fmt in {"sam", "fastx"}: "unknown value of --fmt",
}


BUFFER_CHUNKSIZE = 65536


def interpret_arguments(fmt, min_k, max_k, min_repeats, kmer_counter, motifs):
    """Parse and check arguments"""
    if fmt == "fastx":
        manager, seq_attribute = FastxFile, "sequence"
    else:
        raise NotImplementedError("--fmt={}".format(fmt))
    if motifs:
        print("WARNING: Using --motifs, ignoring -m and -M", file=stderr)
        motif_patterns = OrderedDict([
            [motif, get_circular_pattern(motif)]
            for motif in motifs.split("|")
        ])
    else:
        motif_patterns = None
    return (
        manager, seq_attribute, motif_patterns,
        get_executable("kmer-counter", kmer_counter)
    )


def sweep_specific_motifs(entry_iterator, seq_attribute, motif_patterns, num_reads):
    """Count motif_patterns in each entry from entry_iterator"""
    decorated_iterator = progressbar(
        entry_iterator, total=num_reads, unit="read", desc="Sweeping motifs"
    )
    for entry in decorated_iterator:
        sequence = getattr(entry, seq_attribute, None)
        if sequence:
            motif_counts = array([
                sum(1 for _ in pattern.finditer(sequence, overlapped=True))
                for pattern in motif_patterns.values()
            ])
            if (motif_counts != 0).any():
                yield motif_counts


@lru_cache(maxsize=None)
def get_motif_identity(kmer, min_repeats):
    """Reduce kmer to motif if motif in repeat context"""
    lai = lowest_alpha_inversion(kmer)
    motif = lai[:int(len(lai)/min_repeats)]
    if motif * min_repeats == lai:
        return motif
    else:
        return None


def generate_motif_count_database(entry_iterator, seq_attribute, min_k, max_k, min_repeats, kmer_counter, num_reads):
    """Count incidence of all motifs from min_k to max_k in repeat contexts"""
    motif_count_database, sequence_buffer = [], []
    decorated_iterator = progressbar(
        entry_iterator, total=num_reads, unit="read", desc="Sweeping motifs"
    )
    with TemporaryDirectory() as tempdir:
        temp_fa = path.join(tempdir, "input.fa")
        for entry in decorated_iterator:
            sequence, counts = getattr(entry, seq_attribute, None), {}
            if sequence:
                sequence_buffer.append(sequence)
            if len(sequence_buffer) >= BUFFER_CHUNKSIZE:
                with open(temp_fa, mode="wt") as handle:
                    for i, sequence in enumerate(sequence_buffer):
                        print(">{}\n{}".format(i, sequence), file=handle)
                sequence_buffer = []
                for k in range(min_k, max_k+1):
                    counter_cmd = [
                        kmer_counter, "--double-count-palindromes",
                        "--k", str(k), "--fasta", temp_fa
                    ]
                    raw_counts = check_output(counter_cmd).decode()
                    for line in raw_counts.strip().split("\n"):
                        stats = line.split("\t")[1].split(" ")
                        for stat in stats:
                            try:
                                kmer, count = stat.split(":")
                            except ValueError: # some lines are empty
                                continue
                            motif = get_motif_identity(kmer, min_repeats)
                            if motif:
                                counts[motif] = int(count)
                if counts:
                    motif_count_database.append(counts)
    return motif_count_database


def main(sequencefile, fmt, min_k, max_k, min_repeats, motifs, kmer_counter, num_reads, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    manager, seq_attribute, motif_patterns, kmer_counter = interpret_arguments(
        fmt, min_k, max_k, min_repeats, kmer_counter, motifs,
    )
    with manager(sequencefile) as entry_iterator:
        if motif_patterns:
            counts_iterator = sweep_specific_motifs(
                entry_iterator, seq_attribute, motif_patterns, num_reads,
            )
            for counts in counts_iterator:
                print(*counts, sep="\t", file=file)
        else:
            motif_count_database = generate_motif_count_database(
                entry_iterator, seq_attribute, min_k, max_k, min_repeats,
                kmer_counter, num_reads,
            )
            print("Merging table...", file=stderr)
            motif_count_table = DataFrame(motif_count_database).fillna(0)
            print("Writing table...", file=stderr)
            motif_count_table.astype(int).to_csv(file, sep="\t", index=False)

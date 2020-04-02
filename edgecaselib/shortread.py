from sys import stdout, stderr
from re import compile
from contextlib import contextmanager, ExitStack
from edgecaselib.util import get_executable, progressbar
from edgecaselib.repeatfinder import lowest_alpha_inversion
from functools import lru_cache
from tempfile import TemporaryDirectory
from os import path
from pysam import FastxFile, AlignmentFile
from concurrent.futures import ThreadPoolExecutor
from subprocess import call
from collections import defaultdict
from gzip import open as gzopen
from pandas import DataFrame, Series, concat
from scipy.stats import spearmanr
from numpy import nan
from statsmodels.stats.multitest import multipletests


__doc__ = """edgeCase shortread: experiments with short reads

Usage: {0} shortread [-j integer] [-c integer] [--kmer-counter string]
       {1}           [--target string] [-m integer] [-M integer] [-r integer]
       {1}           [--prefix string] [--fmt string] <sequencefile>

Output:
    TSV-formatted file with motif correlations

Positional arguments:
    <sequencefile>                   name of input BAM/SAM/FASTA/FASTQ file

Options:
    --fmt sam|fastx                  format of input file [default: fastx]
    -j, --jobs [integer]             number of jobs to run in parallel [default: 1]
    -c, --chunks-per-job [integer]   number of chunks of <sequencefile> per job [default: 32]
    --kmer-counter [string]          kmer-counter binary (unless in $PATH)
    --target [string]                only consider reads containing this kmer [default: TTAGGG]
    -m, --min-k [integer]            smallest target repeat length [default: 4]
    -M, --max-k [integer]            largest target repeat length [default: 16]
    -r, --min-repeats [integer]      minimum number of consecutive repeats [default: 2]
    --prefix [string]                prefix for temporary files

Notes:
* reverse-complement of the target motif will also be considered
* kmer-counter is available at https://github.com/alexpreynolds/kmer-counter
* if --prefix is not specified, will store temporary files in $TEMP
"""

__docopt_converters__ = [
    lambda jobs: int(jobs),
    lambda min_k: int(min_k),
    lambda max_k: int(max_k),
    lambda min_repeats: int(min_repeats),
    lambda chunks_per_job: int(chunks_per_job)
]

__docopt_tests__ = {
    lambda min_k, max_k:
        0 < min_k < max_k:
            "not satisfied: 0 < m < M",
    lambda min_repeats:
        min_repeats > 0:
            "--min-repeats must be integer > 0",
    lambda fmt:
        fmt in {"sam", "fastx"}:
            "unknown value of --fmt",
    lambda target, min_k, max_k:
        min_k <= len(target) <= max_k:
            "length of --target is outside of --min-k and --max-k boundaries",
}


ALPHABET = list("AaCcNnnNgGtT")
COMPLEMENTS = dict(zip(ALPHABET, reversed(ALPHABET)))
COMPLEMENT_PATTERN = compile(r'|'.join(COMPLEMENTS.keys()))


@lru_cache(maxsize=None)
def revcomp(sequence):
    matcher = lambda match: COMPLEMENTS[match.group()]
    return COMPLEMENT_PATTERN.sub(matcher, sequence[::-1])


@lru_cache(maxsize=None)
def all_inversions(kmer):
    return {kmer[i:]+kmer[:i] for i in range(len(kmer))}


@lru_cache(maxsize=None)
def is_revcomp_inversion(s1, s2):
    return len(all_inversions(s1) & all_inversions(revcomp(s2))) != 0


@lru_cache(maxsize=None)
def get_motif_identity(kmer, min_repeats=2):
    """Reduce kmer to motif if motif in repeat context"""
    lai = lowest_alpha_inversion(kmer)
    motif = lai[:int(len(lai)/min_repeats)]
    if motif * min_repeats == lai:
        return motif
    else:
        return None


def interpret_arguments(kmer_counter, fmt, prefix, chunks_per_job, jobs):
    """Parse and check arguments"""
    if fmt == "fastx":
        manager, seq_attr = FastxFile, "sequence"
    elif fmt == "sam":
        manager, seq_attr = AlignmentFile, "query_sequence"
    else:
        raise ValueError("Unknown --fmt: {}".format(fmt))
    if prefix is None:
        TempPrefix = TemporaryDirectory
    else:
        def TempPrefix_inner():
            yield path.abspath(prefix)
        TempPrefix = contextmanager(TempPrefix_inner)
    return (
        get_executable("kmer-counter", kmer_counter),
        manager, seq_attr, TempPrefix, chunks_per_job * jobs,
    )


def generate_temp_name(temp_prefix, name):
    """Generate temporary name based on whether temp_prefix is a directory"""
    if path.isdir(temp_prefix):
        return path.join(temp_prefix, name)
    else:
        return temp_prefix + name


def preprocess_input(sequencefile, manager, seq_attr, chunks, temp_prefix):
    """Split input into many small FASTA files"""
    chunked_fastas = [
        generate_temp_name(temp_prefix, str(i)) for i in range(chunks)
    ]
    with FastxFile(sequencefile) as fasta, ExitStack() as chunk_stack:
        chunked_fasta_handles = [
            chunk_stack.enter_context(open(chunkname, mode="wt"))
            for chunkname in chunked_fastas
        ]
        i = 0
        for entry in progressbar(fasta, desc="Splitting input", unit="entry"):
            print(entry, file=chunked_fasta_handles[i])
            i = (i + 1) % chunks
    return chunked_fastas


def get_chunk_child_name(chunkname, *args):
    """Facilitate reproducible naming of temporary files"""
    return chunkname + "-" + "-".join([str(a) for a in args])


def kmer_counter_pipeline(kmer_counter, k, chunkname, chunk_child_name):
    """Unsafely pipe data through standard GNU/Linux tools"""
    pipeline_mask = " ".join([
        "{} --double-count-palindromes --k {} --fasta {}",
        "| cut -f2 | sed 's/$/ ;/g' | tr ' ' '\\n'",
        "| awk -F':' '{{if (($1==\";\") || (substr($1,1,{})==substr($1,{})))",
            "{{print}}}}'",
        "| gzip -3 > {}",
    ])
    cmd = pipeline_mask.format(
        kmer_counter, k * 2, chunkname, k, k + 1, chunk_child_name,
    )
    call(cmd, shell=True)


def prepare_kmer_counter_kwargs(kmer_counter, min_k, max_k, min_repeats, chunked_fastas):
    """Stage thread arguments for kmer-counter"""
    for chunkname in chunked_fastas:
        for k in range(min_k, max_k + 1):
            yield {
                "kmer_counter": kmer_counter,
                "k": k,
                "chunkname": chunkname,
                "chunk_child_name": get_chunk_child_name(
                    chunkname, k, "doubled",
                ),
            }


def run_kmer_counter(kmer_counter, min_k, max_k, min_repeats, chunked_fastas, jobs):
    """Run kmer-counter in multiple threads"""
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        thread_kwargs = list(prepare_kmer_counter_kwargs(
            kmer_counter, min_k, max_k, min_repeats, chunked_fastas,
        ))
        futures = [
            pool.submit(kmer_counter_pipeline, **kwargs)
            for kwargs in thread_kwargs
        ]
        _ = [
            future.result() for future in
            progressbar(futures, desc="Counting kmers", unit="chunk")
        ]
    return [
        kwargs["chunk_child_name"] for kwargs in thread_kwargs
    ]


def generate_count_table(chunked_fastas, k):
    """Combine outputs of kmer-counter for one k into a DataFrame"""
    motif_count_database, counts = [], defaultdict(int)
    for chunkname in chunked_fastas:
        chunk_child_name = get_chunk_child_name(chunkname, k, "doubled")
        with gzopen(chunk_child_name, mode="rt") as raw_kc:
            for line in map(str.strip, raw_kc):
                if line == ";":
                    motif_count_database.append(counts)
                    counts = defaultdict(int)
                else:
                    stat = line.split(":")
                    if len(stat) == 2:
                        counts[get_motif_identity(stat[0])] += int(stat[1])
    return DataFrame(motif_count_database)


def correlate_motifs(count_tables, target, method=spearmanr, adjust="bonferroni", alpha=.05):
    """Combine motif counts and correlate; report significant correlations only"""
    print("Merging...", end=" ", file=stderr, flush=True)
    count_table = concat(count_tables, axis=1)
    print("done", file=stderr, flush=True)
    print("Correlating...", end=" ", file=stderr, flush=True)
    target_lai = lowest_alpha_inversion(target)
    correlation_matrix = count_table.drop(columns=target_lai).corrwith(
        count_table[target_lai],
        method=lambda x, y: Series(method(x, y)),
    ).T
    correlation_matrix.columns = [target, "p"]
    # remove revcomp-to-revcomp pairs and negavite correlations:
    for motif, (r, _) in correlation_matrix.iterrows():
        if (r < 0) or is_revcomp_inversion(target_lai, motif):
            correlation_matrix.loc[motif] = nan
    # apply multiple testing correction:
    correlation_matrix = correlation_matrix.dropna()
    correlation_matrix["p_adjusted"] = multipletests(
        correlation_matrix["p"], method=adjust,
    )[1]
    correlation_matrix["pass"] = (correlation_matrix["p_adjusted"] < alpha)
    correlation_matrix = correlation_matrix.sort_values(
        by=["pass", target], ascending=False,
    )
    # drop rows and columns not containing any significant correlations:
    print("done", file=stderr, flush=True)
    return correlation_matrix.drop(columns="pass")


def main(sequencefile, fmt, target, min_k, max_k, min_repeats, kmer_counter, prefix, chunks_per_job, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    if min_repeats != 2:
        raise NotImplementedError("Only --min-repeats=2 is implemented so far")
    kmer_counter, manager, seq_attr, TempPrefix, chunks = interpret_arguments(
        kmer_counter, fmt, prefix, chunks_per_job, jobs,
    )
    with TempPrefix() as temp_prefix:
        chunked_fastas = preprocess_input(
            sequencefile, manager, seq_attr, chunks, temp_prefix,
        )
        kmer_count_files = run_kmer_counter(
            kmer_counter, min_k, max_k, min_repeats, chunked_fastas, jobs,
        )
        count_tables = [
            generate_count_table(chunked_fastas, k)
            for k in progressbar(
                range(min_k, max_k + 1), unit="k",
                desc="Interpreting motif counts",
            )
        ]
        correlation_matrix = correlate_motifs(
            count_tables, target=target, method=spearmanr,
        )
    correlation_matrix.to_csv(file, sep="\t")

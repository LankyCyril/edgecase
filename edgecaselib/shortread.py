from sys import stdout
from contextlib import contextmanager, ExitStack
from edgecaselib.util import get_executable, progressbar
from edgecaselib.repeatfinder import lowest_alpha_inversion
from functools import lru_cache
from tempfile import TemporaryDirectory
from os import path
from pysam import FastxFile
from threading import Thread
from subprocess import call

__doc__ = """edgeCase shortread: experiments with short reads

Usage: {0} shortread [-j integer] [--bioawk string] [--kmer-counter string]
       {1}           [--target string] [-m integer] [-M integer] [-r integer]
       {1}           [--prefix string] [--chunks integer]
       {1}           [--fmt string] <sequencefile>

Output:
    TSV-formatted file with motif incidence per read

Positional arguments:
    <sequencefile>                name of input BAM/SAM/FASTA/FASTQ file

Options:
    --fmt sam|fastx               format of input file [default: fastx]
    -j, --jobs [integer]          number of jobs to run in parallel [default: 1]
    -c, --chunks [integer]        number of chunks in which to split <sequencefile>
    --kmer-counter [string]       kmer-counter binary (unless in $PATH)
    --bioawk [string]             bioawk binary (unless in $PATH)
    --target [string]             only consider reads containing this kmer [default: TTAGGG]
    -m, --min-k [integer]         smallest target repeat length [default: 4]
    -M, --max-k [integer]         largest target repeat length [default: 16]
    -r, --min-repeats [integer]   minimum number of consecutive repeats [default: 2]
    --prefix [string]             prefix for temporary files

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
    lambda chunks: None if chunks is None else int(chunks),
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


@lru_cache(maxsize=None)
def get_motif_identity(kmer, min_repeats):
    """Reduce kmer to motif if motif in repeat context"""
    lai = lowest_alpha_inversion(kmer)
    motif = lai[:int(len(lai)/min_repeats)]
    if motif * min_repeats == lai:
        return motif
    else:
        return None


def interpret_arguments(kmer_counter, bioawk, prefix, chunks, jobs):
    """Parse and check arguments"""
    if prefix is None:
        TempPrefix = TemporaryDirectory
    else:
        def TempPrefix_inner():
            yield path.abspath(prefix)
        TempPrefix = contextmanager(TempPrefix_inner)
    return (
        get_executable("kmer-counter", kmer_counter),
        get_executable("bioawk", bioawk),
        TempPrefix,
        chunks or jobs
    )


def generate_temp_name(temp_prefix, name):
    """Generate temporary name based on whether temp_prefix is a directory"""
    if path.isdir(temp_prefix):
        return path.join(temp_prefix, name)
    else:
        return temp_prefix + name


def chunk_input(sequencefile, chunks, temp_prefix):
    """Split input into many small FASTA files"""
    chunked_fastas = [
        generate_temp_name(temp_prefix, str(i))
        for i in range(chunks)
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


def kmer_counter_pipeline(kmer_counter, k, chunkname):
    pipeline_mask = " ".join([
        "{} --double-count-palindromes --k {} --fasta {}",
        "| cut -f2 | sed 's/$/ ;/g' | tr ' ' '\\n'",
        "| awk -F':' '{{if (($1==\";\") || (substr($1,1,{})==substr($1,{})))",
            "{{print}}}}'",
        "| gzip -3 > {}.txt.gz"
    ])
    cmd = pipeline_mask.format(
        kmer_counter, k * 2, chunkname, k, k + 1, chunkname
    )
    call(cmd, shell=True)


def main(sequencefile, fmt, target, min_k, max_k, min_repeats, kmer_counter, bioawk, prefix, chunks, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    kmer_counter, bioawk, TempPrefix, chunks = interpret_arguments(
        kmer_counter, bioawk, prefix, chunks, jobs
    )
    with TempPrefix() as temp_prefix:
        chunked_fastas = chunk_input(sequencefile, chunks, temp_prefix)
        threads = [
            Thread(
                target=kmer_counter_pipeline,
                args=(kmer_counter, 7, chunkname)
            )
            for chunkname in chunked_fastas
        ]
        job_offset_iterator = progressbar(
            range(0, len(threads), jobs), desc="Counting kmers", unit="batch"
        )
        for job_offset in job_offset_iterator:
            for thread in threads[job_offset:job_offset+jobs]:
                thread.start()
            for thread in threads[job_offset:job_offset+jobs]:
                thread.join()

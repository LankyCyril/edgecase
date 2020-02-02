#!/usr/bin/env python
from sys import argv, stderr
from re import compile
from tempfile import TemporaryDirectory
from os import path
from urllib.request import urlretrieve
from pysam import FastxFile
from tqdm import tqdm
from itertools import chain
from textwrap import fill


USAGE = """usage:
{0} --local hg38.fasta stong2014.fasta
    generate hg38ext.fa from local files
{0} --remote
    download appropriate assemblies and generate hg38ext.fa
"""

NCBI_FTP_DIR = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405"
HG38_RELEASE = "GCA_000001405.28_GRCh38.p13"
HG38_VERSION = "GRCh38_major_release_seqs_for_alignment_pipelines"
HG38_FASTAGZ = "GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"

CSHLP_URLDIR = "https://genome.cshlp.org/content/suppl/2014/04/16"
CSHLP_DOIDIR = "gr.166983.113.DC1"
CSHLP_SUPPFN = "Supplemental_FileS1.txt"

HG38_URL = "/".join([NCBI_FTP_DIR, HG38_RELEASE, HG38_VERSION, HG38_FASTAGZ])
STONG2014_URL = "/".join([CSHLP_URLDIR, CSHLP_DOIDIR, CSHLP_SUPPFN])

STONG_INFIX = "_1-500K_1_12_12"
STONG2014_SUBHAP_MASKS = {
    "1qtel{}", "2qtel{}", "3ptel{}", "4ptel{}", "5qtel{}", "6qtel{}", "9qtel{}",
    "10qtel{}", "14qtel{}", "16qtel{}", "17ptel{}", "17qtel{}v2", "18qtel{}",
    "19ptel{}",
}

REVCOMP_POSTFIX = "_rc"
TO_REVCOMP_MASKS = {
    "1qtel{}", "2qtel{}", "5qtel{}", "6qtel{}", "9qtel{}", "10qtel{}",
    "14qtel{}", "16qtel{}", "17qtel{}v2", "18qtel{}", "chrUn_KI270745v1",
}

ALPHABET = list("AaCcNnnNgGtT")
COMPLEMENTS = dict(zip(ALPHABET, reversed(ALPHABET)))
COMPLEMENT_PATTERN = compile(r'|'.join(COMPLEMENTS.keys()))


def revcomp(sequence):
    """Reverse-complement a sequence"""
    matcher = lambda match: COMPLEMENTS[match.group()]
    return COMPLEMENT_PATTERN.sub(matcher, sequence[::-1])


def retrieve_hook(count, blocksize, totalsize):
    """Progress bar hook for urlretrieve"""
    if totalsize > 0:
        percent = min(count * blocksize * 100 / totalsize, 100)
        print("{:.2f}% downloaded".format(percent), end="\r", file=stderr)
    else:
        megabytes = count * blocksize / 1e6
        print("{:.3f}Mb downloaded".format(megabytes), end="\r", file=stderr)


def download_assemblies(workdir):
    """Download hg38 and stong2014 assemblies from appropriate URLs"""
    for url, fn in [(HG38_URL, "hg38.fa.gz"), (STONG2014_URL, "stong2014.fa")]:
        local_fn = path.join(workdir, fn)
        print("Downloading into {} from {}".format(local_fn, url), file=stderr)
        urlretrieve(url=url, filename=local_fn, reporthook=retrieve_hook)
        yield local_fn


def parser_iterator(filename, to_revcomp, desc="Parsing", entry_filter=lambda e:True):
    """Parse FASTA and iterate over converted entries"""
    with FastxFile(filename) as fasta:
        for entry in tqdm(fasta, desc=desc, unit="contig"):
            if entry_filter(entry):
                if entry.name in to_revcomp:
                    name = entry.name + REVCOMP_POSTFIX
                    sequence = revcomp(entry.sequence)
                else:
                    name, sequence = entry.name, entry.sequence
                yield name, fill(sequence)


def generate_hg38ext(hg38, stong2014):
    """Generate hg38ext from the hg38 and stong2014 FASTA files"""
    subhaps = {mask.format(STONG_INFIX) for mask in STONG2014_SUBHAP_MASKS}
    to_revcomp = {mask.format(STONG_INFIX) for mask in TO_REVCOMP_MASKS}
    hg38_iterator = parser_iterator(hg38, to_revcomp, desc="Parsing reference")
    stong2014_iterator = parser_iterator(
        stong2014, to_revcomp, desc="Parsing subtelomeres",
        entry_filter=lambda e: (e.name in subhaps)
    )
    for name, sequence_lines in chain(hg38_iterator, stong2014_iterator):
        print(">" + name)
        print(sequence_lines)
    return 0


if __name__ == "__main__":
    # interpret command-line arguments and dispatch to subroutines:
    if (len(argv) == 2) and (argv[1] == "--remote"):
        with TemporaryDirectory() as workdir:
            hg38, stong2014 = download_assemblies(workdir)
            returncode = generate_hg38ext(hg38, stong2014)
    elif (len(argv) == 4) and (argv[1] == "--local"):
        returncode = generate_hg38ext(hg38=argv[2], stong2014=argv[3])
    else:
        print(USAGE.format(__file__).rstrip(), file=stderr)
        returncode = 1
    exit(returncode)

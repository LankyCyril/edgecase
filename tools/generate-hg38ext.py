#!/usr/bin/env python
from sys import argv, stderr
from re import compile
from tempfile import TemporaryDirectory
from os import path
from urllib.request import urlretrieve
from contextlib import contextmanager
from binascii import hexlify
from gzip import open as gzopen
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

STONG_INFIX = "500K_1_12_12"
STONG2014_SUBHAP_MASKS = {
    "1qtel_1-{}", "2qtel_1-{}", "3ptel_1-{}", "4ptel_1-{}", "5qtel_1-{}",
    "6qtel_1-{}", "9qtel_1-{}", "10qtel_1-{}", "14qtel_1-{}", "16qtel_1-{}",
    "17ptel_1_{}", "17qtel_1-{}v2", "18qtel_1-{}", "19ptel_1-{}"
}

REVCOMP_POSTFIX = "_rc"
TO_REVCOMP_MASKS = {
    "1qtel_1-{}", "2qtel_1-{}", "5qtel_1-{}", "6qtel_1-{}", "9qtel_1-{}",
    "10qtel_1-{}", "14qtel_1-{}", "16qtel_1-{}", "17qtel_1-{}v2", "18qtel_1-{}",
    "chrUn_KI270745v1"
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


@contextmanager
def open_fasta(filename):
    """Open FASTA with 'open' if plaintext, 'gzip.open' if gzipped"""
    with open(filename, mode="rb") as bytes_handle:
        is_gzipped = (hexlify(bytes_handle.read(2)) == b"1f8b")
    if is_gzipped:
        yield gzopen(filename, mode="rt")
    else:
        yield open(filename, mode="rt")


def parser_iterator(filename, to_revcomp, desc="Parsing", entry_filter=lambda e:True):
    """Parse FASTA and iterate over converted entries"""
    entry_accepted, entry_needs_revcomp, lines_to_revcomp = False, False, []
    with open_fasta(filename) as fasta:
        for line in tqdm(map(str.strip, fasta), desc=desc, unit="line"):
            if line.startswith(">"):
                if lines_to_revcomp: # generated from previous entry
                    previous_entry_sequence = revcomp("".join(lines_to_revcomp))
                    yield fill(previous_entry_sequence)
                    lines_to_revcomp = []
                name = line[1:].split()[0]
                entry_accepted = entry_filter(name)
                entry_needs_revcomp = (name in to_revcomp)
                if entry_accepted:
                    if entry_needs_revcomp:
                        yield ">" + name + REVCOMP_POSTFIX
                    else:
                        yield ">" + name
            elif entry_accepted:
                if entry_needs_revcomp:
                    lines_to_revcomp.append(line)
                else:
                    yield line


def generate_hg38ext(hg38, stong2014):
    """Generate hg38ext from the hg38 and stong2014 FASTA files"""
    subhaps = {mask.format(STONG_INFIX) for mask in STONG2014_SUBHAP_MASKS}
    to_revcomp = {mask.format(STONG_INFIX) for mask in TO_REVCOMP_MASKS}
    hg38_iterator = parser_iterator(hg38, to_revcomp, desc="Parsing reference")
    stong2014_iterator = parser_iterator(
        stong2014, to_revcomp, desc="Parsing subtelomeres",
        entry_filter=lambda name: (name in subhaps)
    )
    for line in chain(hg38_iterator, stong2014_iterator):
        print(line)
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

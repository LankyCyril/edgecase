#!/usr/bin/env python
from sys import argv, stderr
from re import compile
from tempfile import TemporaryDirectory
from os import path
from urllib.request import urlretrieve
from contextlib import contextmanager
from binascii import hexlify
from gzip import open as gzopen
from itertools import chain
from textwrap import fill
from zlib import decompress
from base64 import decodebytes

try:
    from tqdm import tqdm
except ModuleNotFoundError:
    def tqdm(it, desc=None, *args, **kwargs):
        if desc is not None:
            print(desc, end="...\n", file=stderr)
        return it


USAGE = """usage:
{0} --local hg38.fasta stong2014.fasta
    generate hg38ext.fa from local files and output to stdout
{0} --remote
    download appropriate assemblies, generate hg38ext.fa, and output to stdout
{0} --ecx
    output the edgeCase indeX (hg38ext.fa.ecx) to stdout

NOTE! This tool writes uncompressed data (FASTA or ECX) to stdout.
You should pipe it into a file, for example:
{0} --remote > hg38ext.fa
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

HG38EXT_ECX_GZ = (
b'eNqVWEuPGzcMPs/+jR6LBUQ9yQZBj0WRXgskJ8N1N9nFPms7AfrvS+oxOw+L4wYIRXjJT6LE+Ujp'
b'\np69P+2+nX573p8fd/uVw/3r86A3FD19fj48fEch+OB/3h3P7G0SH/uanu5fz8d/h+LJ/vhveXk/'
b'y\n/2cYDvfH1+fX0yv/+Lx/eNmVv8sEw9vxgdXD0/50Gk6v34+Hu+Hp4eVx+Otpf3h8ejidb4y4ww'
b'CG\n/2UJ5YcsZElDGCarHO6/ORxuh9ubbGKXjrYIzTGbuKWjK0JzzCZ+6eiL0ByzSVg6hiI0x2wSh'
b'5gd\nY3OMRWiO2SQtZ0xFaI7ZBJczYhGaYzah5YxUhOaYTcCsMsBUqeZAMYLleqFJ1bkYrVIIbJWq'
b'czHi\nNMpTi39RoP5epApRjHwHwlepQhQjTq3UINIEIlSpQhSjuNqFWKXqXIzSav9TlapzMcLVzFi'
b'l6pyN\nrFnObE2VqnPOOQtDMGXuMkL9tUiVQsoskjqhrT68Q9gqVYhs+nkZ/OciVM8M/2Xp+aUI1b'
b'Nk3GA9\nko+egZrmlqTrOgi+8Kv1FtAFS01zZsm+PYRQ6ZfQuhACNS2aJQ33EGLlYTLWeIYYtbjk4'
b'x5CJUII\n5Fyg5EbNL2myh1AZEXwwCaOLo5aWfNlDqNQIDp3ziT+GpuGSODsIrjEkOJcwyXE2za0Z'
b'tAfSmBIc\nrz/GDFI0t2bSHkhjTJ7fRs4Falpdw4xReyCNMwG8C97xFjeN1pzaA2msCSYiupTZpWi'
b'0ZtUeSONN\nMJybAEhNI7Pm1R5Io0Z01ifvoSl2TZA9iEaQaGzk/QhNiWua7EEUIxoCRj5U5uuqpP'
b'p7kSoENab1\nnJ4QU1NwzbcdCG8ap8ZIROia4teM24NozBoMGvQ8eVVozbg9CFsZl6nGGaQwanHJv'
b'D0EVzg2JAvJ\ncwJUJS75t+dfso/PIBBiGWjCvLnLZvI+Ptyd75/3L7vn/flw3/z//v72640Pre21'
b'sUg3IV4dgBcQ\naycKaFKWOOlEN72T1o5uemPrSa3LKy8ffKPiLW/SGtMt72Bad8qVKks/IdlNbxh'
b'bVA8C4q2ZEewm\ngN1oUzcB3NiqhhSyjDNi3QTw712ezfeCUsHGLm8TIIzNWszMGGHOQpsAsRJFIi'
b'60IqOZEcgVyR/S\nvHfASRfx3jsUINdfCY79A1gb0I2an/QPOkpZDk3bCBzbCJy1EVvLKdtQWwknm'
b'dm0sMpwBWVex0Ok\npiVzIdUUoFkt98Sc1bRwIeUUIDetogaBmmbnVXQTyE9rmKfYlHkN24QJYx3j'
b'CpqY+apCF9JQPfcY\np+VMTrsqYVbONheUJiUtST2qSpyVtE0YfC9rlusJjJqdlLVNlLKNdvfbH5i'
b'4hww/YLd/Og+eAvO0\nSD87fXki4o9VnosKCKfZ7U1q7R5/ktyeigTdjedlt3rj3X363SYp6nVy5G'
b'+Sbxx5CNP2bWUosG66\nGiuw9r0Z5K9SSLdpdtYMrt1B3HMC//lSJktetuR44N31DCzSzTq59Y54A'
b'fEbN+S1mxO3MLi3893T\nDm6DMZ92sOOj4f3j/FwxXRjEMAO007WGo8rPOfEyDm9nuUZVpLqKy0Ap'
b'CVLjXojS3kMEo/iXSKL4\ncYOqRDJjWzUSuoyTAmLKEid4ajwoqYqmWnI5c1l6xb/EI7mKnESpLGR'
b'3KaDxYLdDQtuB4gtUlRdT\nrBOV5A02B/DghDxk0EFyaCjfC0pNVM5qQbVqaKEDlfgqJgJngHpgkk'
b'b1zkOr76gHUcIK4sq78c96\nKfIte5J/zLsZdNmIOy1A7IJaoiifmAzTd5XKMp0gc27R+D6TAiSyo'
b'6bhlEhRnlLNYK+NdHZj0CLl\ngtIFjUh8nHnACagaKUmqkX1/RwocZBy1pODkSEm4mTeEro101uWr'
b'kfo+qHNyW81DmIDqkUr2URjf\nediW5NJWNafglEilgJD07Fen7/xGogabFFx5/oEy2BmuHm/OQXx'
b'/lYq5lRg1VKFKyFJqOCTwV4fc\naQm2dwDk9bQ3jyOTuEjmIV7RelyehhvePFGnvwHLzhx1Gen/TJ'
b'N3SyIQeKkjF+L40dmxeYXSt8hp\n0HyZBkGQIayLTH9XfMb24wuZIeOkOyyKU6Fq5C5DcJnBqzNlf'
b'lXV444aMN+pQhnmd+CNqFNGnjzq\n8d2Tl1YVUKFq1DFD4BCuDXr2tKOHTApqknPOQ5ig6vFCTk6o'
b'ywAEn6STHzVSkGq4lBFgiNeGO3uL\nUsPldOuj8l9tGdwEdSPcnJFQHSDfGIi736aRglTCBZsRpBu'
b'+Oqfn7zd6xEEBRjJchPOQZsAbQeeE\nHI3JWIvO+6YEFapGLVXyP0oAIII=\n')


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
    """Generate hg38ext from the hg38 and stong2014 FASTA files, write to stdout"""
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


def output_ecx():
    """Write the hg38ext ECX to stdout"""
    print(decompress(decodebytes(HG38EXT_ECX_GZ)).decode().rstrip("\n"))


if __name__ == "__main__":
    # interpret command-line arguments and dispatch to subroutines:
    if (len(argv) == 2) and (argv[1] == "--remote"):
        with TemporaryDirectory() as workdir:
            hg38, stong2014 = download_assemblies(workdir)
            returncode = generate_hg38ext(hg38, stong2014)
    elif (len(argv) == 4) and (argv[1] == "--local"):
        returncode = generate_hg38ext(hg38=argv[2], stong2014=argv[3])
    elif (len(argv) == 2) and (argv[1] == "--ecx"):
        returncode = output_ecx()
    else:
        print(USAGE.format(__file__).rstrip(), file=stderr)
        returncode = 1
    exit(returncode)

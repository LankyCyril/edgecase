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

HG38EXT_ECX_GZ = b'\n'.join([
b'eJyVWEtvIzcMPk/+Ro9FAFFPsos9F8X2WmD3ZLhGdhMkcba2u2j/fUk9JvOwOG6AUMSY/CRqOB8p'
b'/fT1Zf/t/Mvr/vy82x8Pj2+nj95Q/PD17fT8EYHsh8tpf7i03yA69Hc/PRwvp3+H03H/+jB8fzvL'
b'/88wHB5Pb69v5zd++Lp/Ou7K7zLB8P30xOrhZX8+D+e3v0+Hh+Hl6fg8/PmyPzy/PJ0vd0bcYQDD'
b'f1lCeZCFLGkIw2SVw+M3h8P9cH+XTezS0RahOWYTt3R0RWiO2cQvHX0RmmM2CUvHUITmmE3iELNj'
b'bI6xCM0xm6TljKkIzTGb4HJGLEJzzCa0nJGK0ByzCZhVBpgq1RwoRrBcLzSpOhejVQqBrVJ1Lkac'
b'Rnlq8S8K1OdFqhDFyHcgfJUqRDHi1EoNIk0gQpUqRDGKq12IVarOxSit9j9VqToXI1zNjFWqztnI'
b'muXM1lSpOuecszAEU+YuI9SnRaoUUmaR1Alt9eEdwlapQmTTz8vgPxehemb4L0vPL0WoniXjBuuR'
b'fPQM1DS3JF3XQfCFX623gC5YapozS/btIYRKv4TWhRCoadEsabiHECsPk7HGM8SoxSUf9xAqEUIg'
b'5wIlN2p+SZM9hMqI4INJGF0ctbTkyx5CpUZw6JxP/DE0DZfE2UFwjSHBuYRJXmfT3JpBeyCNKcHx'
b'+mPMIEVzaybtgTTG5Plt5FygptU1zBi1B9I4E8C74B1vcdNozak9kMaaYCKiS5ldikZrVu2BNN4E'
b'w7kJgNQ0Mmte7YE0akRnffIemmLXBNmDaASJxkbej9CUuKbJHkQxoiFg5JfKfF2VVJ8XqUJQY1rP'
b'6QkxNaUx8IRvOxDeNE6NkYjQNcWvGbcH0Zg1GDToefKq0JpxexC2Mi5TjTNIYdTiknl7CK5wbEgW'
b'kucEqEpc8m/Pv2Qfv4NAiGWgCfPmLpvJ+/T0cHl83R93r/vL4bH5Px0f/uGW/M6HwrLAyVCkm3Cv'
b'jsFriLUZBTQpS5w0o5veSetIb1w/tubUurz+8uU3Tr4Ng7Q+9TaMYFrLyuUrSz9h3q2dCK1LBeNB'
b'QLxtfFufbwHYjd71xijc2MWGFLKMM87dXId/bwBtPjKU4jY2gJsAYezjYibNCHOC2gSIlUMScQ0W'
b'GRurFHnjTqR5Z4GTHuO9syhYbgsLxx4DrA3oRs1PeowbsWjabeDYbeCs29Cx+ERnJh2Hk1xtWlhl'
b'voIyL/chUtOSuZJ8CtCs5HtiamtauJJ+CpCbFluDQE2z82K7CeSnpc5TbMq81G3ChLHccaFNzI5V'
b'oSspufXqY5wWPnnhVQmzwre5pjQpfkkqV1XirPhtwuB7AbRceWDU7KQA3hhZ2U+7+/V3TNxzhh+w'
b'279cBk+B6Vykn6WBXCnxFyzXSwWK8+3+LrX2kD9SbmdFgu7G87JbPSHvPv1mkzQBdXLk75NPKHkI'
b'03ZvZSiwbroaK7D2vXnkz1P4uGl2inbFHcQ9Z/IfxzJZ8rIlpwPvsWdgkW7W+a13xAuIrz/3TtRr'
b'NyduYXDfLw8vO7gPxnzawY5fDe8fJ+qK+8IghhmgvWNrOKp8/ROv4/B2lmNXRaqruA6UkiA1NoYo'
b'xwGIYBT/EkkUP25olUhmzKtGQtdxUkBMWeIET40HJVXRVEuucS5Lr/iXeCRXkZMolYXsrgU0vtjt'
b'kNB2oPjAVeXVFOtEJXmDzQE8OKEQGXSQHBrK94JSH5V3teBcNbTQgUp8dBOBM0A9MEmjekai1XfU'
b'gyhhBXHl3fhrvRT5lj3JH7NvBl027k4LELuglijKJybD9B6mskwnyJxb9fBgPRMwJLKjpuGUSFGu'
b'Xs1gb410drzQIuWy0gWNSPw684ATUDVSklSjaim3TYGDjKOWFJwcKQk384bQrZHODgD9SFslJN+H'
b'dk7OuHkIE2g9XslBqubgxJbknFc1p+CUeKWMkLTzNyfx/Miivtyk4MqlEZTBznD1eHMmYluCczE3'
b'FKOGKlQJWQoOhwT+5pA7jcH2DoDcufbmcWQSl8o8xBsakOvTcP+bJ+p0OWDZmaMuI/2fafJuSQQC'
b'L9XkShw/Ojs2r1P6FjkNmk/bIAgyhHWp6e+Kz9itQUKud8ZJj1iUK1VrHbnLEFxs8OZMmZ9i9bij'
b'BsxHrFCG+fF4I+qUkSdXgXwa5aVVBVSoGnXMEDiEW4Oe3QbpIZOCmuQ95yFMUPV4IScn1GUAgk/S'
b'z48aKUg1XMoIMMRbw51dX6nhcrr1UflXWwY3Qd0IN2ckVAfI5wbiHrhppCCVcMFmBOmJb87p+dWO'
b'HnFQgJEMl+I8pBnwRtA5IUdjMtai874pQYWqUUuV/A80XjlL', b''])


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

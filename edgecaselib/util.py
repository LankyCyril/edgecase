from sys import stderr
from contextlib import contextmanager, ExitStack
from itertools import chain
from regex import compile
from re import split


MAINCHROMS_ENSEMBL = {str(i) for i in range(1, 23)} | {"X", "Y"}
MAINCHROMS_UCSC = {"chr" + s for s in MAINCHROMS_ENSEMBL}
MAINCHROMS_T2T = {"chrX_fixedBionanoSV_centromereV3"}

ALPHABET = list("ACGT")
COMPLEMENTS = dict(zip(ALPHABET, reversed(ALPHABET)))
MOTIF_COMPLEMENTS = {**COMPLEMENTS, **{"[": "]", "]": "[", ".": "."}}
MOTIF_COMPLEMENT_PATTERN = compile(r'|'.join(MOTIF_COMPLEMENTS.keys()))


@contextmanager
def ReadFileChain(filenames, manager):
    """Chain records from all filenames in list bams, replicating behavior of pysam context managers"""
    with ExitStack() as stack:
        yield chain(*(
            stack.enter_context(manager(filename))
            for filename in filenames
        ))


def motif_revcomp(motif, ignorecase=True):
    """Reverse-complement a regex-like motif; only allows a subset of regex syntax (ACGT, dot, square brackets)"""
    if ignorecase:
        matcher = lambda match: MOTIF_COMPLEMENTS[match.group().upper()]
    else:
        matcher = lambda match: MOTIF_COMPLEMENTS[match.group()]
    try:
        return MOTIF_COMPLEMENT_PATTERN.sub(matcher, motif[::-1])
    except KeyError:
        raise ValueError("Unsupported character(s) in motif: {}".format(motif))


def chromosome_natsort(chrom):
    """Natural order sorting that undestands chr1, 4, chr10, chr14_K*, 7ptel etc"""
    keyoder = []
    for chunk in split(r'(\d+)', chrom): # stackoverflow.com/a/16090640
        if chunk.isdigit():
            keyoder.append(int(chunk))
        elif chunk == "":
            keyoder.append("chr")
        else:
            keyoder.append(chunk.lower())
    return keyoder


def natsorted_chromosomes(chromosomes):
    """Sort chromosome names in natural order"""
    try:
        return sorted(chromosomes, key=chromosome_natsort)
    except Exception as e:
        msg = "natural sorting failed, pages will be sorted alphanumerically"
        print("Warning: " + msg, file=stderr)
        print("The error was: '{}'".format(e), file=stderr)
        return sorted(chromosomes)

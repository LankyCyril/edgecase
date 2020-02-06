from sys import stderr
from regex import compile
from re import split, search, IGNORECASE
from shutil import which
from os import path, access, X_OK
from functools import partial
from tqdm import tqdm


MAINCHROMS_ENSEMBL = {str(i) for i in range(1, 23)} | {"X", "Y"}
MAINCHROMS_UCSC = {"chr" + s for s in MAINCHROMS_ENSEMBL}
MAINCHROMS_T2T = {"chrX_fixedBionanoSV_centromereV3"}

ALPHABET = list("ACGT")
COMPLEMENTS = dict(zip(ALPHABET, reversed(ALPHABET)))
MOTIF_COMPLEMENTS = {**COMPLEMENTS, **{"[": "]", "]": "[", ".": "."}}
MOTIF_COMPLEMENT_PATTERN = compile(r'|'.join(MOTIF_COMPLEMENTS.keys()))


progressbar = partial(
    tqdm, bar_format=(
        "{desc}{percentage:3.0f}% ({n_fmt}/{total_fmt}), " +
        "{elapsed}<{remaining}, {rate_fmt}"
    )
)


def validate_motif(motif):
    """Make sure motif conforms to what we've implemented"""
    return (search(r'^[ACGT\[\]\.]*$', motif, flags=IGNORECASE) is not None)


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


def get_executable(exe_name, suggested_binary, fail_if_none=True):
    """Wrapper to find executable"""
    if suggested_binary is None:
        binary = which(exe_name)
        if (binary is None) and fail_if_none:
            raise OSError("{} not found".format(exe_name))
        else:
            return binary
    else:
        if path.isfile(suggested_binary) and access(suggested_binary, X_OK):
            return suggested_binary
        elif fail_if_none:
            raise OSError("{} not accessible".format(suggested_binary))
        else:
            return None

from contextlib import contextmanager, ExitStack
from itertools import chain
from regex import compile


MAINCHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}
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

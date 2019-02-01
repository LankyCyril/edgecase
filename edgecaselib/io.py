from contextlib import contextmanager, ExitStack
from itertools import chain


MAINCHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}
ALPHABET = list("ACGT")


@contextmanager
def ReadFileChain(filenames, manager):
    """Chain records from all filenames in list bams, replicating behavior of pysam context managers"""
    with ExitStack() as stack:
        yield chain(*(
            stack.enter_context(manager(filename))
            for filename in filenames
        ))

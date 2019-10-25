from sys import stdout
from pysam import AlignmentFile
from shutil import which
from os import path, access, X_OK
from re import search
from subprocess import Popen, PIPE, check_output
from tqdm import tqdm


from contextlib import contextmanager
@contextmanager
def TemporaryDirectory():
    """Temporary TemporaryDirectory plug for development purposes"""
    yield "data/datasets/twins/sandbox/edge-meme"


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


def guess_bg_fmt(background):
    """Decide if `background` is in SAM/BAM format or is a MEME HMM"""
    with open(background, mode="rb") as bg_handle:
        line = next(bg_handle)
    if search(br'^#.*Markov frequencies', line):
        return "hmm"
    else:
        return "sam"


def interpret_args(fmt, fasta_get_markov, bioawk, samtools, meme, background):
    """Parse and check arguments"""
    if fmt == "sam":
        manager = AlignmentFile
    elif fmt == "fastx":
        raise NotImplementedError("--fmt fastx")
    else:
        raise ValueError("Unknown --fmt: '{}'".format(fmt))
    bg_fmt = guess_bg_fmt(background)
    if bg_fmt == "sam":
        if get_executable("fasta-get-markov", fasta_get_markov, False) is None:
            raise ValueError("need fasta-get-markov for background, not found")
        if get_executable("bioawk", bioawk, False) is None:
            raise ValueError("need bioawk for background, not found")
        if get_executable("samtools", samtools, False) is None:
            raise ValueError("need samtools for background, not found")
    return (
        manager,
        get_executable("fasta-get-markov", fasta_get_markov, False),
        get_executable("bioawk", bioawk, False),
        get_executable("samtools", samtools, False),
        get_executable("meme", meme),
        bg_fmt
    )


def convert_background(sam, tempdir, fasta_get_markov, bioawk, samtools, max_order=6):
    """Convert a SAM/BAM file into Markov background for MEME"""
    samtools_view = Popen(
        [samtools, "view", "-F3844", sam], stdout=PIPE
    )
    bioawk_conv = Popen(
        [bioawk, "-c", "sam", '{print ">"$qname; print $seq}'],
        stdin=samtools_view.stdout, stdout=PIPE
    )
    hmm = check_output(
        [fasta_get_markov, "-m", str(max_order)],
        stdin=bioawk_conv.stdout
    )
    bfile = path.join(tempdir, "bfile")
    with open(bfile, mode="wt") as bfile_handle:
        print(hmm.decode(), file=bfile_handle)
    return bfile


def main(readfiles, fmt, fasta_get_markov, bioawk, samtools, meme, minw, background, maxw, jobs=1, file=stdout.buffer, **kwargs):
    # parse arguments
    manager, fasta_get_markov, bioawk, samtools, meme, bg_fmt = interpret_args(
        fmt, fasta_get_markov, bioawk, samtools, meme, background
    )
    with TemporaryDirectory() as tempdir:
        if bg_fmt == "sam": # will need to convert SAM to HMM
            background = convert_background(
                background, tempdir, fasta_get_markov, bioawk, samtools
            )

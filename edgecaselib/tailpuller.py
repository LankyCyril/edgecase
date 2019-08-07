from sys import stdout
from os.path import isfile
from edgecaselib.util import ReadFileChain
from pysam import AlignmentFile
from pandas import read_csv, read_fwf
from functools import reduce
from operator import __or__


def load_index(basename):
    """Load FAI and ECX indices"""
    fai_filename, ecx_filename = basename + ".fai", basename + ".ecx"
    if not isfile(fai_filename):
        raise FileNotFoundError(fai_filename)
    elif not isfile(ecx_filename):
        raise FileNotFoundError(ecx_filename)
    else:
        ecx = read_fwf(ecx_filename, skiprows=1)
        ecx.columns = [c.lstrip("#") for c in ecx.columns]
        fai = read_csv(
            fai_filename, sep="\t",
            usecols=(0,1), index_col=0, names=["rname", "length"]
        )
        return ecx, dict(fai["length"])


def filter_entries(bam_data, ecx, flag_filter):
    """Only pass reads extending past regions specified in the ECX"""
    ecx_rnames = set(ecx["rname"])
    for entry in bam_data:
        passes_flags = (entry.flag & flag_filter == 0)
        if passes_flags and (entry.reference_name in ecx_rnames):
            # measure the None (clipped) stretch on the left:
            p_cigartype, p_clip = entry.cigartuples[0]
            if (p_cigartype != 4) and (p_cigartype != 5):
                p_clip = 0
            # determine location of the start of the read relative to reference:
            p_anchor_pos = entry.reference_start - p_clip
            # find flags in ECX where anchor is to the right of the read start:
            p_flags = set(ecx.loc[
                (ecx["rname"]==entry.reference_name) & (ecx["prime"]==5) &
                (ecx["pos"]>=p_anchor_pos), "flag"
            ])
            # measure the None (clipped) stretch on the right:
            q_cigartype, q_clip = entry.cigartuples[-1]
            if (q_cigartype != 4) and (q_cigartype != 5):
                q_clip = 0
            # determine location of the end of the read relative to reference:
            q_anchor_pos = entry.reference_end + q_clip
            # find flags in ECX where anchor is to the left of the read end:
            q_flags = set(ecx.loc[
                (ecx["rname"]==entry.reference_name) & (ecx["prime"]==3) &
                (ecx["pos"]<=q_anchor_pos), "flag"
            ])
            extra_flags = p_flags | q_flags
            if extra_flags: # update entry flags with edgeCase flags and yield
                entry.flag |= reduce(__or__, extra_flags)
                if q_flags: # add the is_q flag
                    entry.flag |= 0x8000
                yield entry


def main(bams, reference, flag_filter, file=stdout, **kwargs):
    # use header of first input file (NB! fragile):
    with AlignmentFile(bams[0]) as bam:
        print(str(bam.header).rstrip("\n"), file=file)
    # dispatch data to subroutines:
    ecx, reflens = load_index(reference)
    with ReadFileChain(bams, AlignmentFile) as bam_data:
        for entry in filter_entries(bam_data, ecx, flag_filter):
            print(entry.to_string(), file=file)

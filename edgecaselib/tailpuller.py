from sys import stdout
from os.path import isfile
from edgecaselib.util import ReadFileChain
from pysam import AlignmentFile
from pandas import read_fwf
from functools import reduce
from operator import __or__
from copy import deepcopy


def load_index(index_filename):
    """Load FAI and ECX indices"""
    if not isfile(index_filename):
        raise FileNotFoundError(index_filename)
    else:
        ecx = read_fwf(index_filename, skiprows=1)
        ecx.columns = [c.lstrip("#") for c in ecx.columns]
        return ecx


def make_ecx_filter_dict(ecx):
    """Make simple dictionary for filtering purposes from the ECX"""
    slim_ecx = ecx[["rname", "pos", "flag", "prime", "rc"]]
    rnames = set(slim_ecx["rname"])
    return {
        rname: slim_ecx[
                slim_ecx["rname"]==rname
            ].drop(columns="rname").set_index(["prime", "rc"])
        for rname in rnames
    }


def get_anchor_pos(reference_pos, cigartuples, cigarpos):
    """Calculate the position of clipped start/end of read relative to the reference"""
    if not cigartuples: # no CIGAR available
        return None
    # measure the None (clipped) stretch on the left (0) or right (-1):
    cigartype, clip = cigartuples[cigarpos]
    if (cigartype != 4) and (cigartype != 5): # not a soft/hard clip
        clip = 0
    # determine location of start/end of the read relative to reference:
    if cigarpos == 0: # start
        return reference_pos - clip
    elif cigarpos == -1: # end
        return reference_pos + clip
    else:
        raise ValueError("get_anchor_pos(): cigarpos can only be 0 or -1")


def updated_entry(entry, flags, is_q=False):
    """Add ECX flags to entry"""
    new_entry = deepcopy(entry)
    new_entry.flag |= reduce(__or__, flags)
    if is_q:
        new_entry.flag |= 0x8000
    return new_entry


def filter_entries(bam_data, ecxfd, flag_filter):
    """Only pass reads extending past regions specified in the ECX"""
    for entry in bam_data:
        if (entry.flag & flag_filter == 0) and (entry.reference_name in ecxfd):
            ecx_entry = ecxfd[entry.reference_name]
            p_anchor_pos = get_anchor_pos( # pos of start relative to reference
                entry.reference_start, entry.cigartuples, cigarpos=0
            )
            p_flags = set( # ECX flags where anchor is to right of read start
                ecx_entry[ecx_entry["pos"]>=p_anchor_pos].loc[[5, "-"], "flag"]
            )
            if p_flags:
                yield updated_entry(entry, p_flags)
            q_anchor_pos = get_anchor_pos( # pos of end relative to reference
                entry.reference_end, entry.cigartuples, cigarpos=-1
            )
            q_flags = set( # ECX flags where anchor is to left of read end
                ecx_entry[ecx_entry["pos"]<q_anchor_pos].loc[[3, "-"], "flag"]
            )
            if q_flags:
                yield updated_entry(entry, q_flags, is_q=True)


def main(bams, index, flag_filter, file=stdout, **kwargs):
    # use header of first input file (NB! fragile):
    with AlignmentFile(bams[0]) as bam:
        print(str(bam.header).rstrip("\n"), file=file)
    # dispatch data to subroutines:
    ecx = load_index(index)
    with ReadFileChain(bams, AlignmentFile) as bam_data:
        ecxfd = make_ecx_filter_dict(ecx)
        for entry in filter_entries(bam_data, ecxfd, flag_filter):
            print(entry.to_string(), file=file)

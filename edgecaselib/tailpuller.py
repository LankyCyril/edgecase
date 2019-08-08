from sys import stdout
from os.path import isfile
from pysam import AlignmentFile
from pandas import read_csv
from functools import reduce
from operator import __or__
from copy import deepcopy
from tqdm import tqdm


def load_index(index_filename, as_filter_dict=False):
    """Load ECX index; convert to simpler dictionary for filtering if requested"""
    if not isfile(index_filename):
        raise FileNotFoundError(index_filename)
    else:
        ecx = read_csv(index_filename, sep="\t", skiprows=1, escapechar="#")
        if as_filter_dict:
            slim_ecx = ecx[["rname", "pos", "flag", "prime"]]
            rnames = set(slim_ecx["rname"])
            return {
                rname: {
                    5: slim_ecx[
                        (slim_ecx["rname"]==rname) & (slim_ecx["prime"]==5)
                    ].drop(columns=["rname", "prime"]),
                    3: slim_ecx[
                        (slim_ecx["rname"]==rname) & (slim_ecx["prime"]==3)
                    ].drop(columns=["rname", "prime"])
                }
                for rname in rnames
            }
        else:
            return ecx


def get_terminal_pos(entry, cigarpos):
    """Calculate the position of clipped start/end of read relative to the reference"""
    if not entry.cigartuples: # no CIGAR available
        return None
    # measure the clipped stretch on the left (0) or right (-1):
    cigartype, clip = entry.cigartuples[cigarpos]
    if (cigartype != 4) and (cigartype != 5): # not a soft/hard clip
        clip = 0
    # determine location of start/end of the read relative to reference:
    if cigarpos == 0: # start
        return entry.reference_start - clip
    elif cigarpos == -1: # end
        return entry.reference_end + clip
    else:
        raise ValueError("get_terminal_pos(): cigarpos can only be 0 or -1")


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
            # find positions of start and end of read relative to reference:
            p_pos = get_terminal_pos(entry, cigarpos=0)
            q_pos = get_terminal_pos(entry, cigarpos=-1)
            # collect ECX flags where anchor is to right of read start:
            ecx_t_p5 = ecxfd[entry.reference_name][5]
            p_flags = set(ecx_t_p5.loc[ecx_t_p5["pos"]>=p_pos, "flag"])
            if p_flags:
                yield updated_entry(entry, p_flags)
            # collect ECX flags where anchor is to left of read end
            ecx_t_p3 = ecxfd[entry.reference_name][3]
            q_flags = set(ecx_t_p3.loc[ecx_t_p3["pos"]<q_pos, "flag"])
            if q_flags:
                yield updated_entry(entry, q_flags, is_q=True)


def main(bam, index, flag_filter, file=stdout, **kwargs):
    # dispatch data to subroutines:
    ecxfd = load_index(index, as_filter_dict=True)
    with AlignmentFile(bam) as bam_data:
        print(str(bam_data.header).rstrip("\n"), file=file)
        for reference in tqdm(bam_data.references, desc="reference"):
            bam_chunk = bam_data.fetch(reference, None, None)
            for entry in filter_entries(bam_chunk, ecxfd, flag_filter):
                print(entry.to_string(), file=file)

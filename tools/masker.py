#!/usr/bin/env python
from sys import argv
from pysam import AlignmentFile
from regex import finditer, IGNORECASE
from itertools import count
from os import path


def get_circular_motif(motif):
    return r'|'.join({motif[i:]+motif[:i] for i in range(len(motif))})


def get_motifs_pattern(motif_description):
    if path.isfile(motif_description):
        with open(motif_description, mode="rt") as motif_handle:
            return r'|'.join({
                get_circular_motif(line.split()[0]*2)
                for line in motif_handle if line[0] != "#"
            })
    else:
        return r'|'.join({
            get_circular_motif(kmer*2)
            for kmer in motif_description.split("|")
        })

def iterate_bam(filename, f, F):
    with AlignmentFile(filename) as bam:
        for entry in bam:
            if (entry.flag & f == f) and (entry.flag & F == 0) and (entry.seq):
                yield entry


def cut_chunks(seq, positions_to_mask, minlen):
    prev_pos, start, end = -1, 0, 0
    for pos in range(len(seq)):
        if pos not in positions_to_mask:
            if pos == prev_pos + 1:
                end = pos
            else:
                start = pos
            prev_pos = pos
        else:
            if start < end:
                if end - start + 1 >= minlen:
                    yield seq[start:end+1]
            start = end + 1


def main(filename, motif_description, action="mask", f=49152, F=3840, minlen=22):
    minlen = int(minlen)
    motifs_pattern = get_motifs_pattern(motif_description)
    cid = count(1)
    total_masked, total_bases = 0, 0
    for entry in iterate_bam(filename, int(f), int(F)):
        positions_to_mask = set()
        matcher = finditer(
            motifs_pattern, entry.seq, overlapped=True, flags=IGNORECASE
        )
        for match in matcher:
            positions_to_mask |= set(range(match.start(), match.end()))
        if action == "cut":
            for chunk in cut_chunks(entry.seq, positions_to_mask, minlen):
                print(">{}\n{}".format(next(cid), chunk))
        elif action == "measure":
            positions_to_keep = set(range(len(entry.seq)+1)) - positions_to_mask
            for chunk in cut_chunks(entry.seq, positions_to_keep, minlen):
                print(len(chunk))
        elif action == "mask":
            raise NotImplementedError("`action` 'mask'")
        elif action == "count":
            print("{}\t{}".format(
                len(positions_to_mask), len(entry.seq) - len(positions_to_mask)
            ))
        elif action == "total-fraction":
            total_masked += len(positions_to_mask)
            total_bases += len(entry.seq)
        else:
            raise ValueError("Unknown action: '{}'".format(action))
    if action == "total-fraction":
        print(total_masked, total_bases, total_masked/total_bases, sep="\t")
    return 0


if __name__ == "__main__":
    returncode = main(*argv[1:])
    exit(returncode)

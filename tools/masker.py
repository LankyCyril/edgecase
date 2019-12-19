#!/usr/bin/env python
from sys import argv
from pysam import AlignmentFile
from regex import finditer, IGNORECASE
from itertools import count


def get_circular_motif(motif):
    return r'|'.join({motif[i:]+motif[:i] for i in range(len(motif))})


def iterate_bam(filename, f, F):
    with AlignmentFile(filename) as bam:
        for entry in bam:
            if (entry.flag & f == f) and (entry.flag & F == 0) and (entry.seq):
                yield entry


def main(filename, motif_filename, f=49152, F=3840, minlen=22):
    minlen = int(minlen)
    with open(motif_filename, mode="rt") as motif_handle:
        motifs_pattern = r'|'.join({
            get_circular_motif(line.split()[0]*2)
            for line in motif_handle if line[0] != "#"
        })
    cid = count(1)
    for entry in iterate_bam(filename, int(f), int(F)):
        positions_to_mask = set()
        matcher = finditer(
            motifs_pattern, entry.seq, overlapped=True, flags=IGNORECASE
        )
        for match in matcher:
            positions_to_mask |= set(range(match.start(), match.end()))
        prev_pos, start, end = -1, 0, 0
        for pos in range(len(entry.seq)):
            if pos not in positions_to_mask:
                if pos == prev_pos + 1:
                    end = pos
                else:
                    start = pos
                prev_pos = pos
            else:
                if start < end:
                    if end - start + 1 >= minlen:
                        print(">{}".format(next(cid)))
                        print(entry.seq[start:end+1])
                start = end + 1
    return 0


if __name__ == "__main__":
    returncode = main(*argv[1:])
    exit(returncode)

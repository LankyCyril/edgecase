#!/usr/bin/env python
from sys import argv, stdout
from gzip import open as gzopen
from functools import lru_cache
from collections import defaultdict
from tqdm import tqdm
from pandas import DataFrame


""" Generate the input gzipped file with:
./kmer-counter --double-count-palindromes --k 12 --fasta ${INPUT_FILE} \
| cut -f2 | sed 's/$/ ;/g' | tr ' ' '\n' \
| awk -F':' '{if (($1==";") || (substr($1,1,6)==substr($1,7))) {print}}' \
| gzip -3 > ${OUTPUT_FILE}
"""


@lru_cache(maxsize=None)
def lowest_alpha_inversion(kmer):
    return min(kmer[i:]+kmer[:i] for i in range(len(kmer)))


@lru_cache(maxsize=None)
def get_motif_identity(kmer, min_repeats=2):
    lai = lowest_alpha_inversion(kmer)
    motif = lai[:int(len(lai)/min_repeats)]
    if motif * min_repeats == lai:
        return motif


def main(raw_kmercount_filename):
    motif_count_database, counts = [], defaultdict(int)
    with gzopen(raw_kmercount_filename, mode="rt") as raw_kc:
        for line in tqdm(map(str.strip, raw_kc), total=3506715):
            if (line == ";") and counts:
                motif_count_database.append(counts)
                counts = defaultdict(int)
            else:
                stat = line.split(":")
                if len(stat) == 2:
                    counts[get_motif_identity(stat[0])] += int(stat[1])
    DataFrame(motif_count_database).to_csv(stdout, sep="\t", index=False)


if __name__ == "__main__":
    exit(main(argv[1]) or 0)

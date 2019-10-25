#!/usr/bin/env python
from sys import argv
from os import path
from multiprocessing import cpu_count
from subprocess import check_output

def main(*, fasta, bfile, oc, nmotifs, numiter):
    with open(fasta, mode="rt") as fh:
        fasta_data = "".join(fh.readlines()).strip()
    tempfasta = path.join(oc, "original.fa")
    with open(tempfasta, mode="wt") as fh:
        print(fasta_data, file=fh)
    for i in range(numiter):
        meme_output = check_output([
            "meme", "-dna", "-mod", "anr", "-p", str(int(cpu_count() / 2)),
            "-minw", "4", "-maxw", "50", "-minsites", "2",
            "-nmotifs", str(nmotifs),
            "-bfile", bfile, "-text", tempfasta
        ])
        meme_output_file = path.join(oc, "{:05d}.meme".format(i))
        with open(meme_output_file, mode="wt") as mh:
            print(meme_output.decode().strip(), file=mh)
        masked_fasta = check_output([
            "tools/meme-masker.py", tempfasta, meme_output_file,
            "--aggressive"
        ])
        tempfasta = path.join(oc, "{:05d}.fa".format(i))
        with open(tempfasta, mode="wt") as fh:
            print(masked_fasta.decode(), file=fh)
    return 0

if __name__ == "__main__":
    returncode = main(
        fasta=argv[1], bfile=argv[2], oc=argv[3],
        nmotifs=int(argv[4]), numiter=int(argv[5])
    )
    exit(returncode)

from os.path import join
from argparse import Namespace
from edgecaselib import tailpuller, tailchopper, kmerscanner
from gzip import open as gzopen

configfile: "config.yaml"

rule tailpuller:
    input:
        ar=join(config["data_dir"], config["reads_dir"], "{dataset}/AR.sorted.bam"),
        reference=config["reference"]
    output:
        ac=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam")
    run:
        with open(output.ac, mode="wt") as sam:
            tailpuller.main(
                args=Namespace(
                    bams=input.ar,
                    reference=input.reference,
                    prime=int(wildcards.prime)
                ),
                file=sam
            )

rule tailchopper:
    input: sam=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam")
    output: fa=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}OOB.fa")
    run:
        with open(output.fa, mode="wt") as fasta:
            tailchopper.main(
                args=Namespace(bams=input, prime=int(wildcards.prime)),
                file=fasta
            )

rule candidate_densities:
    input: sam=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam")
    output: dat=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC-densities.dat")
    params: kmer="TTAGGG", window_size=120
    threads: 1
    run:
        with gzopen(output.dat, mode="wt") as dat:
            kmerscanner.main(
                args=Namespace(
                    bams=input, num_reads=None,
                    kmer=params.kmer, window_size=params.window_size,
                    head_test=None, tail_test=None, cutoff=None,
                    jobs=threads
                ),
                file=dat
            )

rule all_dataset_tails:
    input:
        targets=expand(
            join(config["data_dir"], config["reads_dir"], "{timepoint}_{subject}/{prime}{kind}"),
            timepoint=config["timepoints"],
            subject=config["subjects"],
            prime=[5, 3],
            kind=["AC.sam", "OOB.fa", "AC-density.dat"]
        )

rule all_dataset_candidate_densities:
    input:
        targets=expand(
            join(config["data_dir"], config["reads_dir"], "{timepoint}_{subject}/{prime}AC-densities.dat"),
            timepoint=config["timepoints"],
            subject=config["subjects"],
            prime=[5, 3]
        )

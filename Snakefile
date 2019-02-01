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
                    bams=[input.ar],
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
                args=Namespace(bams=[input.sam], prime=int(wildcards.prime)),
                file=fasta
            )

rule candidate_densities:
    input:
        sam=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam"),
        motifs=join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-motifs.txt")
    output:
        dat=join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-densities.dat")
    params: window_size=120
    threads: 16
    run:
        with open(input.motifs) as motif_data:
            motifs = {line.strip() for line in motif_data}
        with gzopen(output.dat, mode="wt") as dat:
            for motif in motifs:
                kmerscanner.main(
                    args=Namespace(
                        bams=[input.sam], num_reads=None,
                        motif=motif, window_size=params.window_size,
                        head_test=None, tail_test=None, cutoff=None,
                        jobs=threads
                    ),
                    file=dat
                )

rule all_dataset_tails:
    input:
        targets=expand(
            join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}{kind}"),
            dataset=config["datasets"],
            prime=[5, 3],
            kind=["AC.sam", "OOB.fa"]
        )

rule all_dataset_candidate_densities:
    input:
        targets=expand(
            join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-densities.dat"),
            dataset=config["datasets"],
            prime=[5, 3]
        )

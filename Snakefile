from os.path import join
from argparse import Namespace
from edgecaselib import tailpuller, tailchopper

configfile: "config.yaml"

rule tailpuller:
    input:
        ar=join(config["data_dir"], config["reads_dir"], "{dataset}/AR.sorted.bam"),
        reference=config["reference"]
    output:
        ac=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam")
    run:
        with open(output.ac, "wt") as sam:
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
        with open(output.fa, "wt") as fasta:
            tailchopper.main(
                args=Namespace(bams=input, prime=int(wildcards.prime)),
                file=fasta
            )

rule all_twins_tails:
    input:
        targets=expand(
            join(config["data_dir"], config["reads_dir"], "{timepoint}_{subject}/{prime}{kind}"),
            timepoint=config["timepoints"],
            subject=config["subjects"],
            prime=[5, 3],
            kind=["AC.sam", "OOB.fa"]
        )

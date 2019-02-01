from os.path import join
from argparse import Namespace
from edgecaselib import tailpuller, tailchopper

DATA_DIR = config.get("DATA_DIR", "data")
READS_DIR = config.get("READS_DIR", "processedReads")
REFERENCE = config.get("REFERENCE", "data/references/hs37d5.fa")

rule tailpuller:
    input:
        ar=join(DATA_DIR, READS_DIR, "{dataset}/AR.sorted.bam"),
        reference=REFERENCE
    output:
        ac=join(DATA_DIR, READS_DIR, "{dataset}/{prime}AC.sam")
    run:
        with open(output.ac, "wt") as sam:
            tailpuller.main(
                args=Namespace(
                    bams=input,
                    reference=REFERENCE,
                    prime=int(wildcards.prime)
                ),
                file=sam
            )

rule tailchopper:
    input: sam=join(DATA_DIR, READS_DIR, "{dataset}/{prime}AC.sam")
    output: fa=join(DATA_DIR, READS_DIR, "{dataset}/{prime}OOB.fa")
    run:
        with open(output.fa, "wt") as fasta:
            tailchopper.main(
                args=Namespace(bams=input, prime=int(wildcards.prime)),
                file=fasta
            )

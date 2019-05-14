from os.path import join
from edgecaselib import tailpuller, tailchopper, kmerscanner, densityplot
from edgecaselib.util import motif_revcomp
from gzip import open as gzopen
from pysam import AlignmentFile

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
                bams=[input.ar],
                reference=input.reference,
                prime=int(wildcards.prime),
                names=config["chromosome_names"],
                file=sam
            )

rule tailchopper:
    input: sam=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam")
    output: fa=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}OOB.fa")
    run:
        with open(output.fa, mode="wt") as fasta:
            tailchopper.main(
                bams=[input.sam],
                prime=int(wildcards.prime),
                file=fasta
            )

rule ac_fasta:
    input: join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam")
    output: join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.fa")
    shell: """
        samtools view {input} \
        | bioawk -c sam '{{print ">"$qname; print $seq}}' \
        > {output}
    """

rule ar_fasta:
    input: join(config["data_dir"], config["reads_dir"], "{dataset}/AR.sorted.bam")
    output: join(config["data_dir"], config["reads_dir"], "{dataset}/AR.fa.gz")
    shell: """
        samtools view {input} \
        | bioawk -c sam '{{if ($seq != "*") {{print ">"$qname; print $seq}}}}' \
        | gzip -2 > {output}
    """

rule ar_ib_fasta:
    input:
        ar=join(config["data_dir"], config["reads_dir"], "{dataset}/AR.fa.gz"),
        p5ac=join(config["data_dir"], config["reads_dir"], "{dataset}/5AC.fa"),
        p3ac=join(config["data_dir"], config["reads_dir"], "{dataset}/3AC.fa")
    output:
        fa=join(config["data_dir"], config["reads_dir"], "{dataset}/AR-IB.fa.gz")
    shell: """
        zcat {input.ar} \
        | paste - - \
        | grep -ivFf <(cat {input.p5ac} {input.p3ac} | grep -iF '>') \
        | tr '\\t' '\\n' \
        | gzip -2 > {output.fa}
    """

rule candidate_densities:
    input:
        sam=join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}AC.sam"),
        motifs=join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-motifs.txt")
    output:
        dat=join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-densities.dat")
    params: window_size=120, revcomp=True
    threads: 16
    run:
        with open(input.motifs) as motif_data:
            motifs = {line.strip() for line in motif_data}
        if params.revcomp:
            motifs |= {motif_revcomp(motif) for motif in set(motifs)}
        with gzopen(output.dat, mode="wt") as dat:
            for i, motif in enumerate(motifs):
                if i == 0:
                    no_print_header = False
                else:
                    no_print_header = True
                kmerscanner.main(
                    bams=[input.sam], num_reads=None,
                    motif=motif, window_size=params.window_size,
                    head_test=None, tail_test=None, cutoff=None,
                    jobs=threads,
                    no_print_header=no_print_header, file=dat
                )

rule densityplot:
    input:
        dat=join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-densities.dat")
    output:
        pdf=join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-densities.pdf")
    params: bin_size=100
    run:
        with open(output.pdf, mode="wb") as pdf:
            densityplot.main(
                dat=input.dat,
                bin_size=params.bin_size,
                title="{} {}AC".format(wildcards.dataset, wildcards.prime),
                file=pdf
            )

rule all_dataset_tails:
    input:
        targets=expand(
            join(config["data_dir"], config["reads_dir"], "{dataset}/{prime}{kind}"),
            dataset=config["datasets"],
            prime=[5, 3],
            kind=["AC.sam", "OOB.fa"]
        )

rule all_dataset_backgrounds:
    input:
        targets=expand(
            join(config["data_dir"], config["reads_dir"], "{dataset}/{kind}.fa.gz"),
            dataset=config["datasets"],
            kind=["AR", "AR-IB"]
        )

rule all_dataset_candidate_densities:
    input:
        targets=expand(
            join(config["data_dir"], config["analysis_dir"], "{dataset}/{prime}AC-densities.pdf"),
            dataset=config["datasets"],
            prime=[5, 3]
        )

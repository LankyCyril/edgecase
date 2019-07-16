from edgecaselib import tailpuller, tailchopper, kmerscanner, densityplot
from edgecaselib.util import motif_revcomp
from gzip import open as gzopen

configfile: "config.yaml"

rule tailpuller:
    input: ar="{path}/AR.sorted.bam", reference=config["reference"]
    output: ac="{path}/{prime}AC.sam"
    run:
        with open(output.ac, mode="wt") as sam:
            tailpuller.main(
                bams=[input.ar], reference=input.reference,
                prime=int(wildcards.prime), names=config["chromosome_names"],
                file=sam
            )

rule tailchopper:
    input: sam="{path}/{prime}AC.sam"
    output: fa="{path}/{prime}OOB.fa"
    run:
        with open(output.fa, mode="wt") as fasta:
            tailchopper.main(
                bams=[input.sam], prime=int(wildcards.prime), file=fasta
            )

rule ac_fasta:
    input: "{path}/{prime}AC.sam"
    output: "{path}/{prime}AC.fa"
    shell: """
        samtools view {input} \
        | bioawk -c sam '{{print ">"$qname; print $seq}}' \
        > {output}
    """

rule ar_fasta:
    input: bam="{path}/AR.sorted.bam"
    output: fa="{path}/AR.fa.gz", fai="{path}/AR.fa.gz.fai", gzi="{path}/AR.fa.gz.gzi"
    shell: """
        samtools view -F2304 {input} \
        | bioawk -c sam '{{if ($seq != "*") {{print ">"$qname; print $seq}}}}' \
        | bgzip > {output};
        samtools faidx {output}
    """

rule ar_ib_fasta:
    input: ar="{path}/AR.fa.gz", p5ac="{path}/5AC.fa", p3ac="{path}/3AC.fa"
    output: fa="{path}/AR-IB.fa.gz"
    shell: """
        zcat {input.ar} \
        | paste - - \
        | grep -ivFf <(cat {input.p5ac} {input.p3ac} | grep -iF '>') \
        | tr '\\t' '\\n' \
        | gzip -2 > {output.fa}
    """

rule hmmer:
    input: ar="{path}/AR.fa.gz", model=config["telomere_model"]
    output: tbl="{path}/AR.hmm.tbl"
    params: nhmmer=config.get("nhmmer", "nhmmer")
    threads: 32
    shell: """
        {params.nhmmer} --cpu {threads} \
            --tblout {output.tbl} {input.model} {input.ar} \
        > /dev/null
    """

rule candidate_densities:
    input: sam="{path}/{prime}AC.sam", motifs="{path}/{prime}AC-motifs.txt"
    output: dat="{path}/{prime}AC-densities.dat.gz"
    params: window_size=120, revcomp=True
    threads: 16
    run:
        with open(input.motifs) as motif_data:
            motifs = {line.strip() for line in motif_data}
        if params.revcomp:
            motifs |= {motif_revcomp(motif) for motif in set(motifs)}
        with gzopen(output.dat, mode="wt") as dat:
            for i, motif in enumerate(motifs):
                no_print_header = (i != 0)
                kmerscanner.main(
                    bams=[input.sam], num_reads=None,
                    motif=motif, window_size=params.window_size,
                    head_test=None, tail_test=None, cutoff=None,
                    jobs=threads, no_print_header=no_print_header, file=dat
                )

rule densityplot:
    input: dat="{path}/{prime}AC-densities.dat.gz"
    output: pdf="{path}/{prime}AC-densities.pdf"
    params: bin_size=100
    run:
        with open(output.pdf, mode="wb") as pdf:
            densityplot.main(
                dat=input.dat, bin_size=params.bin_size, file=pdf,
                title="{} {}AC".format(wildcards.dataset, wildcards.prime)
            )

from sys import stdout, stderr
from tempfile import TemporaryDirectory
from pysam import AlignmentFile, FastxFile
from os import path
from edgecaselib.util import get_executable
from edgecaselib.formats import filter_bam
from tqdm import tqdm
from subprocess import check_output
from pandas import read_csv, concat, merge
from numpy import unique
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def interpret_args(fmt, jellyfish):
    """Parse and check arguments"""
    if fmt == "sam":
        manager = AlignmentFile
    elif fmt == "fastx":
        manager = FastxFile
    else:
        raise ValueError("Unsupported --fmt: '{}'".format(fmt))
    return manager, get_executable("jellyfish", jellyfish)


def convert_input(bam, manager, tempdir, samfilters):
    """Convert BAM to fasta for MEME"""
    fasta = path.join(tempdir, "input.fa")
    with manager(bam) as alignment, open(fasta, mode="wt") as fasta_handle:
        for entry in filter_bam(alignment, samfilters, "SAM/BAM -> FASTA"):
            entry_str = ">{}\n{}".format(entry.qname, entry.query_sequence)
            print(entry_str, file=fasta_handle)
    return fasta


def find_repeats(sequencefile, min_k, max_k, no_context, jellyfish, jobs, tempdir):
    """Find all repeats in sequencefile"""
    per_k_reports = []
    for k in tqdm(range(min_k, max_k+1), desc="Sweeping lengths"):
        db = path.join(tempdir, "{}.db".format(k))
        if no_context:
            search_k = str(k)
        else:
            # instead of looking for e.g. TTAGGG, look for TTAGGGTTAGGG,
            # this implies repeat context and larger motifs like TTAGGGA do
            # not confound it:
            search_k = str(k*2)
        check_output([
            jellyfish, "count", "-t", str(jobs), "-s", "2G", "-L", "0",
            "-m", search_k, "-o", db, sequencefile
        ])
        tsv = path.join(tempdir, "{}.tsv".format(k))
        check_output([
            jellyfish, "dump", "-c", "-t", "-L", "0",
            "-o", tsv, db
        ])
        k_report = read_csv(tsv, sep="\t", names=["kmer", "count"])
        if not no_context: # for downstream analyses, roll back to single motif
            doubles_indexer = k_report["kmer"].apply(
                lambda kmer:kmer[:k]==kmer[k:]
            )
            k_report = k_report[doubles_indexer]
            k_report["kmer"] = k_report["kmer"].apply(lambda kmer:kmer[:k])
        k_report["abundance"] = k_report["count"] / k_report["count"].sum()
        k_report["length"] = k
        per_k_reports.append(k_report)
    return concat(per_k_reports, axis=0)


def lowest_alpha_inversion(kmer):
    """Get alphabetically lowest inversion of kmer (e.g., for TTAGGG will return AGGGTT)"""
    return min(kmer[i:]+kmer[:i] for i in range(len(kmer)))


def get_motifs_fisher(single_length_report):
    """Analyze repeat enrichment given the same motif length"""
    lengths = unique(single_length_report["length"].values)
    if len(lengths) != 1:
        raise ValueError("`get_motifs_fisher`: multiple lengths found")
    else:
        k = lengths[0]
    message = "\rCalculating enrichment for k={}... ".format(lengths[0])
    print(message, end="", file=stderr, flush=True)
    total_count = single_length_report["count"].sum()
    N = 4**k
    median_count = single_length_report["count"].median()
    high_indexer = (single_length_report["count"]>=median_count)
    fishers = single_length_report[high_indexer].copy()
    fishers["p"] = fishers["count"].apply(
        lambda c: fisher_exact([[c, total_count-c], [1, N]])[1]
    )
    return fishers


def analyze_repeats(full_report, adj="bonferroni"):
    """Analyze repeat enrichment for multiple lengths and apply multiple testing adjustment"""
    fishers = concat([
        get_motifs_fisher(
            full_report[full_report["length"]==length]
        )
        for length in unique(full_report["length"].values)
    ])
    fishers["motif"] = fishers["kmer"].apply(lowest_alpha_inversion)
    motif_p = fishers[["motif", "p"]].groupby("motif", as_index=False).max()
    motif_p["p_adjusted"] = multipletests(motif_p["p"], method=adj)[1]
    fishers = merge(
        fishers, motif_p[["motif", "p_adjusted"]], on="motif", how="outer"
    )
    ktl = fishers[["count", "abundance", "motif", "length", "p_adjusted"]]
    ktl_grouper = ["motif", "length", "p_adjusted"]
    return ktl.groupby(ktl_grouper, as_index=False).sum()


def format_analysis(filtered_analysis, max_motifs):
    """Make dataframe prettier"""
    formatted_analysis = filtered_analysis.sort_values(
        by=["abundance", "p_adjusted"], ascending=[False, True]
    )
    formatted_analysis = formatted_analysis[
        ["motif", "length", "count", "abundance", "p_adjusted"]
    ]
    formatted_analysis.columns = [
        "#motif", "length", "count", "abundance", "p_adjusted"
    ]
    return formatted_analysis[:max_motifs]


def main(sequencefile, fmt, flags, flags_any, flag_filter, min_quality, min_k, max_k, max_motifs, max_p_adjusted, no_context, jellyfish, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    manager, jellyfish = interpret_args(fmt, jellyfish)
    with TemporaryDirectory() as tempdir:
        if manager == AlignmentFile: # will need to convert SAM to fastx
            samfilters = [flags, flags_any, flag_filter, min_quality]
            sequencefile = convert_input(
                sequencefile, manager, tempdir, samfilters
            )
        full_report = find_repeats(
            sequencefile, min_k, max_k, no_context, jellyfish, jobs, tempdir
        )
    analysis = analyze_repeats(full_report)
    filtered_analysis = analysis[analysis["p_adjusted"]<max_p_adjusted]
    formatted_analysis = format_analysis(filtered_analysis, max_motifs)
    formatted_analysis.to_csv(file, sep="\t", index=False)
    print("Done", file=stderr, flush=True)

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
from scipy.stats import chi2_contingency
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


def find_repeats(sequencefile, min_k, max_k, jellyfish, jobs, tempdir):
    """Find all repeats in sequencefile"""
    per_k_reports = []
    for k in tqdm(range(min_k, max_k+1), desc="Sweeping lengths"):
        db = path.join(tempdir, "{}.db".format(k))
        check_output([
            jellyfish, "count", "-t", str(jobs), "-s", "2G", "-L", "0",
            "-m", str(k), "-o", db, sequencefile
        ])
        tsv = path.join(tempdir, "{}.tsv".format(k))
        check_output([
            jellyfish, "dump", "-c", "-t", "-L", "0",
            "-o", tsv, db
        ])
        k_report = read_csv(tsv, sep="\t", names=["kmer", "count"])
        k_report["length"] = k
        per_k_reports.append(k_report)
    return concat(per_k_reports, axis=0)


def lowest_alpha_inversion(kmer):
    return min(kmer[i:]+kmer[:i] for i in range(len(kmer)))


def analyze_repeats(full_report, adj="fdr_bh"):
    """Analyze repeat enrichment"""
    lengths = unique(full_report["length"].values)
    if adj is None:
        if len(lengths) != 1:
            message = "`{}` can only be directly called with the `adj` argument"
            raise NotImplementedError(message.format("analyze_repeats"))
        message = "\rCalculating for k={}... ".format(lengths[0])
        print(message, end="", file=stderr, flush=True)
        total_count = full_report["count"].sum()
        N = 4**lengths[0]
        median_count = full_report["count"].median()
        chi2s = full_report[full_report["count"]>=median_count].copy()
        chi2s["p"] = chi2s["count"].apply(
            lambda c: chi2_contingency([[c, total_count-c], [1, N]])[1]
        )
        return chi2s
    else:
        chi2s = concat([
            analyze_repeats(
                full_report[full_report["length"]==length], adj=None
            )
            for length in lengths
        ])
        chi2s["motif"] = chi2s["kmer"].apply(lowest_alpha_inversion)
        motif_p = chi2s[["motif", "p"]].groupby("motif", as_index=False).max()
        motif_p["p_adjusted"] = multipletests(motif_p["p"], method="fdr_bh")[1]
        chi2s = merge(
            chi2s, motif_p[["motif", "p_adjusted"]],
            on="motif", how="outer"
        )
        ktl = chi2s[["count", "motif", "length", "p_adjusted"]]
        ktl_grouper = ["motif", "length", "p_adjusted"]
        return ktl.groupby(ktl_grouper, as_index=False).sum()


def format_analysis(filtered_analysis):
    formatted_analysis = filtered_analysis.sort_values(
        by=["count", "p_adjusted"], ascending=[False, True]
    )
    formatted_analysis = formatted_analysis[
        ["motif", "length", "count", "p_adjusted"]
    ]
    formatted_analysis.columns = [
        "#motif", "length", "count", "p_adjusted(length)"
    ]
    return formatted_analysis


def main(sequencefile, fmt, flags, flags_any, flag_filter, min_quality, min_k, max_k, max_p_adjusted, jellyfish, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    manager, jellyfish = interpret_args(fmt, jellyfish)
    with TemporaryDirectory() as tempdir:
        if manager == AlignmentFile: # will need to convert SAM to fastx
            samfilters = [flags, flags_any, flag_filter, min_quality]
            sequencefile = convert_input(
                sequencefile, manager, tempdir, samfilters
            )
        full_report = find_repeats(
            sequencefile, min_k, max_k, jellyfish, jobs, tempdir
        )
    analysis = analyze_repeats(full_report)
    filtered_analysis = analysis[analysis["p_adjusted"]<max_p_adjusted]
    formatted_analysis = format_analysis(filtered_analysis)
    formatted_analysis.to_csv(file, sep="\t", index=False)
    print("Done", file=stderr, flush=True)

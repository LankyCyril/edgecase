from sys import stdout
from tempfile import TemporaryDirectory
from pysam import AlignmentFile, FastxFile
from re import finditer, IGNORECASE
from os import path
from edgecaselib.util import get_executable, progressbar
from edgecaselib.formats import filter_bam
from functools import lru_cache
from subprocess import check_output
from pandas import read_csv, concat
from numpy import unique
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


__doc__ = """edgeCase repeatfinder: de novo repeat discovery

Usage: {0} repeatfinder [-m integer] [-M integer] [-r integer] [-P float]
       {1}              [--jellyfish filename] [--jellyfish-hash-size string]
       {1}              [-n integer] [-j integer] [-f flagspec] [-g flagspec]
       {1}              [-F flagspec] [-q integer] [--fmt string]
       {1}              [--collapse-reverse-complement] <sequencefile>

Output:
    TSV-formatted file with statistics describing enriched motifs

Positional arguments:
    <sequencefile>                      name of input BAM/SAM/FASTA/FASTQ file

Options:
    --fmt sam|fastx                     format of input file [default: sam]
    -m, --min-k [integer]               smallest target repeat length [default: 4]
    -M, --max-k [integer]               largest target repeat length [default: 16]
    -r, --min-repeats [integer]         minimum number of consecutive repeats [default: 2]
    -P, --max-p-adjusted [float]        cutoff adjusted p-value [default: .05]
    --jellyfish [filename]              jellyfish binary (unless in $PATH)
    -s, --jellyfish-hash-size [string]  jellyfish initial hash size [default: 2G]
    -n, --max-motifs [integer]          maximum number of motifs to report
    -j, --jobs [integer]                number of jellyfish jobs (parallel threads) [default: 1]
    -C, --collapse-reverse-complement   collapse counts of reverse complement motifs

Input filtering options:
    -f, --flags [flagspec]              process only entries with all these sam flags present [default: 0]
    -g, --flags-any [flagspec]          process only entries with any of these sam flags present [default: 65535]
    -F, --flag-filter [flagspec]        process only entries with none of these sam flags present [default: 0]
    -q, --min-quality [integer]         process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda min_quality: None if (min_quality is None) else int(min_quality),
    lambda min_k: int(min_k),
    lambda max_k: int(max_k),
    lambda min_repeats: int(min_repeats),
    lambda max_p_adjusted: float(max_p_adjusted),
    lambda jobs: int(jobs),
]

__docopt_tests__ = {
    lambda min_k, max_k: 0 < min_k < max_k: "not satisfied: 0 < m < M",
    lambda min_repeats: min_repeats > 0: "--min-repeats must be integer > 0",
}


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
    """Convert BAM to fasta; count bases"""
    fasta = path.join(tempdir, "input.fa")
    base_count = 0
    with manager(bam) as alignment, open(fasta, mode="wt") as fasta_handle:
        for entry in filter_bam(alignment, samfilters, "SAM/BAM -> FASTA"):
            entry_str = ">{}\n{}".format(entry.qname, entry.query_sequence)
            base_count += len(entry.query_sequence)
            print(entry_str, file=fasta_handle)
    return fasta, base_count


def count_fastx_bases(sequencefile, pattern=r'[acgt]', flags=IGNORECASE, desc="Counting input bases"):
    """Count bases in FASTX file"""
    with FastxFile(sequencefile) as fastx:
        return sum(
            sum(1 for _ in finditer(pattern, entry.sequence, flags=flags))
            for entry in progressbar(fastx, desc=desc, unit="read")
        )


def find_repeats(sequencefile, min_k, max_k, min_repeats, base_count, jellyfish, jellyfish_hash_size, collapse_reverse_complement, jobs, tempdir):
    """Find all repeats in sequencefile"""
    per_k_reports = []
    k_iterator = progressbar(
        range(min_k, max_k+1), desc="Sweeping lengths", unit="k",
    )
    for k in k_iterator:
        db = path.join(tempdir, "{}.db".format(k))
        jellyfish_count_options = [
            jellyfish, "count", "-t", str(jobs), "-s", jellyfish_hash_size,
            "-L", "0", "-m", str(k * min_repeats)
        ]
        if collapse_reverse_complement:
            jellyfish_count_options += ["-C"]
        check_output(jellyfish_count_options + ["-o", db, sequencefile])
        tsv = path.join(tempdir, "{}.tsv".format(k))
        check_output([
            jellyfish, "dump", "-c", "-t", "-L", "0",
            "-o", tsv, db,
        ])
        k_report = read_csv(tsv, sep="\t", names=["kmer", "count"])
        if len(k_report) == 0:
            return None
        repeats_indexer = k_report["kmer"].apply(
            lambda kmer: kmer[:k] * min_repeats == kmer
        )
        k_report = k_report[repeats_indexer]
        k_report["kmer"] = k_report["kmer"].apply(lambda kmer:kmer[:k])
        k_report["abundance"] = k_report["count"] / base_count
        k_report["length"] = k
        per_k_reports.append(k_report)
    return concat(per_k_reports, axis=0)


@lru_cache(maxsize=None)
def lowest_alpha_inversion(kmer):
    """Get alphabetically lowest inversion of kmer (e.g., for TTAGGG will return AGGGTT)"""
    return min(kmer[i:]+kmer[:i] for i in range(len(kmer)))


def custom_alpha_inversion(motif):
    """Get inversion of motif that looks closest to canonical (e.g., for AGGGTTC will return TTCAGGG)"""
    a, c, g, t = [motif.count(letter) for letter in "ACGT"]
    is_g = (g > c) or ((g == c) and (t > a))
    if is_g:
        if t > 0:
            i = motif.find("T")
        else:
            i = motif.find("G")
    else:
        if c > 0:
            i = motif.find("C")
        else:
            i = motif.find("A")
    if i == -1:
        return min(motif[i:]+motif[:i] for i in range(len(motif)))
    else:
        return motif[i:] + motif[:i]


def safe_fisher_exact(count, total_candidate_count, med, total_background_count):
    """Same as running fisher_exact, but returns p=1 if encounters NaNs (i.e., when not enough data to process test)"""
    try:
        return fisher_exact([
            [count, total_candidate_count - count],
            [med, total_background_count - med],
        ])[1]
    except ValueError:
        return 1


def get_motifs_fisher(single_length_report):
    """Analyze repeat enrichment given the same motif length"""
    lengths = unique(single_length_report["length"].values)
    if len(lengths) != 1:
        raise ValueError("`get_motifs_fisher`: multiple lengths found")
    fishery = single_length_report.copy()
    fishery["motif"] = fishery["kmer"].apply(lowest_alpha_inversion)
    fishery_groupby = fishery[["motif", "count", "abundance"]].groupby(
        "motif", as_index=False,
    )
    fishery = fishery_groupby.sum()
    iqr = fishery["count"].quantile(.75) - fishery["count"].quantile(.25)
    whisker = fishery["count"].quantile(.75) + 1.5 * iqr
    candidates, background = (
        fishery[fishery["count"]>=whisker].copy(),
        fishery[fishery["count"]<whisker].copy(),
    )
    total_candidate_count, total_background_count = (
        candidates["count"].sum(), background["count"].sum(),
    )
    med = fishery["count"].median()
    candidates["p"] = candidates["count"].apply(
        lambda count: safe_fisher_exact(
            count, total_candidate_count,
            med, total_background_count,
        )
    )
    candidates["length"] = lengths[0]
    return candidates


def analyze_repeats(full_report, adj="bonferroni"):
    """Analyze repeat enrichment for multiple lengths and apply multiple testing adjustment"""
    candidates = concat([
        get_motifs_fisher(full_report[full_report["length"]==length])
        for length in progressbar(
            unique(full_report["length"].values), unit="k",
            desc="Calculating enrichment",
        )
    ])
    candidates["p_adjusted"] = multipletests(candidates["p"], method=adj)[1]
    return candidates[
        ["motif", "length", "count", "abundance", "p", "p_adjusted"]
    ]


def coerce_and_filter_report(analysis, max_p_adjusted):
    """Collapse functionally synonymous entries like TATATA and TATA"""
    motif_mapper, mapped_motifs = {}, set()
    for motif in analysis.sort_values(by="length")["motif"]:
        length_indexer = analysis["length"]>len(motif)
        for larger_motif in analysis.loc[length_indexer, "motif"]:
            if (larger_motif not in mapped_motifs) and (motif in larger_motif):
                if motif*len(larger_motif) == larger_motif*len(motif):
                    if motif not in motif_mapper:
                        motif_mapper[motif] = {motif}
                    motif_mapper[motif].add(larger_motif)
                    mapped_motifs.add(larger_motif)
    synonyms_to_keep = set()
    for synonyms in motif_mapper.values():
        synonym_data = analysis[
            analysis["motif"].isin(synonyms) &
            (analysis["p_adjusted"]<max_p_adjusted)
        ]
        if len(synonym_data):
            synonyms_to_keep.add(
                synonym_data.sort_values(
                    by="abundance", ascending=False
                ).iloc[0, 0]
            )
    synonyms_to_remove = (
        set.union(set(), *motif_mapper.values()) - synonyms_to_keep
    )
    return analysis[
        (~analysis["motif"].isin(synonyms_to_remove)) &
        (analysis["p_adjusted"]<max_p_adjusted)
    ].copy()


def coerce_to_monomer(motif, min_k):
    """Coerce motif to monomer, e.g. TATA -> TA, CAT -> CAT; this can be used to find functionally synonymous entries too"""
    n = len(motif)
    for i in range(min_k, int(n/2)+1):
        q, r = divmod(n, i)
        if r == 0:
            if motif[:i]*q == motif:
                return motif[:i]
    else:
        return motif


def format_analysis(filtered_analysis, min_k, max_motifs):
    """Make dataframe prettier"""
    filtered_analysis["motif"] = filtered_analysis["motif"].apply(
        custom_alpha_inversion,
    )
    filtered_analysis["monomer"] = filtered_analysis["motif"].apply(
        lambda motif: coerce_to_monomer(motif, min_k=min_k),
    )
    formatted_analysis = filtered_analysis.sort_values(
        by=["abundance", "p_adjusted"], ascending=[False, True],
    )
    formatted_analysis = formatted_analysis[
        ["monomer", "motif", "length", "count", "abundance", "p", "p_adjusted"]
    ]
    formatted_analysis.columns = [
        "#monomer", "motif", "length", "count", "abundance", "p", "p_adjusted",
    ]
    if max_motifs is None:
        return formatted_analysis
    else:
        return formatted_analysis[:max_motifs]


def main(sequencefile, fmt, flags, flags_any, flag_filter, min_quality, min_k, max_k, min_repeats, max_motifs, max_p_adjusted, jellyfish, jellyfish_hash_size, collapse_reverse_complement, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    manager, jellyfish = interpret_args(fmt, jellyfish)
    with TemporaryDirectory() as tempdir:
        if manager == AlignmentFile: # will need to convert SAM to fastx
            samfilters = [flags, flags_any, flag_filter, min_quality]
            sequencefile, base_count = convert_input(
                sequencefile, manager, tempdir, samfilters,
            )
        else:
            base_count = count_fastx_bases(sequencefile)
        full_report = find_repeats(
            sequencefile, min_k, max_k, min_repeats, base_count,
            jellyfish, jellyfish_hash_size, collapse_reverse_complement,
            jobs, tempdir,
        )
    if full_report is None:
        columns = [
            "#monomer", "motif", "length", "count", "abundance",
            "p", "p_adjusted",
        ]
        print(*columns, sep="\t", file=file)
    else:
        analysis = analyze_repeats(full_report)
        filtered_analysis = coerce_and_filter_report(analysis, max_p_adjusted)
        formatted_analysis = format_analysis(
            filtered_analysis, min_k, max_motifs,
        )
        formatted_analysis.to_csv(file, sep="\t", index=False)

from sys import stdout
from tempfile import TemporaryDirectory
from pysam import AlignmentFile, FastxFile
from os import path
from edgecaselib.util import get_executable, progressbar
from edgecaselib.formats import filter_bam
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
    """Convert BAM to fasta; count bases"""
    fasta = path.join(tempdir, "input.fa")
    base_count = 0
    with manager(bam) as alignment, open(fasta, mode="wt") as fasta_handle:
        for entry in filter_bam(alignment, samfilters, "SAM/BAM -> FASTA"):
            entry_str = ">{}\n{}".format(entry.qname, entry.query_sequence)
            base_count += len(entry.query_sequence)
            print(entry_str, file=fasta_handle)
    return fasta, base_count


def count_fastx_bases(sequencefile):
    """Count bases in FASTX file"""
    with FastxFile(sequencefile) as fastx:
        return sum(
            len(entry.sequence) for entry
            in progressbar(fastx, desc="Counting input bases", unit="read")
        )


def find_repeats(sequencefile, min_k, max_k, base_count, no_context, jellyfish, jobs, tempdir):
    """Find all repeats in sequencefile"""
    per_k_reports = []
    k_iterator = progressbar(
        range(min_k, max_k+1), desc="Sweeping lengths", unit="k"
    )
    for k in k_iterator:
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
        k_report["abundance"] = k_report["count"] / base_count
        k_report["length"] = k
        per_k_reports.append(k_report)
    return concat(per_k_reports, axis=0)


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
    elif (not is_g) and (t == a):
        i = -1
    else:
        if c > 0:
            i = motif.find("C")
        else:
            i = motif.find("A")
    if i == -1:
        return min(motif[i:]+motif[:i] for i in range(len(motif)))
    else:
        return motif[i:] + motif[:i]


def get_motifs_fisher(single_length_report):
    """Analyze repeat enrichment given the same motif length"""
    lengths = unique(single_length_report["length"].values)
    if len(lengths) != 1:
        raise ValueError("`get_motifs_fisher`: multiple lengths found")
    else:
        k = lengths[0]
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
        for length in progressbar(
            unique(full_report["length"].values), unit="k",
            desc="Calculating enrichment"
        )
    ])
    fishers["motif"] = fishers["kmer"].apply(lowest_alpha_inversion)
    motif_p = fishers[["motif", "p"]].groupby("motif", as_index=False).max()
    motif_p["p_adjusted"] = multipletests(motif_p["p"], method=adj)[1]
    fishers = merge(
        fishers.drop(columns="p"), motif_p[["motif", "p", "p_adjusted"]],
        on="motif", how="outer"
    )
    ktl = fishers[["count", "abundance", "motif", "length", "p", "p_adjusted"]]
    ktl_grouper = ["motif", "length", "p", "p_adjusted"]
    return ktl.groupby(ktl_grouper, as_index=False).sum()


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


def coerce_to_monomer(motif):
    """Coerce motif to monomer, e.g. TATA -> TA, CAT -> CAT; this can be used to find functionally synonymous entries too"""
    n = len(motif)
    for i in range(1, int(n/2)+1):
        q, r = divmod(n, i)
        if r == 0:
            if motif[:i]*q == motif:
                return motif[:i]
    else:
        return motif


def format_analysis(filtered_analysis, max_motifs):
    """Make dataframe prettier"""
    filtered_analysis["motif"] = filtered_analysis["motif"].apply(
        custom_alpha_inversion
    )
    filtered_analysis["monomer"] = filtered_analysis["motif"].apply(
        coerce_to_monomer
    )
    formatted_analysis = filtered_analysis.sort_values(
        by=["abundance", "p_adjusted"], ascending=[False, True]
    )
    formatted_analysis = formatted_analysis[
        ["monomer", "motif", "length", "count", "abundance", "p", "p_adjusted"]
    ]
    formatted_analysis.columns = [
        "#monomer", "motif", "length", "count", "abundance", "p", "p_adjusted"
    ]
    if max_motifs is None:
        return formatted_analysis
    else:
        return formatted_analysis[:max_motifs]


def main(sequencefile, fmt, flags, flags_any, flag_filter, min_quality, min_k, max_k, max_motifs, max_p_adjusted, no_context, jellyfish, jobs=1, file=stdout, **kwargs):
    # parse arguments:
    manager, jellyfish = interpret_args(fmt, jellyfish)
    with TemporaryDirectory() as tempdir:
        if manager == AlignmentFile: # will need to convert SAM to fastx
            samfilters = [flags, flags_any, flag_filter, min_quality]
            sequencefile, base_count = convert_input(
                sequencefile, manager, tempdir, samfilters
            )
        else:
            base_count = count_fastx_bases(sequencefile)
        full_report = find_repeats(
            sequencefile, min_k, max_k, base_count,
            no_context, jellyfish, jobs, tempdir
        )
    analysis = analyze_repeats(full_report)
    filtered_analysis = coerce_and_filter_report(analysis, max_p_adjusted)
    formatted_analysis = format_analysis(filtered_analysis, max_motifs)
    formatted_analysis.to_csv(file, sep="\t", index=False)

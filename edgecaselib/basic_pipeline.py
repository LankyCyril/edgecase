from os import path, walk
from numpy import inf
from edgecaselib import tailpuller, tailchopper, repeatfinder
from edgecaselib import kmerscanner, densityplot
from edgecaselib.formats import EmptyKmerscanError
from gzip import open as gzopen

__doc__ = """edgeCase basic pipeline: select reads, find enriched motifs, plot

Usage: {0} basic-pipeline -x filename -o dirname [--force] [-j integer]
       {1}                [-m integer] [-n integer] [--target targetspec]
       {1}                [--min-k integer] [--max-k integer]
       {1}                [--jellyfish filename] [--jellyfish-hash-size string]
       {1}                [--max-p-adjusted float] [--window-size integer]
       {1}                [--palette palettespec] [-q integer] <bam>

Output (in --output-dir):
    * tailpuller.sam                 candidate reads
    * tailchopper.sam                telomeric regions
    * repeatfinder-p_arm.tsv         overrepresented motifs on the p arm
    * repeatfinder-q_arm.tsv         overrepresented motifs on the q arm
    * kmerscanner-p_arm.dat.gz       motif densities on the p arm
    * kmerscanner-q_arm.dat.gz       motif densities on the q arm
    * densityplot-p_arm.pdf          plot of motif densities on the p arm
    * densityplot-q_arm.pdf          plot of motif densities on the q arm

Positional arguments:
    <bam>                            name of input BAM/SAM file; must have a .bai index

Required options:
    -x, --index [filename]           location of the reference .ecx index
    -o, --output-dir [dirname]       name of the output directory (must exist)

Options:
    --force                          force overwrite files in --output-dir
    -j, --jobs [integer]             number of parallel jobs (for jellyfish and kmerscanner) [default: 1]
    -m, --max-read-length [integer]  maximum read length to consider when selecting lookup regions
    -n, --max-motifs [integer]       maximum number of motifs to report [default: 4]
    --target [targetspec]            an ECX flag for heads/tails [default: tract_anchor]
    --min-k [integer]                smallest target repeat length [default: 4]
    --max-k [integer]                largest target repeat length [default: 16]
    --max-p-adjusted [float]         cutoff adjusted p-value [default: .05]
    --jellyfish [filename]           jellyfish binary (unless in $PATH)
    --jellyfish-hash-size [string]   jellyfish initial hash size [default: 2G]
    --window-size [integer]          size of the window (for kmerscanner and densityplot) [default: 100]
    --palette [palettespec]          custom palette for plotting motifs

Input filtering options:
    -q, --min-quality [integer]      process only entries with this MAPQ or higher [default: 0]
"""

__docopt_converters__ = [
    lambda jobs: int(jobs),
    lambda max_read_length:
        inf if (max_read_length is None) else int(max_read_length),
    lambda max_motifs: int(max_motifs),
    lambda min_k: int(min_k),
    lambda max_k: int(max_k),
    lambda max_p_adjusted: float(max_p_adjusted),
    lambda window_size: int(window_size),
    lambda min_quality: None if (min_quality is None) else int(min_quality),
]

__docopt_tests__ = {
    lambda bam:
        path.isfile(bam + ".bai"): "BAM index ({}.bai) not found",
    lambda output_dir:
        path.isdir(output_dir): "{} does not exist or is not a directory",
    lambda output_dir, force:
        force or (len(next(walk(output_dir))[2]) == 0):
            "{} already has files in it",
    lambda max_read_length:
        max_read_length > 0: "--max-read-length below 0",
    lambda target:
        target in {"ucsc_mask_anchor", "fork", "tract_anchor"}:
            "unknown value of --target",
    lambda min_k, max_k:
        0 < min_k < max_k: "not satisfied: 0 < m < M",
}


def main(bam, index, output_dir, jobs, max_read_length, max_motifs, target, min_k, max_k, max_p_adjusted, jellyfish, jellyfish_hash_size, window_size, palette, min_quality, **kwargs):
    """basic pipeline: select reads, find enriched motifs, plot"""
    get_filename = lambda fn: path.join(output_dir, fn)
    tailpuller_sam = get_filename("tailpuller.sam")
    with open(tailpuller_sam, mode="wt") as sam:
        tailpuller.main(
            bam, index, flags=0, flags_any=65535, flag_filter=3844,
            min_quality=min_quality, max_read_length=max_read_length, file=sam,
        )
    tailchopper_sam = get_filename("tailchopper.sam")
    with open(tailchopper_sam, mode="wt") as sam:
        tailchopper.main(
            tailpuller_sam, index, target=target, flags=target, flags_any=65535,
            flag_filter=3844, min_quality=min_quality, file=sam,
        )
    for arm, f, F in [("p", target, "is_q|3840"), ("q", target+"|is_q", 3840)]:
        repeatfinder_tsv = get_filename("repeatfinder-{}_arm.tsv".format(arm))
        with open(repeatfinder_tsv, mode="wt") as tsv:
            repeatfinder.main(
                tailchopper_sam, fmt="sam",
                flags=f, flags_any=65535, flag_filter=F,
                min_quality=min_quality, min_k=min_k, max_k=max_k,
                max_motifs=max_motifs, max_p_adjusted=max_p_adjusted,
                no_context=False, jellyfish=jellyfish,
                jellyfish_hash_size=jellyfish_hash_size, jobs=jobs, file=tsv,
            )
        kmerscanner_dat = get_filename("kmerscanner-{}_arm.dat.gz".format(arm))
        with gzopen(kmerscanner_dat, mode="wt") as dat:
            kmerscanner.main(
                tailpuller_sam, flags=f, flags_any=65535, flag_filter=F,
                min_quality=min_quality, motif_file=repeatfinder_tsv,
                head_test=None, tail_test=None, cutoff=None, num_reads=None,
                window_size=window_size, jobs=jobs, file=dat
            )
        densityplot_pdf = get_filename("densityplot-{}_arm.pdf".format(arm))
        with open(densityplot_pdf, mode="wb") as pdf:
            try:
                densityplot.main(
                    kmerscanner_dat, gzipped=True, index=index, flags=f,
                    flags_any=65535, flag_filter=F, min_quality=min_quality,
                    bin_size=window_size, exploded=False, zoomed_in=False,
                    palette=palette, title=None, file=pdf
                )
            except EmptyKmerscanError:
                pass
    return 0

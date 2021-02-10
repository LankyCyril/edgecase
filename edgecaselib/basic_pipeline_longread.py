from os import path, walk
from numpy import inf
from edgecaselib import tailpuller, tailchopper, repeatfinder
from edgecaselib import kmerscanner, densityplot
from edgecaselib.formats import EmptyKmerscanError
from gzip import open as gzopen

__doc__ = """edgeCase basic pipeline: select reads, find enriched motifs, plot

Usage: {0} basic-pipeline-longread -x filename -o dirname [--force] [-j integer]
       {1}       [-M integer] [-m integer] [--target targetspec] [-q integer]
       {1}       [-n integer] [--min-k integer] [--max-k integer]
       {1}       [--jellyfish filename] [--jellyfish-hash-size string]
       {1}       [--max-p-adjusted float] [--bin-size integer]
       {1}       [--n-boot integer] [--palette palettespec]
       {1}       [--title string] <bam>

Output (in --output-dir):
    * tailpuller.sam                       candidate reads
    * tailchopper.sam                      telomeric regions
    * repeatfinder-p_arm.tsv               overrepresented motifs on the p arm
    * repeatfinder-q_arm.tsv               overrepresented motifs on the q arm
    * kmerscanner-p_arm.dat.gz             motif densities on the p arm
    * kmerscanner-q_arm.dat.gz             motif densities on the q arm
    * densityplot-p_arm.pdf                plot of motif densities on the p arm
    * densityplot-q_arm.pdf                plot of motif densities on the q arm

Positional arguments:
    <bam>                                  name of input BAM/SAM file; must have a .bai index

Required options:
    -x, --index [filename]                 location of the reference .ecx index
    -o, --output-dir [dirname]             name of the output directory (must exist)

Options:
    --force                                force overwrite files in --output-dir
    -j, --jobs [integer]                   number of parallel jobs (for jellyfish and kmerscanner) [default: 1]
    -M, --max-read-length [integer]        maximum read length to consider when selecting lookup regions *
    --min-map-overlap [integer]            minimum overlap of reference to consider read as mapped [default: 1] **
    -m, --min-candidate-overlap [integer]  minimum overlap of subtelomere to consider read as candidate [default: 1] ***
    --target [targetspec]                  an ECX flag for heads/tails [default: tract_anchor]
    -n, --max-motifs [integer]             maximum number of motifs to report [default: 4]
    --min-k [integer]                      smallest target repeat length [default: 4]
    --max-k [integer]                      largest target repeat length [default: 16]
    --max-p-adjusted [float]               cutoff adjusted p-value [default: .05]
    --jellyfish [filename]                 jellyfish binary (unless in $PATH)
    --jellyfish-hash-size [string]         jellyfish initial hash size [default: 2G]
    --bin-size [integer]                   size of the window (for kmerscanner and densityplot) [default: 10]
    --n-boot [integer]                     number of bootstrap iterations for plotting [default: 1000]
    --palette [palettespec]                custom palette for plotting motifs
    --title [string]                       title for the density plots

Input filtering options:
    -q, --min-quality [integer]      process only entries with this MAPQ or higher [default: 0] ****

Notes:
   * Suggested value of --max-read-length for PacBio: 200000;
     if not specified, will assume +infinity (will be slow).
  ** Suggested value of --min-map-overlap for PacBio: 500;
     if not specified, will assume 1.
 *** Suggested value of --min-candidate-overlap for PacBio: 3000;
     if not specified, will assume 1.
**** Depending on the aligner used, MAPQ of secondary reads may have been set to
     zero regardless of real mapping quality; use this filtering option with
     caution.
"""

__docopt_converters__ = [
    lambda jobs: int(jobs),
    lambda max_read_length:
        inf if (max_read_length is None) else int(max_read_length),
    lambda min_map_overlap:
        1 if (min_map_overlap is None) else int(min_map_overlap),
    lambda min_candidate_overlap:
        1 if (min_candidate_overlap is None) else int(min_candidate_overlap),
    lambda max_motifs: int(max_motifs),
    lambda min_k: int(min_k),
    lambda max_k: int(max_k),
    lambda max_p_adjusted: float(max_p_adjusted),
    lambda bin_size: int(bin_size),
    lambda n_boot: int(n_boot),
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
        target in {"mask_anchor", "fork", "tract_anchor"}:
            "unsupported value of --target",
    lambda min_k, max_k:
        0 < min_k < max_k: "not satisfied: 0 < m < M",
}


get_tailpuller_kws = lambda min_quality, max_read_length, min_map_overlap, min_candidate_overlap, target, file: dict(
    max_read_length=max_read_length, target={target},
    min_map_overlap=min_map_overlap,
    min_candidate_overlap=min_candidate_overlap,
    flags=[0], flag_filter=[0], min_quality=min_quality,
    file=file,
)


get_tailchopper_kws = lambda target, min_quality, file: dict(
    target=target, flags=[target], flag_filter=[0],
    min_quality=min_quality, file=file,
)


get_repeatfinder_kws = lambda f, F, q, max_motifs, min_k, max_k, max_p_adjusted, jellyfish, jellyfish_hash_size, jobs, file: dict(
    fmt="sam", collapse_reverse_complement=False, flags=f, flag_filter=F,
    min_quality=q, min_repeats=2, max_motifs=max_motifs,
    min_k=min_k, max_k=max_k, no_context=False, max_p_adjusted=max_p_adjusted,
    jellyfish=jellyfish, jellyfish_hash_size=jellyfish_hash_size,
    jobs=jobs, file=file,
)


get_kmerscanner_kws = lambda f, F, q, motif_file, bin_size, jobs, file: dict(
    fmt="sam", flags=f, flag_filter=F,
    min_quality=q, motif_file=motif_file, bin_size=bin_size,
    head_test=None, tail_test=None, cutoff=None, num_reads=None,
    jobs=jobs, file=file,
)


get_densityplot_kws = lambda f, F, q, bin_size, n_boot, palette, title, arm, file: dict(
    gzipped=True, flags=f, flag_filter=F, min_quality=q,
    exploded=False, zoomed_in=False,
    bin_size=bin_size, n_boot=n_boot,
    title=None if title is None else f"{title}, {arm} arm",
    palette=palette, file=file,
)


def main(bam, index, output_dir, jobs, max_read_length, min_map_overlap, min_candidate_overlap, max_motifs, target, min_k, max_k, max_p_adjusted, jellyfish, jellyfish_hash_size, bin_size, n_boot, palette, title, min_quality, **kwargs):
    """basic pipeline: select reads, find enriched motifs, plot"""
    get_filename = lambda fn: path.join(output_dir, fn)
    tailpuller_sam = get_filename("tailpuller.sam")
    with open(tailpuller_sam, mode="wt") as sam:
        tailpuller.main(bam, index, **get_tailpuller_kws(
            min_quality, max_read_length, min_map_overlap,
            min_candidate_overlap, target, sam,
        ))
    tailchopper_sam = get_filename("tailchopper.sam")
    with open(tailchopper_sam, mode="wt") as sam:
        tailchopper.main(tailpuller_sam, index, **get_tailchopper_kws(
            target, min_quality, sam,
        ))
    for arm, f, F in [("p", [target], ["is_q"]), ("q", [target, "is_q"], [0])]:
        repeatfinder_tsv = get_filename(f"repeatfinder-{arm}_arm.tsv")
        with open(repeatfinder_tsv, mode="wt") as tsv:
            repeatfinder.main(tailchopper_sam, **get_repeatfinder_kws(
                f, F, min_quality, max_motifs, min_k, max_k, max_p_adjusted,
                jellyfish, jellyfish_hash_size, jobs, tsv,
            ))
        kmerscanner_dat = get_filename("kmerscanner-{}_arm.dat.gz".format(arm))
        with gzopen(kmerscanner_dat, mode="wt") as dat:
            kmerscanner.main(tailpuller_sam, **get_kmerscanner_kws(
                f, F, min_quality, repeatfinder_tsv, bin_size, jobs, dat,
            ))
        densityplot_pdf = get_filename(f"densityplot-{arm}_arm.pdf")
        with open(densityplot_pdf, mode="wb") as pdf:
            try:
                densityplot.main(kmerscanner_dat, index, **get_densityplot_kws(
                    f, F, min_quality, bin_size, n_boot, palette, title,
                    arm, pdf,
                ))
            except EmptyKmerscanError:
                pass
    return 0

from os import path, walk
from numpy import inf

__doc__ = """edgeCase basic pipeline: select reads, find enriched motifs, plot

Usage: {0} basic-pipeline -x filename -o dirname [-j integer] [-m integer]
       {1}                [-n integer] [--target targetspec] [--min-k integer]
       {1}                [--max-k integer] [--max-p-adjusted float]
       {1}                [--jellyfish filename] [--jellyfish-hash-size string]
       {1}                [--window-size integer] [--palette palettespec]
       {1}                [-q integer] <bam>

Output (in --output-dir):
    * tailpuller.bam                 candidate reads
    * tailchopper.bam                telomeric regions
    * repeatfinder-p_arm.tsv         overrepresented motifs on the p arm
    * repeatfinder-q_arm.tsv         overrepresented motifs on the q arm
    * densityplot-p_arm.pdf          plot of motif densities on the p arm
    * densityplot-q_arm.pdf          plot of motif densities on the q arm

Positional arguments:
    <bam>                            name of input BAM/SAM file

Required options:
    -x, --index [filename]           location of the reference .ecx index
    -o, --output-dir [dirname]       name of the output directory (must exist)

Options:
    -j, --jobs [integer]             number of parallel jobs (for jellyfish and kmerscanner) [default: 1]
    -m, --max-read-length [integer]  maximum read length to consider when selecting lookup regions
    -n, --max-motifs [integer]       maximum number of motifs to report [default: 4]
    --target [targetspec]            an ECX flag (cut relative to reference) or 'cigar' [default: tract_anchor]
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
    lambda output_dir:
        len(next(walk(output_dir))[2]) == 0: "{} already has files in it",
    lambda max_read_length:
        max_read_length > 0: "--max-read-length below 0",
    lambda target:
        target in {"ucsc_mask_anchor", "fork", "tract_anchor", "cigar"}:
            "unknown value of --target",
    lambda min_k, max_k:
        0 < min_k < max_k: "not satisfied: 0 < m < M",
}


def main(bam, index, output_dir, jobs, max_read_length, max_motifs, target, min_k, max_k, max_p_adjusted, jellyfish, jellyfish_hash_size, window_size, palette, min_quality, **kwargs):
    return 0

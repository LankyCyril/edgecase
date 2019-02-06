from sys import stdout
from re import compile, IGNORECASE
from edgecaselib.util import ReadFileChain, MAINCHROMS
from tqdm import tqdm
from pysam import FastxFile, AlignmentFile
from pandas import read_csv, DataFrame
from itertools import takewhile, filterfalse


def get_anchors(reference):
    """Get coordinates of hard-masked bounds at each end of each main chromosome"""
    if reference.endswith(".tsv"): # assume precomputed anchors
        return read_csv(reference, sep="\t", index_col=0)
    else:
        pattern = compile(r'[^n]', flags=IGNORECASE)
        anchor_data = {}
        bar = tqdm(
            desc="Finding anchors", total=len(MAINCHROMS), unit="chromosome"
        )
        with FastxFile(reference) as genome:
            for entry in genome:
                if entry.name in MAINCHROMS:
                    bar.update()
                    bound_5prime = pattern.search(entry.sequence).span()[0]
                    bound_3prime = (
                        len(entry.sequence) -
                        pattern.search(entry.sequence[::-1]).span()[0]
                    )
                    anchor_data[entry.name] = bound_5prime, bound_3prime
        return DataFrame(data=anchor_data, index=["5prime", "3prime"]).T


def is_good_entry(entry):
    """Simple filter"""
    if entry.is_unmapped or entry.is_secondary or entry.is_supplementary:
        return False
    elif entry.reference_name in MAINCHROMS:
        return True
    else:
        return False


def filter_entries(bam_data, anchors, prime):
    """Only pass reads extending past anchors"""
    isnone = lambda p: p is None
    for entry in bam_data:
        if is_good_entry(entry):
            positions = entry.get_reference_positions(full_length=True)
            left_clip = sum(
                True for _ in takewhile(isnone, positions)
            )
            left_mappos = next(filterfalse(isnone, positions))
            right_clip = sum(
                True for _ in takewhile(isnone, reversed(positions))
            )
            right_mappos = next(filterfalse(isnone, reversed(positions)))
            if prime not in {5, 3}:
                raise ValueError("`prime` can only be 5 or 3")
            elif prime == 5:
                anchor = anchors.loc[entry.reference_name, "5prime"]
                if left_mappos - left_clip < anchor:
                    yield entry
            elif prime == 3:
                anchor = anchors.loc[entry.reference_name, "3prime"]
                if right_mappos + right_clip > anchor:
                    yield entry


def main(bams, reference, prime, file=stdout, **kwargs):
    # use header of first input file (NB! fragile):
    with AlignmentFile(bams[0]) as bam:
        print(str(bam.header).rstrip("\n"), file=file)
    # dispatch data to subroutines:
    anchors = get_anchors(reference)
    with ReadFileChain(bams, AlignmentFile) as bam_data:
        for entry in filter_entries(bam_data, anchors, prime):
            print(entry.to_string(), file=file)

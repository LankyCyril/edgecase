from sys import stdout, stderr
from pysam import AlignmentFile, FastxFile
from os import path
from re import search
from subprocess import PIPE, check_output, run
from edgecaselib.formats import filter_bam
from edgecaselib.formats import MEME_HEADER_FMT, MEME_REGEX_HEADER_FMT
from edgecaselib.util import get_executable, circularize_motif
from pandas import DataFrame, concat
from uuid import uuid4
from regex import finditer, IGNORECASE
from itertools import count
from tempfile import TemporaryDirectory


def guess_bg_fmt(background):
    """Decide if `background` is in SAM/BAM format or is a MEME HMM"""
    if background is None:
        return None
    with open(background, mode="rb") as bg_handle:
        line = next(bg_handle)
    if search(br'^#.*Markov frequencies', line):
        return "hmm"
    else:
        return "sam"


def interpret_args(fmt, fasta_get_markov, meme, background):
    """Parse and check arguments"""
    if fmt == "sam":
        manager = AlignmentFile
    elif fmt == "fastx":
        manager = FastxFile
    else:
        raise ValueError("Unsupported --fmt: '{}'".format(fmt))
    bg_fmt = guess_bg_fmt(background)
    if bg_fmt == "sam":
        error_mask = "No {} found, needed to generate background"
        if get_executable("fasta-get-markov", fasta_get_markov, False) is None:
            raise ValueError(error_mask.format("fasta-get-markov"))
    return (
        manager,
        get_executable("fasta-get-markov", fasta_get_markov, False),
        get_executable("meme", meme),
        bg_fmt
    )


def interpret_lengths(lengths):
    """Convert a lengths string into an array of minw and maxw"""
    lengths_array = []
    for span in lengths.split(","):
        fields = span.split("-")
        if len(fields) > 2:
            raise ValueError("Incorrect syntax: '{}'".format(lengths))
        else:
            lengths_array.append([int(fields[0]), int(fields[-1])])
    return lengths_array


def convert_background(sam, tempdir, fasta_get_markov, max_order=6, max_entries=1000000):
    """Convert a SAM/BAM file into Markov background for MEME"""
    print("SAM/BAM -> HMM...", file=stderr, flush=True)
    fasta_lines = []
    i = 0
    with AlignmentFile(sam) as alignment:
        for entry in alignment:
            if i >= max_entries:
                break
            elif entry.seq:
                fasta_lines.extend([">"+entry.qname, entry.seq])
                i += 1
    hmm = run(
        [fasta_get_markov, "-m", str(max_order)],
        input="\n".join(fasta_lines), encoding="ascii", stdout=PIPE
    )
    bfile = path.join(tempdir, "bfile")
    with open(bfile, mode="wt") as bfile_handle:
        print(hmm.stdout, file=bfile_handle)
    print("...done", file=stderr, flush=True)
    return bfile


def convert_input(bam, manager, tempdir, samfilters):
    """Convert BAM to fasta for MEME"""
    fasta = path.join(tempdir, "input.fa")
    with manager(bam) as alignment, open(fasta, mode="wt") as fasta_handle:
        for entry in filter_bam(alignment, samfilters, "SAM/BAM -> FASTA"):
            print(
                ">{}\n{}".format(entry.qname, entry.query_sequence),
                file=fasta_handle
            )
    return fasta


def alpha_inversion(motif):
    """Get alphanumerically smallest inversion of motif"""
    return min(motif[i:]+motif[:i] for i in range(len(motif)))


def parse_meme_output(meme_txt, expect_nmotifs):
    """Parse brief MEME output"""
    with open(meme_txt) as meme_report_handle:
        meme_report = meme_report_handle.readlines()
    meme_df_as_list = []
    for line in meme_report:
        header_matcher = search(MEME_HEADER_FMT, line)
        if header_matcher:
            meme_df_as_list.append([
                alpha_inversion(header_matcher.group(1)),
                int(header_matcher.group(2)), float(header_matcher.group(3))
            ])
    if len(meme_df_as_list) == 0:
        return None, True
    stop_flag = (len(meme_df_as_list) < expect_nmotifs)
    meme_df = DataFrame(data=meme_df_as_list, columns=["motif", "meme_id", "e"])
    meme_df["regex"] = None
    for i, line in enumerate(meme_report):
        regex_header_matcher = search(MEME_REGEX_HEADER_FMT, line)
        if regex_header_matcher:
            _id = int(regex_header_matcher.group(1))
            _regex = meme_report[i+2].strip()
            meme_df.loc[meme_df["meme_id"]==_id, "regex"] = _regex
    meme_df = meme_df.drop(columns="meme_id")
    return meme_df.groupby(["motif", "regex"], as_index=False).min(), stop_flag


def find_and_mask_motifs(current_readfile, tempdir, meme_output_dir, expect_nmotifs):
    """Parse the MEME output"""
    meme_report, stop_flag = parse_meme_output(
        path.join(meme_output_dir, "meme.txt"), expect_nmotifs
    )
    if meme_report is None:
        return None, current_readfile, True
    circular_motifs = {
        circularize_motif(motif) for motif in meme_report["regex"]
    }
    maskable_regex = r'|'.join(circular_motifs).lower()
    next_readfile = path.join(tempdir, str(uuid4()) + ".fa")
    with FastxFile(current_readfile) as unmasked:
        with open(next_readfile, mode="wt") as masked:
            for entry in unmasked:
                seq_as_list = list(entry.sequence)
                match_iter = finditer(
                    maskable_regex, entry.sequence,
                    overlapped=True, flags=IGNORECASE
                )
                for match in match_iter:
                    s, e = match.span()
                    seq_as_list[s:e] = ["N"] * (e - s)
                print(">" + entry.name, file=masked)
                print("".join(seq_as_list), file=masked)
    return meme_report, next_readfile, stop_flag


def run_meme(meme, jobs, readfile, background, lengths_array, nmotifs, evt, tempdir):
    """Run the MEME binary with preset parameters"""
    current_readfile = readfile
    meme_reports = []
    for minw, maxw in lengths_array:
        for it in count():
            meme_output_dir = path.join(
                tempdir, str(minw)+"-"+str(maxw)+"."+str(it)
            )
            if background is None:
                check_output([
                    meme, "-p", str(jobs), "-mod", "anr",
                    "-dna", "-minsites", "2", "-nmotifs", str(nmotifs),
                    "-minw", str(minw), "-maxw", str(maxw), "-evt", str(evt),
                    "-brief", "0", "-oc", meme_output_dir, current_readfile
                ])
            else:
                check_output([
                    meme, "-p", str(jobs), "-bfile", background, "-mod", "anr",
                    "-dna", "-minsites", "2", "-nmotifs", str(nmotifs),
                    "-minw", str(minw), "-maxw", str(maxw), "-evt", str(evt),
                    "-brief", "0", "-oc", meme_output_dir, current_readfile
                ])
            next_meme_report, current_readfile, stop_flag = find_and_mask_motifs(
                current_readfile, tempdir, meme_output_dir, nmotifs
            )
            if next_meme_report is not None:
                meme_reports.append(next_meme_report)
            if stop_flag:
                break
    return concat(meme_reports, axis=0)


def main(readfile, fmt, flags, flags_any, flag_filter, min_quality, fasta_get_markov, background, meme, lengths, nmotifs, evt, jobs=1, file=stdout, **kwargs):
    # parse arguments
    manager, fasta_get_markov, meme, bg_fmt = interpret_args(
        fmt, fasta_get_markov, meme, background
    )
    lengths_array = interpret_lengths(lengths)
    with TemporaryDirectory() as tempdir:
        if bg_fmt == "sam": # will need to convert SAM to HMM
            background = convert_background(
                background, tempdir, fasta_get_markov
            )
        if manager != FastxFile: # will need to convert SAM to fastx
            samfilters = [flags, flags_any, flag_filter, min_quality]
            readfile = convert_input(readfile, manager, tempdir, samfilters)
        combined_meme_report = run_meme(
            meme, jobs, readfile, background, lengths_array,
            nmotifs, evt, tempdir
        )
    combined_meme_report.to_csv(file, sep="\t", index=False)

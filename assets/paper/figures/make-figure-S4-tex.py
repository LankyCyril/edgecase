#!/usr/bin/env python
from sys import argv
from os import path
from tqdm import tqdm
from subprocess import check_output
from re import search
from collections import OrderedDict, defaultdict


SOURCE_DIR = "./Figure_S4"
SUBJECTS = [f"HG00{i}" for i in range(1, 8)]

CHROMS = [
    "chr2",
    "3ptel_1-500K_1_12_12",
    "4ptel_1-500K_1_12_12",
    "chr5",
    "chr9",
    "chr12",
    "17ptel_1_500K_1_12_12",
]

pad = lambda s, n: "\hspace*{"+str(n*3.9)+"pt}"+s+"\hspace*{"+str(n*3.9)+"pt}"
ADJUSTED_NAMES = {
    "chr2": pad("2p (chr2)", 2),
    "3ptel_1-500K_1_12_12": pad("3p (3ptel_1...)", 2),
    "4ptel_1-500K_1_12_12": pad("4p (4ptel_1...)", 30),
    "chr5": pad("5p (chr5)", 1),
    "chr9": pad("9p (chr9)", 0),
    "chr12": pad("12p (chr12)", 0),
    "17ptel_1_500K_1_12_12": pad("17p (17ptel_1...)", 13),
}

TEX_LEGEND = r'''\begin{textblock}{13}($X,$Y)
\includegraphics[width=$WIDTHin,keepaspectratio]{Figure_4/legend.pdf}
\end{textblock}'''
LEGEND_WIDTH = 4.1
EXTRA_MARGIN = 3
Y_MINOR_SKIP = .020
Y_MAJOR_SKIP = .08

TEX_HEADER = r'''\documentclass{article}
\usepackage[paperheight=$HEIGHTin,paperwidth=$WIDTHin,margin=0in]{geometry}
\usepackage[sfdefault]{roboto}
\usepackage{graphicx}
\usepackage{tikz}
\usepackage{ulem}
    \renewcommand{\ULdepth}{7pt}
\usepackage[absolute,overlay]{textpos}
    \setlength{\TPHorizModule}{1in}
    \setlength{\TPVertModule}{1in}
\makeatletter
    \newcommand*{\textoverline}[1]{$\overline{\hbox{#1\vphantom{\"A}}}\m@th$}
    \makeatother
\begin{document}'''

TEX_IMAGE = r'\begin{textblock}{13}($X,$Y)\includegraphics{$PDF}\end{textblock}'
TEX_YLABEL = r'\begin{textblock}{13}($X,$Y)\rotatebox{90}{\Large{\textoverline{$T}}}\end{textblock}'

TEX_FOOTER = r'\end{document}'


def get_pdfs(source_dir, chroms, subjects):
    for chrom in chroms:
        for subject in subjects:
            pdf = f"{source_dir}/{chrom}-{subject}.pdf"
            if path.isfile(pdf):
                yield pdf, (chrom, subject)


def get_pdf_sizes(pdfs):
    for pdf in tqdm(pdfs, desc="Identifying sizes"):
        imagemagick_bytes = check_output(["identify", "-verbose", pdf])
        for line in imagemagick_bytes.decode().split("\n"):
            if line.strip().startswith("Print size:"):
                matcher = search(r'([0-9.]+)x([0-9.]+)', line)
                if matcher:
                    w, h = float(matcher.group(1)), float(matcher.group(2))
                    yield pdf, (w, h)
                    break
        else:
            yield pdf, (None, None)


def get_ylabels(pdfs, pdf_sizes, combined_width):
    chrom2xs, chrom2ys = defaultdict(list), defaultdict(list)
    y = .05
    for pdf, (chrom, subject) in pdfs.items():
        chrom2ys[chrom].append(y if (y != .05) else .47)
        width, height = pdf_sizes[pdf]
        if "chr2-HG001" not in pdf:
            chrom2xs[chrom].append(width+len(pdf)/22)
        y += height + (
            Y_MAJOR_SKIP if (
                (subject == "HG007") or
                ((chrom == "chr5") and (subject == "HG006")) or
                ((chrom == "chr9") and (subject == "HG005"))
            )
            else Y_MINOR_SKIP
        )
    for chrom, xs in chrom2xs.items():
        x = min(xs) - .3
        y = min(chrom2ys[chrom])
        yield x, y, ADJUSTED_NAMES[chrom].replace("_", r'\_')


def main(argv):
    pdfs = OrderedDict(get_pdfs(SOURCE_DIR, CHROMS, SUBJECTS))
    pdf_sizes = OrderedDict(get_pdf_sizes(pdfs))
    combined_width = max(width for width, _ in pdf_sizes.values())
    combined_height = (
        sum(height+Y_MINOR_SKIP for _, height in pdf_sizes.values()) +
        (Y_MAJOR_SKIP-Y_MINOR_SKIP) * len(CHROMS)
    )
    print((TEX_HEADER
        .replace("$WIDTH", format(combined_width+EXTRA_MARGIN, ".3f"))
        .replace("$HEIGHT", format(combined_height, ".3f"))
    ))
    y = .05
    for pdf, (chrom, subject) in pdfs.items():
        width, height = pdf_sizes[pdf]
        print((TEX_IMAGE
            .replace("$X", str(
                (-.11 if (subject=="HG007") else 0)
            ))
            .replace("$Y", format(y, ".3f")).replace("$PDF", pdf)
        ))
        y += height + (
            Y_MAJOR_SKIP if (
                (subject == "HG007") or
                ((chrom == "chr5") and (subject == "HG006")) or
                ((chrom == "chr9") and (subject == "HG005"))
            )
            else Y_MINOR_SKIP
        )
    for x, y, ylabel in get_ylabels(pdfs, pdf_sizes, combined_width):
        print((TEX_YLABEL
            .replace("$X", str(x)).replace("$Y", str(y)).replace("$T", ylabel)
        ))
    print((TEX_LEGEND
        .replace("$WIDTH", format(LEGEND_WIDTH, ".3f"))
        .replace("$X", format(
            combined_width+EXTRA_MARGIN-LEGEND_WIDTH*1.05, ".3f",
        ))
        .replace("$Y", ".4")
    ))
    print(TEX_FOOTER)
    return 0


if __name__ == "__main__":
    returncode = main(argv)
    exit(returncode)

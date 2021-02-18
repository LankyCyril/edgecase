#!/usr/bin/env python
from sys import argv
from tqdm import tqdm
from subprocess import check_output
from re import search
from collections import OrderedDict


PDFS = [
    "Figure_5/chr7.pdf",
    "Figure_5/chr8.pdf",
    "Figure_5/chr11.pdf",
    "Figure_5/chr12.pdf",
    "Figure_5/14qtel_1-500K_1_12_12_rc.pdf",
    "Figure_5/chr15.pdf",
    "Figure_5/18qtel_1-500K_1_12_12_rc.pdf",
]


TEX_LEGEND = r'''\begin{textblock}{13}($X,$Y)
\includegraphics[width=$WIDTHin,keepaspectratio]{Figure_5/legend.pdf}
\end{textblock}'''
LEGEND_WIDTH = 3.6
EXTRA_MARGIN = 1.5
LEFT_SHIFT = .2
UPWARD_SHIFT = .12


TEX_HEADER = r'''\documentclass{article}
\usepackage[paperheight=$HEIGHTin,paperwidth=$WIDTHin,margin=0in]{geometry}
\usepackage[sfdefault]{roboto}
\usepackage{graphicx}
\usepackage{tikz}
\usepackage[absolute,overlay]{textpos}
    \setlength{\TPHorizModule}{1in}
    \setlength{\TPVertModule}{1in}
\begin{document}'''

TEX_IMAGE = r'\begin{textblock}{13}($X,$Y)\includegraphics{$PDF}\end{textblock}'

TEX_FOOTER = r'\end{document}'


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


def main(argv):
    pdf_sizes = OrderedDict(get_pdf_sizes(PDFS))
    combined_width = max(width for width, _ in pdf_sizes.values())
    combined_height = (
        sum(height for _, height in pdf_sizes.values()) -
        UPWARD_SHIFT * (len(pdf_sizes) - 1)
    )
    header = (TEX_HEADER
        .replace("$WIDTH", format(combined_width+EXTRA_MARGIN, ".3f"))
        .replace("$HEIGHT", format(combined_height, ".3f"))
    )
    print(header)
    y = 0
    for pdf, (width, height) in pdf_sizes.items():
        section = (TEX_IMAGE
            .replace("$X", format(combined_width-width-LEFT_SHIFT, ".3f"))
            .replace("$Y", format(y, ".3f"))
            .replace("$PDF", pdf)
        )
        print(section)
        y += height - UPWARD_SHIFT
    legend_section = (TEX_LEGEND
        .replace("$WIDTH", format(LEGEND_WIDTH, ".3f"))
        .replace("$X", format(
            combined_width+EXTRA_MARGIN-LEGEND_WIDTH*1.05, ".3f",
        ))
        .replace("$Y", "0")
    )
    print(legend_section)
    print(TEX_FOOTER)
    return 0


if __name__ == "__main__":
    returncode = main(argv)
    exit(returncode)

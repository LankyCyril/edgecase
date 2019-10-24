#!/usr/bin/env python
from sys import argv
from re import search
from io import StringIO
from pandas import read_fwf, concat
from pysam import FastxFile

SECTION_SEPARATOR = "*"*80
SUBSECTION_SEPARATOR = "-"*80

def parse_meme(*, filename=None, lines=None, separator=SECTION_SEPARATOR):
    if lines is None:
        if filename is None:
            raise ValueError("Nothing to parse")
        else:
            with open(filename, mode="rt") as memehandle:
                lines = list(map(str.strip, memehandle))
    elif filename is not None:
        raise ValueError("Both `filename` and `lines` supplied")
    separator_positions = {
        i for i, line in enumerate(lines) if line == separator
    }
    empty_positions = {
        i for i, line in enumerate(lines) if line == ""
    }
    header_positions = sorted({
        i for i in separator_positions
        if (i+2 in separator_positions) and (i+1 not in empty_positions)
    })
    return {
        lines[head_pre+1]: [
            line for line in lines[head_pre+3:next_head_pre]
            if line not in {"", separator}
        ]
        for head_pre, next_head_pre
        in zip(header_positions, header_positions[1:])
    }

def parse_meme_nested(memetxt):
    return {
        section_name: parse_meme(
            lines=section_lines, separator=SUBSECTION_SEPARATOR
        )
        for section_name, section_lines
        in parse_meme(filename=memetxt, separator=SECTION_SEPARATOR).items()
    }

def parse_sections(meme_sections):
    meme_position_data_list = []
    for section_name, section_data in meme_sections.items():
        if search(r'^MOTIF\s*[ACGT]+.+MEME', section_name):
            motif = search(r'^MOTIF\s*([ACGT]+)', section_name).group(1)
            for subsection_name, subsection_data in section_data.items():
                if search(r'sites sorted by position', subsection_name):
                    clean_lines = StringIO("\n".join((
                        line for line in subsection_data
                        if search(motif, line)
                    )))
                    motif_data = read_fwf(
                        clean_lines, usecols=(0, 1, 2, 4),
                        names=["prefix", "pos", "pval", "motif"]
                    )
                    meme_position_data_list.append(motif_data)
    return concat(meme_position_data_list, axis=0, ignore_index=True)

def generate_prefix_to_name(meme_position_data, fasta_records):
    prefixes = set(meme_position_data["prefix"])
    names = set(fasta_records.keys())
    prefix_to_name = {}
    for prefix in prefixes:
        for name in names:
            if name.startswith(prefix):
                prefix_to_name[prefix] = name
                names.remove(name)
                break
    for prefix in prefixes:
        for name in names:
            if name.startswith(prefix):
                raise KeyError("Ambiguous prefix: '{}'".format(prefix))
    return prefix_to_name

def mask_records(fasta_records, meme_position_data, prefix_to_name, aggressive=False):
    masked_records = {}
    for prefix, name in prefix_to_name.items():
        position_data = meme_position_data.loc[
            meme_position_data["prefix"]==prefix, ["pos", "motif"]
        ]
        if aggressive:
            sequence = fasta_records[name]
            for _, (_, motif) in position_data.iterrows():
                sequence = sequence.replace(motif, "N"*len(motif))
        else:
            sequence_as_list = list(fasta_records[name])
            for _, (pos, motif) in position_data.iterrows():
                sequence_as_list[pos:len(motif)] = ["N"]*len(motif)
            sequence = "".join(sequence_as_list)
        masked_records[name] = sequence
    return masked_records

def main(*, fasta, memetxt, aggressive=False):
    fasta_records = {entry.name: entry.sequence for entry in FastxFile(fasta)}
    meme_sections = parse_meme_nested(memetxt)
    meme_position_data = parse_sections(meme_sections)
    prefix_to_name = generate_prefix_to_name(meme_position_data, fasta_records)
    masked_records = mask_records(
        fasta_records, meme_position_data, prefix_to_name, aggressive=aggressive
    )
    for name, sequence in masked_records.items():
        print(">{}\n{}".format(name, sequence))
    return 0

if __name__ == "__main__":
    aggressive = (len(argv) > 3) and (argv[3] == "--aggressive")
    returncode = main(fasta=argv[1], memetxt=argv[2], aggressive=aggressive)
    exit(returncode)

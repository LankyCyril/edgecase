#!/usr/bin/env python
from sys import argv
from pandas import read_csv
from subprocess import Popen, PIPE


def get_header_height(filename):
    with open(filename, mode="rt") as handle:
        for i, line in enumerate(handle):
            if line[0] != "#":
                break
    return i


def RunCounter(iterable):
    previous_value = None
    current_count = 1
    for value in iterable:
        if previous_value is None:
            previous_value = value
        elif value == previous_value:
            current_count += 1
        else:
            yield (previous_value, current_count)
            current_count = 1
            previous_value = value
    yield (previous_value, current_count)


def init_table(ncols):
    return (
        r'\begin{samepage} \begin{table}[h!] \small \begin{tabular}{' +
        "l"*ncols + "}\n\\hline"
    )


def preformat_header(tsv, level):
    name_counts = RunCounter(tsv.columns.get_level_values(level))
    for name, count in name_counts:
        realname = (
            "" if name.startswith("Unnamed:") else name.replace("_", r'\_')
        )
        if count == 1:
            yield r'\textbf{'+realname+"}"
        else:
            yield r'\multicolumn{'+str(count)+r'}{l}{\textbf{'+realname+"}}"


def preformat_row(row):
    for label, value in row.items():
        is_p = False
        if isinstance(label, (tuple, list)):
            if len(set(label) & {"p", "p_adjusted"}):
                is_p = True
        elif label in {"p", "p_adjusted"}:
            is_p = True
        if is_p:
            try:
                float_value = float(value)
            except:
                yield value.replace("_", r'\_')
            else:
                if float(value) == 0.0:
                    yield r'<1e-300'
                elif float(value) == 1.0:
                    yield "1.00"
                else:
                    yield format(float(value), ".2e")
        else:
            try:
                float_value = float(value)
            except ValueError:
                yield value.replace("_", r'\_')
            else:
                yield format(float_value, ".6f")


def columnize(list_of_iterables):
    raw_rows = []
    for iterable in list_of_iterables:
        if isinstance(iterable, str) and (iterable == r'\hline'):
            raw_rows.append(iterable)
        else:
            raw_rows.append("\t&\t".join(iterable)+" \\\\\n")
    column_t = Popen(
        ["column", "-t", "-s\t"], stdin=PIPE, stdout=PIPE,
        universal_newlines=True,
    )
    return column_t.communicate(input="\n".join(raw_rows))[0].rstrip("\n")


def finish_table():
    return "\n".join([
        r'\hline', r'\end{tabular}', r'\caption{}',
        r'\label{}', r'\end{table}', r'\end{samepage}',
    ])


def main(filename):
    header_height = get_header_height(filename)
    tsv = read_csv(
        filename, sep="\t", escapechar="#",
        header=list(range(header_height)),
        dtype=str,
    )
    print(init_table(ncols=tsv.shape[1]))
    print(columnize((
        [preformat_header(tsv, n) for n in range(header_height)] +
        [r'\hline'] +
        [preformat_row(row) for _, row in tsv.iterrows()]
    )))
    print(finish_table())
    return 0


if __name__ == "__main__":
    returncode = main(filename=argv[1])
    exit(returncode)

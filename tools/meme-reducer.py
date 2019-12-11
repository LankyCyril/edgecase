#!/usr/bin/env python
from sys import argv
from json import loads, dumps
from re import search


def parse_html(memehtml):
    html_pre_lines, json_lines, html_post_lines = [], [], []
    is_inside_pre, is_inside_json, is_inside_post = True, False, False
    with open(memehtml, mode="rt") as html:
        for line in map(str.strip, html):
            if line.startswith("<script"):
                if is_inside_post:
                    html_post_lines.append(line)
                else:
                    is_inside_pre, is_inside_json = False, True
            elif line.startswith("</script"):
                if is_inside_post:
                    html_post_lines.append(line)
                else:
                    is_inside_pre, is_inside_json, is_inside_post = (
                        False, False, True
                    )
            elif is_inside_json and (not line.startswith("//")):
                if line == "var data = {":
                    json_lines.append("{")
                elif line == "};":
                    json_lines.append("}")
                else:
                    json_lines.append(line)
            elif is_inside_pre:
                html_pre_lines.append(line)
            elif is_inside_post:
                html_post_lines.append(line)
    return (
        "\n".join(html_pre_lines),
        loads(" ".join(json_lines)),
        "\n".join(html_post_lines)
    )


def isin_motifs(known_motifs, motif):
    if motif in known_motifs:
        return motif
    else:
        for known_motif in known_motifs:
            if len(motif) == len(known_motif):
                if search(motif, known_motif*2):
                    return known_motif
        else:
            return False


def get_unique_motifs(motifs_json):
    known_motifs = set()
    for motif_json in motifs_json:
        if not isin_motifs(known_motifs, motif_json["id"]):
            known_motifs.add(motif_json["id"])
    return sorted(known_motifs, key=lambda x:len(x), reverse=True)


def get_longest_unique_motifs(unique_motifs):
    longest_unique_motifs = set()
    unique_motifs_set = set(unique_motifs)
    for motif in unique_motifs:
        for known_motif in unique_motifs:
            if known_motif in unique_motifs_set:
                if len(motif) >= len(known_motif):
                    if search(known_motif, motif):
                        longest_unique_motifs.add(motif)
                        unique_motifs_set.remove(known_motif)
    return sorted(longest_unique_motifs, key=lambda x:len(x), reverse=True)


def filter_json(meme_json, unique_motifs):
    filtered_json = {}
    for key, value in meme_json.items():
        if key != "motifs":
            filtered_json[key] = value
        else:
            motifs_json, filtered_motifs_json = value, []
            for motif_json in motifs_json:
                known_motif = isin_motifs(unique_motifs, motif_json["id"])
                if known_motif:
                    filtered_motifs_json.append(motif_json)
                    unique_motifs.remove(known_motif)
            filtered_json["motifs"] = filtered_motifs_json
    return filtered_json


def main(memehtml):
    html_pre, meme_json, html_post = parse_html(memehtml)
    unique_motifs = get_unique_motifs(meme_json["motifs"])
    unique_motifs = get_longest_unique_motifs(unique_motifs)
    filtered_json = filter_json(meme_json, unique_motifs)
    print(html_pre)
    print("<script>\nvar data = ", end="")
    print(dumps(filtered_json, indent=4), end="")
    print(";\n</script>")
    print(html_post)
    return 0


if __name__ == "__main__":
    returncode = main(argv[1])
    exit(returncode)

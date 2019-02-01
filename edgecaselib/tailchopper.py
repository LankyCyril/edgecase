from types import SimpleNamespace

def chop(entry, anchors, target):
    """Return only part of sequence extending past anchor"""
    if target[0] == "5":
        chop_pos = anchors.loc[entry.reference_name, "5prime"]
        return SimpleNamespace(
            name=entry.query_name,
            sequence=entry.query_sequence[:chop_pos]
        )
    if target[0] == "3":
        chop_pos = anchors.loc[entry.reference_name, "3prime"]
        return SimpleNamespace(
            name=entry.query_name,
            sequence=entry.query_sequence[chop_pos:]
        )

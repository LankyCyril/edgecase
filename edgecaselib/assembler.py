from sys import stdout, stderr
from edgecaselib.util import natsorted_chromosomes
from edgecaselib.formats import filter_bam, load_index
from pysam import AlignmentFile
from collections import OrderedDict, defaultdict
from tqdm import tqdm
from pandas import Series
from numpy import percentile
from networkx import DiGraph, compose
from functools import reduce


def collect_sequence_unimers(sequence, k, desc=None):
    """Collect identities and positions of unique k-mers in a sequence"""
    pos2unimer, unimer2pos = OrderedDict(), {}
    repeated_unimers = set()
    if desc:
        pos_iterator = tqdm(range(len(sequence)-k), desc=desc)
    else:
        pos_iterator = range(len(sequence)-k)
    for pos in pos_iterator:
        unimer = sequence[pos:pos+k].upper()
        if unimer not in repeated_unimers:
            if unimer in unimer2pos:
                repeated_unimers.add(unimer)
                del pos2unimer[unimer2pos[unimer]]
                del unimer2pos[unimer]
            else:
                unimer2pos[unimer] = pos
                pos2unimer[pos] = unimer
    return pos2unimer


def build_unimer_database(bam, samfilters, chromosome, unimer_size):
    """Collect identities and positions of unique k-mers in all reads"""
    with AlignmentFile(bam) as alignment:
        decorated_entry_iterator = filter_bam(
            alignment.fetch(chromosome), samfilters,
            desc="Collecting unimers for reads mapped to {}".format(chromosome)
        )
        return {
            entry.qname: collect_sequence_unimers(entry.seq, unimer_size)
            for entry in decorated_entry_iterator
        }


def filter_unimer_database(unimer_database, cutoff_percentile=95):
    """Remove infrequent unimers from database"""
    decorated_count_iterator = tqdm(
        unimer_database.items(), desc="Counting unimers", unit="read"
    )
    unimer_counts = defaultdict(int)
    for name, p2u in decorated_count_iterator:
        for unimer in p2u.values():
            unimer_counts[unimer] += 1
    unimer_count_table = Series(unimer_counts)
    print("Total unimers:", len(unimer_count_table), file=stderr)
    count_filter = (
        unimer_count_table >= percentile(unimer_count_table, cutoff_percentile)
    )
    frequent_unimers = set(unimer_count_table[count_filter].index)
    print("Retained frequent unimers:", len(frequent_unimers), file=stderr)
    decorated_filter_iterator = tqdm(
        unimer_database.items(), desc="Filtering unimers", unit="read"
    )
    return {
        name: {p: u for p, u in p2u.items() if u in frequent_unimers}
        for name, p2u in decorated_filter_iterator
    }


def assemble_related_reads(bam, chromosome, unimer_size, samfilters, output_prefix):
    """Perform local assembly on one chromosome arm"""
    raw_unimer_database = build_unimer_database(
        bam, samfilters, chromosome, unimer_size
    )
    unimer_database = filter_unimer_database(raw_unimer_database)
    decorated_readgraph_iterator = tqdm(
        unimer_database.items(), desc="Building read paths", unit="read"
    )
    paths = []
    for name, p2u in decorated_readgraph_iterator:
        path, nodes = DiGraph(), list(p2u.values())
        path.add_edges_from(zip(nodes, nodes[1:]))
        paths.append(path)
    raw_assembly_graph = reduce(
        compose, tqdm(paths, desc="Gluing read paths", unit="read")
    )
    print("Nodes in raw assembly graph:", len(raw_assembly_graph.nodes()))
    print("Edges in raw assembly graph:", len(raw_assembly_graph.edges()))


def main(bam, index, flags, flags_any, flag_filter, min_quality, unimer_size, chromosomes, output_prefix, jobs=1, file=stdout, **kwargs):
    # parse and check arguments:
    if chromosomes:
        chromosome_subset = set(chromosomes.split("|"))
    else:
        chromosome_subset = set()
    if index:
        indexed_chromosomes = set(load_index(index)["rname"])
    else:
        indexed_chromosomes = set()
    target_chromosomes = natsorted_chromosomes(
        chromosome_subset & indexed_chromosomes
    )
    if not target_chromosomes:
        raise ValueError("No usable chromosomes found")
    else:
        for chromosome in target_chromosomes:
            assemble_related_reads(
                bam, chromosome, unimer_size,
                [flags, flags_any, flag_filter, min_quality],
                output_prefix
            )

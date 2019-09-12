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
from itertools import count


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


def compose_raw_assembly_graph(unimer_database):
    """String reads through retained unimers and glue paths together"""
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
    msg = "Raw assembly graph: N={}, E={}".format(
        len(raw_assembly_graph.nodes()), len(raw_assembly_graph.edges())
    )
    print(msg, file=stderr)
    return raw_assembly_graph


def reduce_transitive_edges(G):
    """Reduce transitive edges"""
    for _ in tqdm(count(), desc="Reducing transitive edges", unit="pass"):
        node_pool = list(G.nodes())
        for node in node_pool:
            if (G.in_degree(node) == 1) and (G.out_degree(node) == 1):
                parent = next(G.predecessors(node))
                child = next(G.successors(node))
                G.remove_edge(parent, node)
                G.remove_edge(node, child)
                assert G.degree(node) == 0
                G.remove_node(node)
                G.add_edge(parent, child)
        if len(G.nodes()) == len(node_pool):
            return G


def simplify_assembly_graph(assembly_graph):
    """Reduce transitive edges and maybe something more"""
    reduced_assembly_graph = reduce_transitive_edges(assembly_graph)
    msg = "Simplified assembly graph: N={}, E={}".format(
        len(reduced_assembly_graph.nodes()), len(reduced_assembly_graph.edges())
    )
    print(msg, file=stderr)
    return reduced_assembly_graph


def assemble_related_reads(bam, chromosome, unimer_size, samfilters, output_prefix):
    """Perform local assembly on one chromosome arm"""
    raw_unimer_database = build_unimer_database(
        bam, samfilters, chromosome, unimer_size
    )
    unimer_database = filter_unimer_database(raw_unimer_database)
    raw_assembly_graph = compose_raw_assembly_graph(unimer_database)
    assembly_graph = simplify_assembly_graph(raw_assembly_graph)


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

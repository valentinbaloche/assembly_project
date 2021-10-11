#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import random
import sys
from operator import itemgetter

from random import randint
random.seed(9001)
import matplotlib as plt
import networkx as nx


__authors__ = "Valentin Baloche & Alix de Thoisy"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Valentin Baloche & Alix de Thoisy"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Valentin Baloche & Alix de Thoisy"
__email__ = "valentin.baloche@gmail.fr, alixdet@protonmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open (fastq_file, 'r') as ff:
        lines = ff.readlines()
        for i in range(len(lines)) : # try using next
            if i%4 == 1:
                yield lines[i][:-1]


def cut_kmer(read, kmer_size):
    for i in range(len(read) - kmer_size + 1):
        yield read[i:i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            kmer_dict.setdefault(kmer, 0)
            kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        graph.add_edge(kmer[:-1], kmer[1:], weight=kmer_dict[kmer])
    return graph


def get_starting_nodes(graph):
    starting_nodes = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    sink_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, starting_nodes, sink_nodes):
    """
    prend un graphe, une liste de noeuds d'entrée et une liste de noeuds de sortie
    et retourne une liste de tuples (contig, longueur du contig)
    """
    contigs = []
    for starting_node in starting_nodes:
        for sink_node in sink_nodes:
            if nx.has_path(graph, starting_node, sink_node):
                for path in nx.all_simple_paths(graph, starting_node, sink_node):
                    contig = path[0]
                    for i in range(1, len(path)):
                        contig += path[i][-1]
                    contigs.append((contig, len(contig)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    with open(output_file, 'w+') as o_file:
        i = 0
        for contig in contigs_list:
            o_file.write('>contig_%i len=%i\n%s\n' % (i,
                                                      contig[1],
                                                      fill(contig[0])))
            i += 1

"""
3. Simplification du graphe de De Bruijn
"""
def path_average_weight(graph, path):
    total_weight = 0
    i = 0
    for edge in graph.subgraph(path).edges(data=True):
        total_weight += edge[2]['weight']
        i += 1
    return total_weight/i


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node:
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph


def select_best_path(graph,
                     path_list,
                     path_length,
                     weight_avg_list,
                     delete_entry_node=False,
                     delete_sink_node=False):
    if std(weight_avg_list) > 0:
        del path_list[weight_avg_list.index(max(weight_avg_list))]
        remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    elif std(weight_avg_list) == 0:
        if std(path_length) > 0:
            del path_list[path_length.index(max(path_length))]
            remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

        elif std(path_length) == 0:
            del path_list[randint(0, len(path_list)-1)]
            remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph


def std(data):
    total = 0
    for d in data:
        total += d
    mean = total/len(data)
    squared_dev = 0
    for d in data:
        squared_dev += (d - mean)**2
        print (d)
    return (squared_dev/(len(data) - 1))**0.5


def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = []
    path_length = []
    weight_avg_list = []

    for path in nx.all_simple_paths(graph, ancestor_node, descendant_node):
        path_list.append(path)
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))

    select_best_path(graph,
                     path_list,
                     path_length,
                     weight_avg_list,
                     delete_entry_node = False,
                     delete_sink_node = False)

    return graph


def simplify_bubbles(graph):
    potential_ancestors = []
    potential_descendants = []

    #looking for potential ancestors and descendants:
    for node in graph.nodes:
        if len(list(graph.successors(node))) > 1:
            potential_ancestors.append(node)

        if len(list(graph.predecessors(node))) > 1:
            potential_descendants.append(node)

    #looking for bubbles:
    for pot_anc in potential_ancestors:
        for pot_desc in potential_descendants:
            if len(list(nx.all_simple_paths(graph, pot_anc, pot_desc))) > 1:
                solve_bubble(graph, pot_anc, pot_desc)

    return graph


def solve_entry_tips(graph, starting_nodes):
    tip_nodes = []
    path_list = []
    path_length = []
    weight_avg_list = []

    for node in graph.nodes:
        if len(list(graph.predecessors(node))) > 1:
            tip_nodes.append(node)

    i = 0
    for node in tip_nodes:
        for start_node in starting_nodes:
            if nx.has_path(graph, start_node, node):
                for path in nx.all_simple_paths(graph, start_node, node):
                    path_list.append(path)
                path_length.append(len(path_list[i]))
                weight_avg_list.append(path_average_weight(graph, path_list[i]))
                i += 1

        select_best_path(graph,
                         path_list,
                         path_length,
                         weight_avg_list,
                         delete_entry_node=True,
                         delete_sink_node=False)

    return graph

def solve_out_tips(graph, sink_nodes):
    tip_nodes = []
    path_list = []
    path_length = []
    weight_avg_list = []

    for node in graph.nodes:
        if len(list(graph.successors(node))) > 1:
            tip_nodes.append(node)

    i = 0
    for node in tip_nodes:
        for sink_node in sink_nodes:
            if nx.has_path(graph, node, sink_node):
                for path in nx.all_simple_paths(graph, node, sink_node):
                    path_list.append(path)
                path_length.append(len(path_list[i]))
                weight_avg_list.append(path_average_weight(graph, path_list[i]))
                i += 1

        select_best_path(graph,
                         path_list,
                         path_length,
                         weight_avg_list,
                         delete_entry_node=False,
                         delete_sink_node=True)
    return graph



def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()

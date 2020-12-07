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
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Maelle Boucher"
__copyright__ = "Universite CYTech"
__credits__ = ["Maelle Boucher"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Maelle Boucher"
__email__ = "maelle.boucher09@gmail.com"
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
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
        """
        prend un seul argument correspondant au fichier fastq
        retourne un générateur de séquences
    """
        with open(fastq_file, "r") as filin:
            for line in filin:
                yield next(filin).strip()
                next(filin)
                next(filin)



def cut_kmer(read, kmer_size):
        """
            prend une séquence, une taille de k-mer
            retourne un générateur de k-mer
    """
        for i in range(len(read)-(kmer_size-1)):
            yield read[i:i+kmer_size]



def build_kmer_dict(fastq_file, kmer_size):
    """
    Contruction d'un dict de kmers
    prend un fichier fastq, une taille k- mer et retourne un dictionnaire ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer
    """
        kmer_dict = {}
        kmer_list = []
        read_list = read_fastq(fastq_file)
        for read in read_list:
                kmer_list = kmer_list + (list(cut_kmer(read, kmer_size)))
        for kmer in set(kmer_list):
                kmer_dict[kmer] = kmer_list.count(kmer)
        return kmer_dict



def build_graph(kmer_dict):
        """
        Creation d'un graph a partir des kmers
    Les noeuds correspondent aux prefix et au suffix des kmers
        Retourne : un graph
    """
        graph = nx.DiGraph()
        for kmer in kmer_dict:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph.add_node(prefix)
        graph.add_node(suffix)
        graph.add_edge(prefix, suffix, weight = kmer_dict[kmer])
    return graph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    pass

def get_sink_nodes(graph):
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    pass

def save_contigs(contigs_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    main()

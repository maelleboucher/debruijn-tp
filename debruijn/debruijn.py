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


# Création du graphe de de Bruijn
# Identification des k-mer unique 
def read_fastq(fastq_file):
        """
        Prend un seul argument correspondant au fichier fastq
        Retourne un générateur de séquences
    """
        if(isfile(fastq_file)): 
            with open(fastq_file, "r") as filin:
                for line in filin:
                    yield next(filin).strip()
                    next(filin)
                    next(filin)


def cut_kmer(read, kmer_size):
        """
        Prend une séquence, une taille de k-mer
        Retourne un générateur de k-mer
    """
        for i in range(len(read)-(kmer_size-1)):
            yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
        """
    Contruction d'un dict de kmers
    Prend un fichier fastq, une taille k-mer et retourne un dictionnaire ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer
    """
        kmer_dict = {}
        kmer_list = []
        read_list = read_fastq(fastq_file)
        for read in read_list:
                kmer_list = kmer_list + (list(cut_kmer(read, kmer_size)))
        for kmer in set(kmer_list):
                kmer_dict[kmer] = kmer_list.count(kmer)
        return kmer_dict


#Construction de l’arbre de de Bruijn
def build_graph(kmer_dict):
        """
        Creation d'un arbre a partir des kmers, les noeuds correspondent aux prefix et au suffix des kmers
    """
        graph = nx.DiGraph()
        for kmer in kmer_dict:
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph.add_node(prefix)
            graph.add_node(suffix)
            graph.add_edge(prefix, suffix, weight = kmer_dict[kmer])
        return graph


# Simplification du graphe de de Bruijn
# Résolution des bulles
def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
        """ 
        Detection des chemins indesirables et des noeuds a supprimer
    """
        start= 1
        end = 1
        if delete_entry_node:
            start = 0
        if delete_sink_node:
            end = 0
        for path in path_list:
            for i in range(start, len(path)-end):
                graph.remove_node(path[i])
        return graph


def std(data):
        """ 
        Calcul de l'ecart type
    """
        return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=False, delete_sink_node=False):
            """ 
        Selection du meilleur chemin
    """
            if weight_avg_list[0] > weight_avg_list[1]:
                remove_paths(graph, [path_list[1]], delete_entry_node, delete_sink_node)
            elif weight_avg_list[0] < weight_avg_list[1]:
                remove_paths(graph, [path_list[0]], delete_entry_node, delete_sink_node)
            else:
                if path_length[0] > path_length[1]:
                    remove_paths(graph, [path_list[1]], delete_entry_node, delete_sink_node)
                elif path_length[0] < path_length[1]:
                     remove_paths(graph, [path_list[0]], delete_entry_node, delete_sink_node)
                else:
                     remove_paths(graph, [path_list[random.randint(0,1)]], delete_entry_node, delete_sink_node)
            return graph



def path_average_weight(graph, path):
        """ 
        Calcul d'un poids moyen dans un graph
    """
        weight = []
        for i in range(len(path)-1):
            weight.append(graph.get_edge_data(path[i],
            path[i+1])['weight'])
        return statistics.mean(weight)


def solve_bubble(graph, ancestor_node, descendant_node):
        """ 
        Supression d'une bulle
    """
        paths = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
        while len(paths) >= 2:
            path_size_list = [len(paths[0]), len(paths[1])]
            average_weight_list = [path_average_weight(graph, paths[0]),
                               path_average_weight(graph, paths[1])]
        graph = select_best_path(graph, paths, path_size_list, average_weight_list)
        paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
        return graph


def simplify_bubbles(graph):
        """ 
        Nettoyage d'un graph et suppression de bulles
    """
        starting_nodes = get_starting_nodes(graph)
        ending_nodes = get_sink_nodes(graph)
        for starting_node in starting_nodes:
            for ending_node in ending_nodes:
                 paths = list(nx.all_simple_paths(graph, source=starting_node, target=ending_node))
                 while len(paths) > 1:
                     for i, node in enumerate(paths[0]):
                         if not node in paths[1]:
                            ancetre = paths[0][i - 1]
                            break
                     for node in paths[0][paths[0].index(ancetre)+1:]:
                        if node in paths[1]:
                            descendant_node = node
                        break
                     graph = solve_bubble(graph, ancetre, descendant_node)
                 paths = list(nx.all_simple_paths(graph, source=starting_node, target=ending_node))
        return graph



#Détection des pointes 
def solve_entry_tips(graph, starting_nodes):
    """
    Retourne graphe sans chemin d’entrée indésirable
    """
    


def solve_out_tips(graph, ending_nodes):
    """
    Retourne graphe sans chemin de sortie indésirable
    """


# Parcours du graphe de De Bruijn
def get_starting_nodes(graph):
        """
        Obtention de noeuds d'entree sans predecesseurs
        Retourne : une liste de noeuds d'entree
    """
        starting_nodes = []
        for node in graph.nodes:
            if not graph.in_edges(node):
                starting_nodes.append(node)
        return starting_nodes


def get_sink_nodes(graph):
    """
    Obtention de noeuds de sortie sans successeurs
    Retourne : une liste de noeuds de sortie
    """
    ending_nodes = []
    for node in graph.nodes:
        if not graph.out_edges(node):
            ending_nodes.append(node)
    return ending_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
        """
        Creation d'une liste de contigs en concatenant les kmers d'un chemin
        Retourne d'une lite de tuple (contig, longueur du contig)
    """
        contigs = []
        for starting_node in starting_nodes:
            for ending_node in ending_nodes:
                for paths in nx.all_simple_paths(graph, source=starting_node, target=ending_node):
                    contig  = paths[0]
                    for path in paths[1:]:
                        contig  = contig  + path[-1]
                        contigs.append( (contig, len(contig)) )
        return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """
     Ecrit un fichier de sortie contenant les contigs selon le format fasta
    """
    with open(output_file, "w") as final_file:
        for contig in contigs_list:
            final_file.write(">contig_Numero{} len = {}".format(contigs_list.index(contig),contig))
            final_file.write(fill(contig[0]))



#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    file = read_fastq(args.fastq_file)

    # Lecture du fichier et construction du graphe
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    # Résolution des bulles
    graph = simplify_bubbles(graph)

    # Résolution des pointes d’entrée et de sortie
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)

    # Ecriture du/des contigs
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    list_contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(list_contigs, args.output_file)

if __name__ == '__main__':
    main()
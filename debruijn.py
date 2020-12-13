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
from networkx import *
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
from statistics import stdev
__author__ = "Marcel Mounsi"
__copyright__ = "Cy Tech"
__credits__ = ["Marcel Mounsi"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Marcel Mounsi"
__email__ = "mounsimarc@eisti.eu"
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
	with open(fastq_file,"r") as fichier:
		for i in fichier:
			yield next(fichier).strip()
			next(fichier)
			next(fichier)


def cut_kmer(read, kmer_size):
	for i in read:
		yield read[read.index(i):read.index(i)+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
	D = {}
	L = []
	for element in read_fastq(fastq_file):
		for thing in cut_kmer(element,kmer_size):
			L.append(thing)
	for i in L:
		D[i]=L.count(i)-1
	return D


def build_graph(kmer_dict):
	L=[]
	G = nx.DiGraph()
	for key,value in kmer_dict.items():
		if key[:-1] not in L:
			G.add_node(key[:-1])
		if key[1:] not in L:
			G.add_node(key[1:])
		G.add_edge(key[:-1],key[1:],weight=value)
	return G



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
	i=0
	for element in path_list:
		while i<len(element)-1:
			graph.remove_edge(element[i],element[i+1])
			i=i+1
	if delete_entry_node == True:
		for element in path_list:
			graph.remove_node(element[0])
	if delete_sink_node == True:
		for element in path_list:
			graph.remove_node(element[-1])
	graph.remove_nodes_from(list(nx.isolates(graph)))
	return graph

def std(data):
    return stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
	new_graph = remove_paths(graph,[], delete_entry_node, delete_sink_node)
	#get only the heaviest paths
	L=[]
	for element in weight_avg_list:
  		if element != max(weight_avg_list):
			new_graph = remove_paths(graph, [path_list[weight_avg_list.index(element)]], delete_entry_node, delete_sink_node)
			L.append(element)    			
			del path_length[weight_avg_list.index(element)]
	for element in L:
		del path_list[weight_avg_list.index(element)]
	#get only the longest path
	for thing in path_length:
		if thing != max(path_length):
			new_graph = remove_paths(graph, [path_list[path_length.index(thing)]], delete_entry_node, delete_sink_node)
			del path_list[path_length.index(thing)]
	#random choice
	while len(path_list)>1:
		r = random.choice(path_list)
		new_graph = remove_paths(graph,r, delete_entry_node, delete_sink_node)
		del path_list[weight_avg_list.index(element)]
	return new_graph

def path_average_weight(graph, path):
	i=0
	L=[]
	while i<len(path)-1:
  		L.append(graph[path[i]][path[i+1]]["weight"])
  		i=i+1
	return (sum(L)/len(L))

def solve_bubble(graph, ancestor_node, descendant_node):
	path_list =[]
	path_length = []
	weight_avg_list = []
	for path in nx.all_simple_paths(graph,ancestor_node, descendant_node):
		path_list.append(path)
	for element in path_list:
		path_length.append(len(element))
	 	weight_avg_list.append(path_average_weight(graph,element))
	graph = select_best_path(graph,path_list,path_length,weight_avg_list)
	#graph_3 = select_best_path(graph_3, [[2, 4, 5], [2, 8, 9, 5]],[1, 4], [13, 10])
	return graph

def simplify_bubbles(graph):
	des=[]
	anc=[]
	for element in graph.nodes:
		if len(list(graph.predecessors(element)))>1:
			des.append(element)
		if len(list(graph.successors(element)))>1:
			anc.append(element)
	for u in des:
		for v in anc:
			solve_bubble(graph,v,u)
	return graph

def solve_entry_tips(graph, starting_nodes):
	path_list = []
	path_length = []
	weight_avg_list = []
	#on va chercher tous les chemins du graph
	for element in starting_nodes:
		for thing in get_sink_nodes(graph):
			for path in nx.all_simple_paths(graph, source=element, target=thing):
				path_list.append(path)
	#on associe aux chemins les poids et les longueurs 
	for element in path_list:
		path_length.append(len(element))
		weight_avg_list.append(path_average_weight(graph,element))
	select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=False, delete_sink_node=False)

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
	L = []
	for element in graph.nodes:
  		if list(graph.predecessors(element)) == []:
    			L.append(element)
	return L

def get_sink_nodes(graph):
	L = []
	for element in graph.nodes:
  		if list(graph.successors(element)) == []:
    			L.append(element)
	return L

def get_contigs(graph, starting_nodes, ending_nodes):
	L = []
	B = []
	for input in starting_nodes:
  		for output in ending_nodes:
    			for element in next(nx.all_simple_paths(graph,source = input, target = output)):
      				L.append(element[:-1])
			L.append(output[1:])
    			new = ''.join(L)
    			L[:] = []
    			B.append((new,len(new)))
	return B

def save_contigs(contigs_list, output_file):
 	output_file = 'eva71_hundred_reads.fq'
	fichier = open(output_file,"w")
 	i = 0
	for element in contigs_list:
		fichier.write(fill(element[0]+"\n",80))
		fichier.write(fill(">"+"contig_"+str(i)+" "+str(element[1])+"\n",80))
		i = i+1


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

#==============================================================
# Main program
#==============================================================
def main():
	"""
	Main program function
	"""
	args = get_arguments() 	
	dico = build_kmer_dict(args.fastq_file, args.kmer_size)
	new_graph = build_graph(dico)
	start_nodes = get_starting_nodes(new_graph)
	end_nodes = get_sink_nodes(new_graph)
	contigs = get_contigs(new_graph,start_nodes,end_nodes)
	save_contigs(contigs,"debruijn/contigs.txt")

if __name__ == '__main__':
    main()

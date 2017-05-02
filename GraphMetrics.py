from __future__ import division
import numpy as np
import pandas as pd
from igraph import *
import time
from collections import defaultdict
import itertools
import datetime

def metrics(file_path):
	ori_igraph = Graph.Read_Ncol(file_path, directed=True)
	n_ori = ori_igraph.vcount()
	e_ori = ori_igraph.ecount()

	er_igraph = Graph.Erdos_Renyi(n=n_ori, p=e_ori/((n_ori*(n_ori-1))/2), directed=False, loops=False)
	n_er = er_igraph.vcount()
	e_er = er_igraph.ecount()

	ff_igraph = Graph.Forest_Fire(n=n_ori, fw_prob=0.354, bw_factor=0.56, ambs=3, directed=False)
	n_ff = ff_igraph.vcount()
	e_ff = ff_igraph.ecount()

	ff1_igraph = Graph.Forest_Fire(n=int(n_ori/10), fw_prob=0.354, bw_factor=0.56, ambs=3, directed=False)
	n_ff1 = ff1_igraph.vcount()
	e_ff1 = ff1_igraph.ecount()

	ori_density = ori_igraph.density(loops=False)
	er_density = er_igraph.density(loops=False)
	ff_density = ff_igraph.density(loops=False)
	ff1_density = ff1_igraph.density(loops=False)

	ori_average_path_length = ori_igraph.average_path_length(directed=False, unconn=True)
	er_average_path_length = er_igraph.average_path_length(directed=False, unconn=True)
	ff_average_path_length = ff_igraph.average_path_length(directed=False, unconn=True)
	ff1_average_path_length = ff1_igraph.average_path_length(directed=False, unconn=True)

	ori_global_clustering_coefficient = ori_igraph.transitivity_undirected(mode="zero")
	er_global_clustering_coefficient = er_igraph.transitivity_undirected(mode="zero")
	ff_global_clustering_coefficient = ff_igraph.transitivity_undirected(mode="zero")
	ff1_global_clustering_coefficient = ff1_igraph.transitivity_undirected(mode="zero")

	np.set_printoptions(threshold='nan')
	utc_datetime = datetime.datetime.utcnow()
	formatted_string = utc_datetime.strftime("%Y-%m-%d-%H%MZ")
	filename = 'metrics_flickr-edges.txt'

	with open(filename, "a") as myfile:
		myfile.write("Original Graph")
		myfile.write("\nNodes: %s" %str(n_ori))
		myfile.write("\nEdges: %s" %str(e_ori))
		myfile.write("\nThe density of the graph is %f" %ori_density)
		myfile.write("\nThe average path length of the graph is %f" %ori_average_path_length)
		myfile.write("\nThe global clustering coefficient of the graph is %f" %ori_global_clustering_coefficient)
		myfile.write("\n\nErdos-Renyi Graph")
		myfile.write("\nNodes: %s" %str(n_er))
		myfile.write("\nEdges: %s" %str(e_er))
		myfile.write("\nThe density of the graph is %f" %er_density)
		myfile.write("\nThe average path length of the graph is %f" %er_average_path_length)
		myfile.write("\nThe global clustering coefficient of the graph is %f" %er_global_clustering_coefficient)
		myfile.write("\n\nForest Fire Graph")
		myfile.write("\nNodes: %s" %str(n_ff))
		myfile.write("\nEdges: %s" %str(e_ff))
		myfile.write("\nThe density of the graph is %f" %ff_density)
		myfile.write("\nThe average path length of the graph is %f" %ff_average_path_length)
		myfile.write("\nThe global clustering coefficient of the graph is %f" %ff_global_clustering_coefficient)
		myfile.write("\n\nSecond Forest Fire Graph - FF/10")
		myfile.write("\nNodes: %s" %str(n_ff1))
		myfile.write("\nEdges: %s" %str(e_ff1))
		myfile.write("\nThe density of the graph is %f" %ff1_density)
		myfile.write("\nThe average path length of the graph is %f" %ff1_average_path_length)
		myfile.write("\nThe global clustering coefficient of the graph is %f" %ff1_global_clustering_coefficient)
	return "Success"

original_graph_path = "/home/03962/tg832615/python/graphs/soc-Slashdot0811.txt"
start_time = time.time()
metrics(original_graph_path)
print("--- %s seconds ---" % (time.time() - start_time))

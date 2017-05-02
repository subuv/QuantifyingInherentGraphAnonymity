from __future__ import division
import numpy as np
import pandas as pd
from igraph import *
import time
import networkx as nx
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
from collections import defaultdict
import itertools
import datetime
from multiprocessing.pool import ThreadPool
import csv
import editdistance

linkage_covariance_matrix_dict = {}
# %matplotlib inline

def linkage_covariance(ori_igraph, order_no = 1, topk = 10, num_processes = 4, graph_name='Original'):
	n_ori = ori_igraph.vcount()
	e_ori = ori_igraph.ecount()
	d1 = defaultdict(list)
	
	ori_edgelist = ori_igraph.get_edgelist()

	for k, v in ori_edgelist:
		d1[k].append(v)
		d1[v].append(k)

	global neisets_ori 
	neisets_ori = dict((k, tuple(v)) for k, v in d1.iteritems())

	nodes = neisets_ori.keys()

	global nodes_array
	nodes_array = np.array(nodes)

	global linkage_covariance_matrix_dict
	linkage_covariance_matrix_dict = {}

	global n
	n = ori_igraph.vcount()

	pool = ThreadPool(processes=num_processes)

	pool.map(linkage_covariances, split_list(nodes_array.tolist(), int(n/num_processes-1)))

	pool.close()

	pool.join()
	topk_val = int(topk)/100
	linkage_covariance_matrix_ori = {}
	
	for k, v in linkage_covariance_matrix_dict.iteritems():
		linkage_covariance_matrix_ori[k] = '-'.join((np.array(v)[np.argsort(np.array(v))[::-1]]).tolist())

	results_list = pd.DataFrame(
		{'lcm': linkage_covariance_matrix_ori.values()
		})

	identical = results_list.groupby('lcm').size()
	
	identical_group = pd.DataFrame(
		{'lcm_groups': identical.values
		})
	identical_grp = identical_group.groupby('lcm_groups').size()
	
	group_x, group_y = getGroup(identical_grp)
	
	percent_unique = 0

	if group_x[0] == 1:
		percent_unique = 100*(group_y[0]/n_ori)
	
	np.set_printoptions(threshold='nan')
	utc_datetime = datetime.datetime.utcnow()
	formatted_string = utc_datetime.strftime("%Y-%m-%d-%H%MZ")
	filename = 'output_%s_%s_%s.txt'%('sweden', graph_name, formatted_string)
	
	with open(filename, "a") as myfile:
		myfile.write("\nThere are %d nodes and %d edges in the graph" %(n_ori, e_ori))
		myfile.write("\n%.3f%% of nodes have unique linkage covariance values" %(percent_unique))
		myfile.write("\n%.3f%% of nodes have identical linkage covariance values" %(100-percent_unique))
		myfile.write("\n\nIdentical Group:\n %s"%identical_grp)
		myfile.write("\n\nGroupX and GroupY: \n %s, %s"%(group_x, group_y))

	return group_x, group_y

def split_list(lcm_ori_list, N):
	result_list=[]
	for i in range(0, len(lcm_ori_list), N):
		result_list.append(lcm_ori_list[i:i+N])
	return result_list

def linkage_covariances(nodes):
	for node in nodes:
		neis = nodes_array[np.where(nodes_array > node)].tolist()
		if node not in linkage_covariance_matrix_dict.keys():
			linkage_covariance_matrix_dict[node] = []
		for nei in neis:
			if nei not in linkage_covariance_matrix_dict.keys():
				linkage_covariance_matrix_dict[nei] = []
			common_neis_ori = set(neisets_ori[node]) & set(neisets_ori[nei])
			lcm = len(common_neis_ori)
			if lcm > 0:
				linkage_covariance_matrix_dict[node].append(str(lcm))
				linkage_covariance_matrix_dict[nei].append(str(lcm))

def getGroup(identical_grp):
	group_axes = []
	for k, v in sorted(identical_grp.iteritems()):
		group_axes.append((int(k), int(v)))

	return np.array(sorted(group_axes)).T

def main():
	start_time = time.time()
	original_graph_path = "/home/03962/tg832615/python/graphs/terrornet_edges"
	ori_igraph = Graph.Read_Ncol(original_graph_path, directed=True)
        n_ori = ori_igraph.vcount()
        e_ori = ori_igraph.ecount()
        order = 1
        topk = 1
        num_processes = 50
        identical_lc_ori = linkage_covariance(ori_igraph, order, topk, num_processes)

	#Forest Fire Normal
	ff_nxgraph = Graph.Forest_Fire(n=n_ori, fw_prob = 0.35, bw_factor = 0.56, ambs = 2, directed=False)
	n_ff = ff_nxgraph.vcount()
	e_ff = ff_nxgraph.ecount()
	group_xff, group_yff = linkage_covariance(ff_nxgraph, order, topk, num_processes, 'FF')

	#Forest Fire FF/10
	ff1_nxgraph = Graph.Forest_Fire(n=int(n_ori/10), fw_prob = 0.354, bw_factor = 0.56, ambs = 3, directed=False)
	n_ff1 = ff1_nxgraph.vcount()
	e_ff1 = ff1_nxgraph.ecount()
	group_xff1, group_yff1 = linkage_covariance(ff1_nxgraph, order, topk, num_processes, 'FFBy10')
	
	print("--- %s seconds ---" % (time.time() - start_time))
main()

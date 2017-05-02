from __future__ import division
import numpy as np
import pandas as pd
import logging
from igraph import *
import time
import matplotlib.pyplot as plt
from collections import defaultdict
import itertools
import datetime
import copy
from multiprocessing.pool import ThreadPool
from scipy.stats.stats import pearsonr
from scipy.stats import norm
from sklearn import metrics

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def linkage_covariance(file_path, order_no = 1, topk = 10, num_processes = 4):
	ori_igraph = Graph.Read_Ncol(file_path, directed=True)
	logging.info('Parsing the original graph to get the number of nodes')
	
	logging.info('Creating a dummy linkage covariance matrix original')

	d1 = defaultdict(list)
	
	#ori_edgelist = ori_igraph.neighborhood(vertices=None, order=int(order_no), mode="in")
	ori_edgelist = ori_igraph.get_edgelist()

	#for dummy in ori_edgelist:
	#	d1[dummy.pop(0)] = dummy

	for k, v in ori_edgelist:
		d1[k].append(v)

	global neisets_ori 
	neisets_ori = dict((k, tuple(v)) for k, v in d1.iteritems())

	logging.debug("Neighbor Sets: %s" %neisets_ori)

	nodes = neisets_ori.keys()

	global nodes_array
	nodes_array = np.array(nodes)

	global linkage_covariance_matrix_dict
	linkage_covariance_matrix_dict = {}

	global n
	n = ori_igraph.vcount()

	#names = np.array(ori_igraph.vs["name"])
	
	logging.debug("Nodes: %s" %nodes)
	
	logging.info("Number of nodes: %s" %n)
	
	logging.info("Computing the linkage covariance values")

	pool = ThreadPool(processes=num_processes)

	pool.map(linkage_covariances, split_list(nodes, int(len(nodes)/num_processes-1)))
	
	pool.close()

	pool.join()

	linkage_covariance_matrix_ori = linkage_covariance_matrix_dict.values()

	#for i in xrange(0, n):
	#	linkage_covariance_matrix_ori[i] = np.array(linkage_covariance_matrix_ori[i])[np.argsort(linkage_covariance_matrix_ori[i])[::-1]]

	logging.debug("Linkage covariance:\n %s" %linkage_covariance_matrix_ori)

	logging.info("Linkage covariance done")

	metrics={'identical': {}, 'similar': {}}
	#variance = np.var(linkage_covariance_matrix_ori)

	for v1 in nodes:
		neis = nodes_array[np.where(nodes_array > v1)].tolist()
		for v2 in neis:
			l1 = linkage_covariance_matrix_dict[v1]
			l2 = linkage_covariance_matrix_dict[v2]
			if(len(l1) == len(l2)):
				l1 = np.array(l1)[np.argsort(l1)[::-1]]
				l2 = np.array(l1)[np.argsort(l1)[::-1]]
				array_equal = np.array_equal(l1, l2)
				array_close = np.allclose(np.array(l1),np.array(l2))
				if array_equal:
					if sum(l1) not in metrics['identical']:
						metrics['identical'][sum(l1)] = []
					metrics['identical'][sum(l1)].append((v1, v2))
				if array_close:
					if sum(l1) not in metrics['similar']:
						metrics['similar'][sum(l1)] = []
					metrics['similar'][sum(l1)].append((v1, v2))

	logging.info("The identical linkage covariance values are %s" %metrics['identical'].values())
	logging.info("The similar linkage covariance values are %s" %metrics['similar'].values())
	
	identical_group = []
	identical_nodes = []
	similar_group = []
	similar_nodes = []
	for k in metrics['identical'].values():
		identical_group.append(len(k))
		for v1, v2 in k:
			identical_nodes.append(v1)
			identical_nodes.append(v2)
	for k in metrics['similar'].values():
		similar_group.append(len(k))
		for v1, v2 in k:
			similar_nodes.append(v1)
			similar_nodes.append(v2)

	# t = set([i for sub in metrics['identical'] for i in sub])
	percent_identical = 100*(len(set(identical_nodes))/n)

	# t = set([i for sub in metrics['similar'] for i in sub])
	percent_similar = 100*(len(set(similar_nodes))/n)

	logging.info("%.3f%% of nodes have identical linkage covariance values" %percent_identical)
	logging.info("%.3f%% of nodes have similar linkage covariance values" %percent_similar)

	results_list = pd.DataFrame(
		{'identical_group_size': identical_group,
		 'similar_group_size': similar_group
		})
    
	#Identical Nodes in linkage covariance values
	identical = results_list.groupby('identical_group_size').size()
	identical = np.array(identical.values)
	identical = np.insert(identical, 0,0)
	print identical

	fig, ax = plt.subplots(1)

	ax.set_ylim(ymin=0)
	ax.plot(np.cumsum(identical)/sum(identical), color='k', alpha=0.5, label = "Real Graph")
	ax.legend(loc="best", prop={'size':5})
	#plt.plot(sorted_result_list['linkage_covariance'].values, color='k', alpha=0.5)
	#plt.plot(sorted_result_list['leverage_centrality'].values, color='r', alpha=0.5)
	#plt.plot(sorted_result_list['clustering_coefficient'].values, color='b', alpha=0.5)
	
	return "Success"

def split_list(lcm_ori_list, N):
	result_list=[]
	for i in range(0, len(lcm_ori_list), N):
		result_list.append(lcm_ori_list[i:i+N])
	return result_list

def top_k_indices(a, N):
    return np.argsort(a)[::-1][:N]

def h_dist(file_path, order=1, topk=10, num_processes=4):
	df_rlc_original = linkage_covariance(file_path, order, topk, int(num_processes))
	print df_rlc_original

def linkage_covariances(nodes):
	for node in nodes:
		neis = nodes_array[np.where(nodes_array > node)].tolist()
		if node not in linkage_covariance_matrix_dict:
			linkage_covariance_matrix_dict[node] = []
		for nei in neis:
			if nei not in linkage_covariance_matrix_dict:
				linkage_covariance_matrix_dict[nei] = []
			common_neis_ori = set(neisets_ori[node]) & set(neisets_ori[nei])
			lcm = len(common_neis_ori)/n
			if lcm > 0:
				linkage_covariance_matrix_dict[node].append(lcm)
				linkage_covariance_matrix_dict[nei].append(lcm)

def main():
	print "lcm_revamped: Let us calculate the linkage covariance of the graph"
	original_graph_path = raw_input("Path of the graph file original?\n")
	order = raw_input("Order of hops?\n")
	topk = raw_input("What percentage of hub nodes do you want to see?\n")
	num_processes = raw_input("How many threads do you want to use?\n")
	start_time = time.time()
	h_dist(original_graph_path, order, topk, num_processes)
	print("--- %s seconds ---" % (time.time() - start_time))

main()
plt.show()

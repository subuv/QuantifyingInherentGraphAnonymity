from __future__ import division
import numpy as np
import pandas as pd
import logging, sys
from igraph import *
import time
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
from scipy.stats.stats import pearsonr
from joblib import Parallel, delayed
from scipy.signal import argrelextrema

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def linkage_covariance(file_path, order_no):
	ori_igraph = Graph.Read_Ncol(file_path, directed=True)
	logging.info('Parsing the original graph to get the number of nodes')
	
	n_ori = ori_igraph.vcount()
	e_ori = ori_igraph.ecount()
	ori_igraph_degree = ori_igraph.vs.degree(mode=ALL, loops=False)

	logging.info(n_ori)
	
	er_igraph = Graph.Erdos_Renyi(n=n_ori, m=e_ori, directed=False, loops=False)

	logging.info('Parsing the er graph to get the number of nodes')

	n_er = er_igraph.vcount()

	logging.info(n_er)
	
	random_igraph = Graph.Degree_Sequence(out=ori_igraph_degree, method="simple")

	logging.info('Parsing the random graph to get the number of nodes')

	n_random = random_igraph.vcount()

	logging.info(n_random)
	
	logging.info('Creating a dummy linkage covariance matrix original')

	linkage_covariance_matrix_ori = np.zeros((n_ori, n_ori))
	linkage_covariance_matrix_er = np.zeros((n_er, n_er))
	linkage_covariance_matrix_random = np.zeros((n_random, n_random))

	logging.info('Collecting the neighbor sets of the original graph for the linkage covariance matrix')

	ori_igraph = Graph.Read_Edgelist(file_path, directed=True)

	#ori_igraph.delete_vertices((6,11))
	
	#neisets_ori = [set(ori_igraph.neighbors(i)) for i in xrange(n_ori)]
	neisets_ori = Parallel(n_jobs=2, backend="threading")(delayed(set)(ori_igraph.neighbors(i)) for i in xrange(n_ori))
	#neisets_er = [set(er_igraph.neighbors(i)) for i in xrange(n_er)]
	neisets_er = Parallel(n_jobs=2, backend="threading")(delayed(set)(er_igraph.neighbors(i)) for i in xrange(n_er))
	#neisets_random = [set(random_igraph.neighbors(i)) for i in xrange(n_random)]
	neisets_random = Parallel(n_jobs=2, backend="threading")(delayed(set)(random_igraph.neighbors(i)) for i in xrange(n_er))

	logging.info('Populating the linkage covariance matrix original')
	
	for v1, v2 in np.vstack(np.triu_indices_from(linkage_covariance_matrix_ori,k=0)).T:
		if (v1 != v2):
			common_neis_ori = neisets_ori[v1].intersection(neisets_ori[v2])
			#linkage_covariance_matrix_ori[v1, v2] = (len(common_neis_ori)/n_ori) - (len(neisets_ori[v1]) * len(neisets_ori[v2]))/n_ori**2
			linkage_covariance_matrix_ori[v1, v2] = (len(common_neis_ori)/n_ori)
			linkage_covariance_matrix_ori[v2, v1] = linkage_covariance_matrix_ori[v1, v2]
	
	logging.info('Populating the linkage covariance matrix ER')
	
	for v1, v2 in np.vstack(np.triu_indices_from(linkage_covariance_matrix_er,k=0)).T:
		if (v1 != v2):
			try:
				common_neis_er = neisets_er[v1].intersection(neisets_er[v2])
				linkage_covariance_matrix_er[v1, v2] = (len(common_neis_er)/n_er) - (len(neisets_er[v1]) * len(neisets_er[v2]))/n_er**2
				linkage_covariance_matrix_er[v2, v1] = linkage_covariance_matrix_er[v1, v2]
			except:
				continue

	logging.info('Populating the linkage covariance matrix random')
	
	for v1, v2 in np.vstack(np.triu_indices_from(linkage_covariance_matrix_random,k=0)).T:
		if (v1 != v2):
			try:
				common_neis_random = neisets_random[v1].intersection(neisets_random[v2])
				linkage_covariance_matrix_random[v1, v2] = (len(common_neis_random)/n_random) - (len(neisets_random[v1]) * len(neisets_random[v2]))/n_random**2
				linkage_covariance_matrix_random[v2, v1] = linkage_covariance_matrix_random[v1, v2]
			except:
				continue

	logging.info('Linkage covariance matrix done')
	
	logging.debug('The linkage_covariance_matrix original!\n')
	logging.debug(linkage_covariance_matrix_ori)
	
	logging.debug('The linkage_covariance_matrix ER!\n')
	logging.debug(linkage_covariance_matrix_er)
	
	logging.debug('The linkage_covariance_matrix random!\n')
	logging.debug(linkage_covariance_matrix_random)
	
	logging.info('Sorting each row in linkage covariance matrix in descending order')
	
	for i in range(0, n_ori):
		linkage_covariance_matrix_ori[linkage_covariance_matrix_ori[i, :][::-1].sort()]
	
	for i in range(0, n_er):
		linkage_covariance_matrix_er[linkage_covariance_matrix_er[i, :][::-1].sort()]
	
	for i in range(0, n_random):
		linkage_covariance_matrix_random[linkage_covariance_matrix_random[i, :][::-1].sort()]
	
	logging.info('Sorted the matrix successfully')
	
	logging.debug('RLC Original!\n')
	logging.debug(linkage_covariance_matrix_ori)

	logging.debug('RLC ER!\n')
	logging.debug(linkage_covariance_matrix_er)

	logging.debug('RLC Random!\n')
	logging.debug(linkage_covariance_matrix_random)

	logging.info('Use a sort key and sort the matrix by each column in descending order')
	
	rlc_ori = linkage_covariance_matrix_ori.copy()
	rlc_er = linkage_covariance_matrix_er.copy()
	rlc_random = linkage_covariance_matrix_random.copy()

	rlc_ori = rlc_ori[rlc_ori[:, 0].argsort()][::-1]
	rlc_er = rlc_er[rlc_er[:, 0].argsort()][::-1]
	rlc_random = rlc_random[rlc_random[:, 0].argsort()][::-1]

	logging.info('Sorted the matrix by each column in ascending order')
	
	logging.debug('RLC original after making it a matrix with sorted order!\n')
	logging.debug(rlc_ori)

	logging.debug('RLC ER after making it a matrix with sorted order!\n')
	logging.debug(rlc_er)
	
	logging.debug('RLC Random after making it a matrix with sorted order!\n')
	logging.debug(rlc_random)

	logging.info('Normalize the matrix now!')

	rownorm_ori = np.linalg.norm(rlc_ori, axis=1)
	rownorm_er = np.linalg.norm(rlc_er, axis=1)
	rownorm_random = np.linalg.norm(rlc_random, axis=1)
	
	for i in range(0, n_ori):
		rlc_ori[i,] = rlc_ori[i,]/rownorm_ori[i]
		rlc_er[i,] = rlc_er[i,]/rownorm_er[i]
		rlc_random[i,] = rlc_random[i,]/rownorm_random[i]

	rlc_er[np.isnan(rlc_er)]=0
	rlc_random[np.isnan(rlc_random)]=0

	logging.debug('Normalized sum of squares of link covariances of the original graph!\n')
	logging.debug(rlc_ori)

	logging.debug('Normalized sum of squares of link covariances of the ER graph!\n')
	logging.debug(rlc_er)

	logging.debug('Normalized sum of squares of link covariances of the random graph!\n')
	logging.debug(rlc_random)

	logging.info('Normalized the matrix successfully')

	pearson_distance_ori = list()
	for i in range(1, n_ori):
		if math.isnan(pearsonr(rlc_ori[i,:], rlc_ori[i,:])[0]) == False:
			pearson_distance_ori.append(pearsonr(rlc_ori[0,:], rlc_ori[i,:])[0])

	pearson_distance_er = list()
	for i in range(1, n_er):
		if math.isnan(pearsonr(rlc_ori[i,:], rlc_er[i,:])[0]) == False:
			pearson_distance_er.append(pearsonr(rlc_er[0,:], rlc_er[i,:])[0])

	pearson_distance_random = list()
	for i in range(1, n_random):
		if math.isnan(pearsonr(rlc_ori[i,:], rlc_random[i,:])[0]) == False:
			pearson_distance_random.append(pearsonr(rlc_random[0,:], rlc_random[i,:])[0])

	np.set_printoptions(threshold='nan')

	logging.info("The pearson distance between first row and various rows of the original graph!\n")
	logging.info(np.sort(pearson_distance_ori))

	logging.info("The pearson distance between first row and various rows of the ER graph!\n")
	logging.info(np.sort(pearson_distance_er))

	logging.info("The pearson distance between first row and various rows of the Random graph!\n")
	logging.info(np.sort(pearson_distance_random))

	with open("output.txt", "a") as myfile:
		myfile.write("Original\n")
		myfile.write(str(np.sort(pearson_distance_ori)))
		myfile.write("\nER\n")
		myfile.write(str(np.sort(pearson_distance_er)))
		myfile.write("\nRandom\n")
		myfile.write(str(np.sort(pearson_distance_random)))

	return "Success"

def h_dist(file_path, order):
	df_rlc_original = linkage_covariance(file_path, order)
	print df_rlc_original

def main():
	print "rankedlinkcharacteristicvectorL2normed_er_random: Let us calculate the linkage covariance of the graph"
	original_graph_path = raw_input("Path of the graph file original?\n")
	order = raw_input("Order of hops?\n")
	start_time = time.time()
	h_dist(original_graph_path, order)
	print("--- %s seconds ---" % (time.time() - start_time))

main()
# plt.show()

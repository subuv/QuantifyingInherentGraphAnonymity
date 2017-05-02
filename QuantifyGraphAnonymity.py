from __future__ import division
import numpy as np
import pandas as pd
import logging, sys
from igraph import *
import time
#import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

def linkage_covariance(file_path, order_no):
	terrorist_igraph = Graph.Read_Ncol(file_path, directed=True)
	
	logging.info('Parsing the graph to get the number of nodes')
	n = terrorist_igraph.vcount()

	logging.info(n)
	
	#LC <- matrix(0,n,n)
	# matrix of zeroes
	logging.info('Creating a dummy linkage covariance matrix')
	linkage_covariance_matrix = np.zeros((n, n))
	logging.info('Collecting the neighbor sets for the linkage covariance matrix')
	neisets = [set(terrorist_igraph.neighbors(i)) for i in xrange(n)]
	logging.info('Populating the linkage covariance matrix')
	for v1, v2 in terrorist_igraph.get_edgelist():
		common_neis = neisets[v1].intersection(neisets[v2])
		linkage_covariance_matrix[v1, v2] = (len(common_neis)/n) - (len(neisets[v1]) * len(neisets[v2]))/n**2
		linkage_covariance_matrix[v2, v1] = linkage_covariance_matrix[v1, v2]

	for v1, v2 in np.transpose(np.nonzero(linkage_covariance_matrix == 0)):
		if v1 != v2:
			common_neis = neisets[v1].intersection(neisets[v2])
			linkage_covariance_matrix[v1, v2] = (len(common_neis)/n) - (len(neisets[v1]) * len(neisets[v2]))/n**2

	logging.info('Linkage covariance matrix done')
	
	logging.debug('The linkage_covariance_matrix!\n')
	logging.debug(linkage_covariance_matrix)
	
	# create sorted LC matrix, first sort within each row
	#RLC.o <- matrix(0,n,n)
	#rownames(RLC.o) <- c(1:63)
	#for (i in 1:n) {RLC.o[i,] <- sort(LC[i,],decreasing=T)}

	logging.info('Sorting each row in linkage covariance matrix in descending order')
	
	rlc = np.sort(linkage_covariance_matrix, axis=1)
	rlc = np.fliplr(rlc)

	logging.info('Sorted the matrix successfully')
	
	logging.debug('RLC!\n')
	logging.debug(rlc)

	df_rlc = pd.DataFrame(rlc)
	
	logging.debug('RLC Data Framed!\n')
	logging.debug(df_rlc)
	
	logging.info('Use a sort key and sort the matrix by each column in descending order')
	
	rlc = rlc[rlc[:, 0].argsort()][::-1]

	logging.debug('RLC after making it a matrix with sorted order!\n')
	logging.debug(rlc)

	logging.info('Sorted the matrix by each column in ascending order')
	
	#now normalize by the L2 modulus
	#for (i in 1:n) {
	#  RLC.o[i,] <- RLC.o[i,]/sqrt(sum(RLC.o[i,]^2))
	#}

	logging.info('Normalize the matrix now!')

	for i in range(0, n):
		rlc[i,] = rlc[i,]/np.sqrt(np.sum(np.square(rlc[i,])))

	logging.debug('Normalized sum of squares of link covariances of the input graph!\n')
	logging.debug(rlc)

	logging.info('Normalized the matrix successfully')

	rlc_pd_df = pd.DataFrame(data=rlc, index=None, columns=None)
	
	return rlc_pd_df

def h_dist(file_path1, file_path2, order):
	df_rlc_original = linkage_covariance(file_path1, order)
	print df_rlc_original
	df_rlc_anon = pd.DataFrame(data=linkage_covariance(file_path2, order), index=df_rlc_original.index.copy(), columns=df_rlc_original.index.copy())
	df_rlc_anon.fillna(0, inplace=True)
	print df_rlc_anon
	pearson_coefficient = np.corrcoef(df_rlc_original[:,0], df_rlc_anon[:, 0])
	#df_rlc = df_rlc_anon.subtract(df_rlc_original, fill_value=0)
	#df_rlc = df_rlc.sum(axis=0)
	#df_rlc = df_rlc.sum()
	#df_rlc_abs = df_rlc.abs()
	#print df_rlc
	print pearson_coefficient
	#plt.figure()
	#df_rlc_abs.plot()

def main():
	print "rankedlinkcharacteristicvectorL2normed_el: Let us calculate the linkage covariance of the graph"
	original_graph_path = raw_input("Path of the graph file 1?\n")
	anon_graph_path = raw_input("Path of the graph file 2?\n")
	order = raw_input("Order of hops?\n")
	start_time = time.time()
	h_dist(original_graph_path, anon_graph_path, order)
	print("--- %s seconds ---" % (time.time() - start_time))

main()
#plt.show()

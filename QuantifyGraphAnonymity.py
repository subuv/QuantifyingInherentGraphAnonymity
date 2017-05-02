import numpy as np
import pandas as pd
import logging, sys
from igraph import *
import time
import matplotlib.pyplot as plt
from pympler.tracker import SummaryTracker
import multiprocessing as mp
import itertools
from scipy import spatial

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

num_cpu = mp.cpu_count()
print ("The number of processors is: ", num_cpu)

def linkage_covariance(file_path, order_no):
	terrorist_igraph = Graph.Read(file_path, format="edgelist", directed=True)
	es = EdgeSeq(terrorist_igraph)
	logging.info('Got the file details. Reading file...')

	logging.info('Getting the neighborhood for the order...')
	#tracker = SummaryTracker()
	#notn = terrorist_igraph.neighborhood(vertices=None, order=int(order_no), mode="all")
	#tracker.print_diff()
	logging.info('Got the neighborhood for the order successfully')
	#logging.debug(notn)

	logging.info('Stringing up the edgelist from the neighborhood')
	terror_edge_array = [e.tuple for e in terrorist_igraph.es]
	#tuple(terror_edge_array)
	#print terror_edge_array
	logging.info('Stringed up the edgelist from the neighborhood successfully')
	logging.debug(terror_edge_array)
	
	logging.info('Creating a graph from the edgelist of neighborhood')
	terrorist_igraph = Graph.TupleList(directed=True, edges = terror_edge_array)
	logging.info('Created a graph from the edgelist of neighborhood successfully')
	
	logging.info('Computing the Adjacency matrix of the input graph')
	terrorist_adjacency = terrorist_igraph.get_adjacency()
	logging.info('Computed the Adjacency matrix of the input graph successfully')

	# print the header row
	logging.debug('Adjacency Matrix of the graph!\n')
	logging.debug(terrorist_adjacency)

	logging.info('Adjacency matrix in readable type')
	d = np.array(terrorist_adjacency.data)
	d = d.astype(float)

	#hop = order in call below
	# ax <- neighborhood(g,order=1,neighborhood.type="total",mode="graph",partial=T)
 	
 	#Use the shape command to get the dimensions of the matrix
	logging.info('Parsing the graph to get the number of nodes')
	[m, n] = d.shape
	# header_row = list()
	# logging.info('Creating a header row for the number of nodes')
	# for i in range(0,n):
	# 	header_row.append(i+1)
	
	logging.info('Graph parsed successfully and header row created!')
	
	# print the header row
	# logging.debug('Header Row!\n')
	# logging.debug(header_row)

	logging.info('Insert the header row above')
	#colnames(ax) <- c(1:63)
	input_graph = d
	#print input_graph
	# TO-DO ------ not working!!! - Start
	#np.insert(input_graph, 1, header_row, axis=0)
	# TO-DO ------ not working!!! - End

	logging.info('Header row inserted!!')
	
	logging.info('Finding the sum of all the values in each row')
	#ax.r <- rowSums(ax)
	input_graph_rowsums = input_graph.sum(axis=1)

	logging.info('Found the sum!')
	logging.debug(input_graph_rowsums)
	#LC <- matrix(0,n,n)
	# matrix of zeroes
	linkage_covariance_matrix = np.zeros((n, n))
	
	logging.info('Find the linkage covariance and compute the linkage covariance matrix')

	# arglist = []
	# try :
	#     arglist.append((linkage_covariance_matrix, input_graph, input_graph_rowsums, n))
	#     pool = mp.Pool(num_cpu)
	#     result = pool.apply_async(compute_lcm, args=arglist)
	# except Exception as e:
	# 	logging.error(e)
	# 	logging.debug("There was some exception when parallelizing the matrix computation")
	# 	pool.terminate()
	# finally:
	# 	pool.close()
	# 	pool.join()

	# logging.info('Linkage covariance matrix done successfully by parallelization')

	# linkage_covariance_matrix = result.get()
	linkage_covariance_matrix = np.cov(input_graph, numpy.transpose(input_graph), ddof=0)[0][1]
	np.fill_diagonal(linkage_covariance_matrix, 0)
	
	logging.debug('The linkage_covariance_matrix!\n')
	logging.debug(linkage_covariance_matrix)
	logging.debug('The input_graph!\n')
	logging.debug(input_graph)

	# create sorted LC matrix, first sort within each row
	#RLC.o <- matrix(0,n,n)
	#rownames(RLC.o) <- c(1:63)
	#for (i in 1:n) {RLC.o[i,] <- sort(LC[i,],decreasing=T)}

	logging.info('Sorting each row in linkage covariance matrix in descending order')
	
	idx = linkage_covariance_matrix[:, 0].argsort()
	rlc = np.take(linkage_covariance_matrix, idx, axis=1)
	#rlc = np.sort(linkage_covariance_matrix, axis=1)
	rlc = np.fliplr(rlc)

	logging.info('Sorted the matrix successfully')
	
	logging.debug('RLC!\n')
	logging.debug(rlc)

	sort_key = list()
	#ascending_key = list()
	for i in range(0,n):
		sort_key.append(i)
		idx = rlc[i,:].argsort()
		#ascending_key.append(False)

	rlc = np.take(rlc, idx, axis=0)
	df_rlc = pd.DataFrame(rlc)
	
	logging.debug('RLC Data Framed!\n')
	logging.debug(df_rlc)
	
	logging.info('Use a sort key and sort the matrix by each column in ascending order')
	
	# idx = rlc[:, ].argsort()
	# print idx
	# rlc = np.take(rlc, idx, axis=1)
	
	#rlc[rlc[:, sort_key].argsort()]

	# df_rlc.sort_values(sort_key, ascending=ascending_key, inplace=True)
	# rlc = df_rlc.as_matrix()
	
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
	logging.debug(df_rlc)

	logging.info('Normalized the matrix successfully')

	rlc_pd_df = pd.DataFrame(data=rlc, index=None, columns=None)
	
	return rlc_pd_df

def h_dist(file_path1, file_path2, order):
	df_rlc_original = linkage_covariance(file_path1, order)
	df_rlc_anon = linkage_covariance(file_path2, order)
	df_rlc_abs = df_rlc_anon.subtract(df_rlc_original, fill_value=0)
	#distance = spatial.distance.cdist(df_rlc_original, df3, 'euclidean')
	df_rlc_abs = df_rlc_abs.abs()
	sh_df = shannon(df_rlc_abs)
	#sh_df = sh_df.cumsum()
	plt.figure()
	sh_df.plot()
	# logging.info("Distance computed")
	# print distance
	# logging.info(distance)

def shannon(df):
	entropy = -(np.sum((df * np.log(df))/np.log(2.0), dtype=np.float64))
	return entropy

def compute_lcm(args):
	logging.info("Inside parallelization code")
	start = time.time()
	try:
		lc_matrix, i_gr, i_gr_r, dim = args
		for i in range(0, dim-1):
			for j in range(1, dim):
				#LC[i,j] <- sum(ax[i,]*ax[j,])/n - ax.r[i]*ax.r[j]/n^2 ; LC[j,i] <- LC[i,j]
				lc_matrix[i,j] = np.cov(i_gr[i,],i_gr[j,], ddof = 0)[0][1]
				lc_matrix[j,i] = lc_matrix[i,j]
		print("--- %s seconds ---" % (time.time() - start))
	except Exception as e:
		logging.error("Unexpected exception inside parallel processing: ", e)
	finally:
		logging.info("Done parallelizing. Returning the object now!")
		return lc_matrix

def main():
	print "rankedlinkcharacteristicvectorL2normed_el: Let us calculate the linkage covariance of the graph"
	original_graph_path = raw_input("Path of the graph file 1?\n")
	anon_graph_path = raw_input("Path of the graph file 2?\n")
	order = raw_input("Order of hops?\n")
	mp.freeze_support()
	start_time = time.time()
	h_dist(original_graph_path, anon_graph_path, order)
	print("--- %s seconds ---" % (time.time() - start_time))

main()
plt.show()

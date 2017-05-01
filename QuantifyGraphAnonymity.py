import numpy as np
import pandas as pd
import logging, sys

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

def linkage_covariance(file_path):
	#Delimiter
	tabDelimiter = "\t"

	logging.info('Got the file details. Reading file...')
	
	#Use the load text method of numpy with skiprows to skip the header(1st row)
	d = np.loadtxt(file_path, delimiter="\t", skiprows=1)

	logging.info('File read successfully')

	logging.debug('The input file has been processed for info!\n')
	logging.debug(d)

	#Use the shape command to get the dimensions of the matrix
	logging.info('Parsing the graph to get the number of nodes')
	[m, n] = d.shape
	header_row = list()
	logging.info('Creating a header row for the number of nodes')
	for i in range(0,n):
		header_row.append(i+1)
	
	logging.info('Graph parsed successfully and header row created!')
	
	# print the header row
	logging.debug('Header Row!\n')
	logging.debug(header_row)

	logging.info('Insert the header row above')
	#colnames(ax) <- c(1:63)
	input_graph = d
	# TO-DO ------ not working!!! - Start
	#np.insert(input_graph, 1, header_row, axis=0)
	# TO-DO ------ not working!!! - End

	logging.info('Header row inserted!!')
	
	logging.info('Finding the sum of all the values in each row')
	#ax.r <- rowSums(ax)
	input_graph_rowsums = input_graph.sum(axis=1)

	logging.info('Found the sum!')
	#LC <- matrix(0,n,n)
	# matrix of zeroes
	linkage_covariance_matrix = np.zeros((n, n))

	logging.info('Find the linkage covariance and compute the linkage covariance matrix')
	
	#for (i in 1:(n-1)) for (j in (i+1):n) 
	for i in range(0, n-1):
		for j in range((i+1), n):
			#LC[i,j] <- sum(ax[i,]*ax[j,])/n - ax.r[i]*ax.r[j]/n^2 ; LC[j,i] <- LC[i,j]
			linkage_covariance_matrix[i,j] = sum(input_graph[i,]*input_graph[j,])/n - input_graph_rowsums[i]*input_graph_rowsums[j]/n**2
			linkage_covariance_matrix[j,i] = linkage_covariance_matrix[i,j]

	logging.info('Linkage covariance matrix done')
	
	logging.debug('The linkage_covariance_matrix!\n')
	logging.debug(linkage_covariance_matrix)
	logging.debug('The input_graph!\n')
	logging.debug(input_graph)

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

	sort_key = list()
	ascending_key = list()
	for i in range(0,n):
		sort_key.append(i)
		ascending_key.append(False)

	df_rlc = pd.DataFrame(rlc)
	
	logging.debug('RLC Data Framed!\n')
	logging.debug(df_rlc)
	
	logging.info('Use a sort key and sort the matrix by each column in ascending order')
	
	df_rlc.sort_values(sort_key, ascending=ascending_key, inplace=True)
	rlc = df_rlc.as_matrix()
	
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
	
	return rlc

def main():
	print "rankedlinkcharacteristicvectorL2normed: Let us calculate the linkage covariance of the graph"
	path = raw_input("Path of the graph file?\n")
	print linkage_covariance(path)

main()
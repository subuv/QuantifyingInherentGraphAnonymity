library(sna)
library(network)
library(igraph)

g <- read.graph("/Users/sriram/hpdc_dk2/graphs/terrornet_edges", format = "edgelist", directed = TRUE)
ax.start <- get.adjacency(g)
ax <- ax.start
colnames(ax) <- c(1:63)
n <- nrow(ax)
ax.r <- rowSums(ax)

#create matrix of link covariance values
LC <- matrix(0,n,n)
for (i in 1:(n-1)) for (j in (i+1):n) 
{
  LC[i,j] <- sum(ax[i,]*ax[j,])/n ; LC[j,i] <- LC[i,j]
}

head(LC)
#create sorted LC matrix, first sort within each row
RLC.o <- matrix(0,n,n)
rownames(RLC.o) <- c(1:63)
for (i in 1:n) {RLC.o[i,] <- sort(LC[i,],decreasing=T)}
head(RLC.o)

#and the sort each row by col1 to coln in decreasing LC value
RLC.o <- RLC.o[do.call(order,-as.data.frame(RLC.o)),]
head(RLC.o)

#now normalize by the L2 modulus
for (i in 1:n)
{
  RLC.o[i,] <- RLC.o[i,]/sqrt(sum(RLC.o[i,]^2))
}
head(RLC.o)

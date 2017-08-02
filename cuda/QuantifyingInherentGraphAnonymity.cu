#include <stdio.h>
#include <igraph/igraph.h>
#include <string.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <sys/time.h>
#include <algorithm>
#include <stdlib.h>
#include <iterator>
#include <vector>
// #include "./cuPrintf.cu"
igraph_neimode_t OUTALL;
//CUDA ERROR
void checkCudaError(cudaError_t e, const char* in) {
    if (e != cudaSuccess) {
        printf("CUDA Error: %s, %s \n", in, cudaGetErrorString(e));
        exit(EXIT_FAILURE);
    }
}

__global__ void TEST(int n, int* d_x, int* d_y, int* dy_num, int* d_num) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;     
    
    if(i<n)
    {
        // d_x[i] = (int *)malloc(i * sizeof(int));
        // int *x;

        int l = 0;
        d_num[i] = 0;
        __syncthreads();
        for(int j = 0; j< i; j++)
        {
            if(i%150)
            {
                d_x[i*n + j] = i*j;
                d_y[i*n + j] = i + j;
                // d_num[i] = d_num[i] + 1;
                l++;
            }
        }
        __syncthreads();
        d_num[i] = l;
        // if(i==0)
        // {
        //  for(int j=1; j<n; j++)
        //  {
        //      d_num[j] = d_num[j] + d_num[j - 1];
        //      // printf("%d - %d\t", j, d_num[j]);
        //      __syncthreads();
        //  }
        //  // printf("%d - %d\t", i, d_num[i]);
        // }

    }

}


void TEST_PREP() {
    int n = 63;
    int *x, *d_x, *y, *d_y;
    x =(int *) malloc(n * n * sizeof(int));
    y =(int *) malloc(n * n * sizeof(int));

    int *num, *d_num, *numy, *dy_num;
    num =(int *) malloc(n * sizeof(int));
    numy =(int *) malloc(n * sizeof(int));

    checkCudaError(cudaMalloc((void**)&d_x, n * n *sizeof(int)), "Malloc Error d_x");
    checkCudaError(cudaMalloc((void**)&d_y, n * n *sizeof(int)), "Malloc Error d_x");

    checkCudaError(cudaMalloc((void**)&d_num, n* sizeof(int)), "Malloc Error");
    checkCudaError(cudaMalloc((void**)&dy_num, n* sizeof(int)), "Malloc Error");
    
    cudaMemcpy(d_x, x, n * n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, n * n * sizeof(int), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_num, num, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dy_num, numy, n * sizeof(int), cudaMemcpyHostToDevice);
    
    int numThreads = 32;

   dim3 DimGrid(ceil(n/numThreads), 1, 1);   
   if (n%numThreads) 
   {
        DimGrid.x++;   
   }
   dim3 DimBlock(numThreads, 1, 1);

    TEST<<<DimGrid,DimBlock>>>(n, d_x, d_y, dy_num,d_num);
    checkCudaError(cudaGetLastError(), "Checking Last Error, Kernel Launch");
    
    cudaMemcpy(num, d_num, n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(numy, dy_num, n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(x, d_x, n * n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(y, d_y, n * n * sizeof(int), cudaMemcpyDeviceToHost);
    
    for(int i = 0; i < n; i++)
    {
        printf("\n---%d\n", num[i]);
        for(int j=0;j<num[i];j++)
            printf("%d\t", x[i*n + j]);
            // int start = 0;
        // if(i>0)
        // {
        //  start = num[i-1];
        // }
        // printf("%d , %d - %d:\t", num[i], start, i);
        // for(int j=start;j<num[i] - 1;j++)
        // {
        //  printf("%d\t", x[j]);
        // }
        // printf("\n");
        
    }
    
    // return 0;
}


__host__ __device__  int compare (const void * a, const void * b)
{
// return (*(int*)a > *(int*)b) - (*(int*)b > *(int*)a);
  return ( *(int*)a - *(int*)b );
}

__host__ __device__  int commonNeighbor(int arr1[], int arr2[], int m, int n)
{
  int i = 0, j = 0;
  int numCommon = 0;
  while (i < m && j < n)
  {
    // printf("%d <--> %d\n", arr1[i], arr2[j]);
    if(arr1[i] < 0)
    {
        i++;
        continue;
    }
    if(arr2[j] < 0)
    {
        j++;
        continue;
    }
    // else
    {
        if (arr1[i] < arr2[j])
          i++;
        else if (arr2[j] < arr1[i])
          j++;
        else /* if arr1[i] == arr2[j] */
        {
          // printf(" %d ", arr2[j++]);
          j++;
          i++;
          numCommon++;
        }
    }
  }
  return numCommon;
}

__host__ __device__  __device__ int equalArray(int *a1, int *a2, int n_vertices)
{
    for(int i = 0; i < n_vertices; i++)
    {
        // printf("%d \? %d\t", a1[i], a2[i]);
        if(a1[i] != a2[i])
            return 0;
    }
    return 1;
}

/*
void LCM_Test(igraph_t &graph, igraph_neimode_t OUTALL, int numThreads)
{
    int n_vertices = igraph_vcount(&graph);
    igraph_adjlist_t al;
    igraph_adjlist_init(&graph, &al, OUTALL);
    igraph_adjlist_simplify(&al);

    long int **adjList;
    int *sizeAdj;

    adjList = (long int **) calloc(n_vertices, sizeof(long int *));
    sizeAdj = (int *) calloc(n_vertices, sizeof(int));
    for (int i = 0; i < n_vertices; i++) {
        igraph_vector_t *adjVec = igraph_adjlist_get(&al, i);

        adjList[i] = (long int *) calloc(igraph_vector_size(adjVec), sizeof(long int *));
        sizeAdj[i] = (int) igraph_vector_size(adjVec);
        for(int k = 0; k< igraph_vector_size(adjVec); k++)
        {
            adjList[i][k] = (long int) VECTOR(*adjVec)[k];
        }
    }

    for(int i = 0; i< n_vertices; i++)
    {
        qsort(adjList[i], sizeAdj[i], sizeof(long int), compare);
    }
}
*/
__global__ void Get_LCMSize_Kernel(int *d_adjList, int *d_sizeAdj, int *d_LCMSize, int n_vertices)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;     
    if(i<n_vertices)
    {
        int indexUsed = 0;
        int iStart = 0, iEnd = 0;
        int k = 0;

        if(i > 0)
        {       
            k = d_sizeAdj[i-1];
        }

        iEnd = d_sizeAdj[i];

        __syncthreads();

        for(int j = 0; j < n_vertices; j++) {
            if(i==j)
                continue;
            iStart = k;
            int jStart = 0, jEnd = 0;

            if(j > 0)
                jStart = d_sizeAdj[j-1];
            jEnd = d_sizeAdj[j];
            
            int compVec = 0;

            while (iStart < iEnd && jStart < jEnd)
            {
                    if(d_adjList[iStart] < d_adjList[jStart])
                        iStart++;
                    else if (d_adjList[jStart] < d_adjList[iStart])
                        jStart++;
                    else // if arr1[i] == arr2[j] 
                    {
                        jStart++;
                        iStart++;
                        compVec++;
                        break;
                    }
            }

            if (compVec > 0)
            {
                indexUsed++;
            }
        }
    
        __syncthreads();
        d_LCMSize[i] = indexUsed;
        // __syncthreads();
    
    }

}

__global__ void Get_LCM_Kernel(int *d_adjList, int *d_sizeAdj, int *d_lcmMatrix, int *d_LCMSize, int n_vertices)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;     
    if(i<n_vertices)
    {
        int indexUsed = 0, indexOffset = 0;
        int iStart = 0, iEnd = 0;
        int k = 0;

        if(i > 0)
        {       
            k = d_sizeAdj[i-1];
            indexOffset = d_LCMSize[i-1];
        }

        iEnd = d_sizeAdj[i];
        
        for(int j = indexOffset; j<iEnd; j++)
        {
            d_lcmMatrix[j] = 0;
        }

        __syncthreads();

        for(int j = 0; j < n_vertices; j++) {
            if(i==j)
                continue;
            iStart = k;
            int jStart = 0, jEnd = 0;

            if(j > 0)
                jStart = d_sizeAdj[j-1];
            jEnd = d_sizeAdj[j];
            
            int compVec = 0;

            while (iStart < iEnd && jStart < jEnd)
            {
                    if(d_adjList[iStart] < d_adjList[jStart])
                        iStart++;
                    else if (d_adjList[jStart] < d_adjList[iStart])
                        jStart++;
                    else // if arr1[i] == arr2[j] 
                    {
                        jStart++;
                        iStart++;
                        compVec++;
                    }
            }

            if (compVec > 0)
            {
                atomicAdd((int*)&d_lcmMatrix[indexUsed + indexOffset], compVec);
                // d_lcmMatrix[indexUsed + indexOffset] = compVec;
                indexUsed++;
            }
            // __syncthreads();
        }
    
        // __syncthreads();
        // d_LCMSize[i] = indexUsed;
        // __syncthreads();
    
    }

}

__global__ void LCM_Hist_Kernel(int *d_lcmMatrix, int *d_LCMSize, int *d_histogram, int n_vertices)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int count = 0, countMax = -1;
  
    if(i<n_vertices)
    {
        int iStart = 0;
        if(i>0)
            iStart = d_LCMSize[i - 1]; //Offset
        count = 0;
        int iSize = d_LCMSize[i] - iStart;

        for(int j = 0; j < n_vertices; j++) {
            int jStart = 0;
            if(j>0)
                jStart = d_LCMSize[j - 1]; //Offset

            int jSize = d_LCMSize[j] - jStart;
            if(iSize != jSize)
                continue;
      
            int eq = 1;
            for(int k = 0; k < iSize; k++)
            {
                if(d_lcmMatrix[iStart + k] != d_lcmMatrix[jStart + k])
                {
                    eq = 0;
                    break;
                }
            }
            if(eq == 1)
            {               
                count++;
            }
        }

        if(countMax < count)
            countMax = count;
        atomicAdd((int*)&d_histogram[count], 1);
        // d_histogram[count]++;
    }
}

int LCM_Kernel_Prep(igraph_t &graph, igraph_neimode_t OUTALL, int numThreads)
{
    printf("Allocating Adjacency List\n");
    int n_vertices = igraph_vcount(&graph);
    igraph_adjlist_t al;
    igraph_adjlist_init(&graph, &al, OUTALL);
    igraph_adjlist_simplify(&al);

    int **adjList2D;
    int totalSize = 0;

    int *adjList, *d_adjList;
    int *sizeAdj, *d_sizeAdj;

    int *lcmMatrix, *d_lcmMatrix;

    int *d_LCMSize, *LCMSize, *LCMSize_Calc;
    
    adjList2D = (int **) calloc(n_vertices, sizeof(int *));
    sizeAdj = (int *) malloc(n_vertices * sizeof(int));
    LCMSize = (int *) malloc(n_vertices * sizeof(int));
    LCMSize_Calc = (int *) malloc(n_vertices * sizeof(int));
    memset(LCMSize, 0, n_vertices*sizeof(int));
    memset(LCMSize_Calc, 0, n_vertices*sizeof(int));
    printf("Computing Adjacency List - %d vertices...\n", n_vertices);

    for (int i = 0; i < n_vertices; i++) {
        igraph_vector_t *adjVec = igraph_adjlist_get(&al, i);

        // igraph_vector_t adjVec;
        // igraph_vector_init(&adjVec, 0);
        // igraph_neighbors(&graph, &adjVec, i, OUTALL);

        adjList2D[i] = (int *) malloc(igraph_vector_size(adjVec) * sizeof(int));
        sizeAdj[i] = (int) igraph_vector_size(adjVec);
        totalSize += sizeAdj[i];
        for(int k = 0; k< igraph_vector_size(adjVec); k++)
        {
            adjList2D[i][k] = (int) VECTOR(*adjVec)[k];
        }
    }

    for(int i = 0; i< n_vertices; i++)
    {
        qsort(adjList2D[i], sizeAdj[i], sizeof(int), compare);
    }
    
    adjList = (int *) malloc(totalSize * sizeof(int));
    int l = -1;
    for (int q = 0; q < n_vertices; q++)
    {
        for (int t = 0; t < sizeAdj[q]; t++)
        {
            l++;
            adjList[l] = adjList2D[q][t];
        }
    }
    for(int i = 0; i< n_vertices; i++)
    {
        free(adjList2D[i]);
        if(i>0)
        {
            sizeAdj[i] += sizeAdj[i - 1];
        }
    }
    
    free(adjList2D);
    // memset(LCMSize, 0, n_vertices*sizeof(int));
    printf("%d-%d\n", totalSize, sizeAdj[n_vertices-1]);
    printf("Got Adj List...\n Allocating on gpu mem...");
    checkCudaError(cudaMalloc((void**)&d_adjList, totalSize * sizeof(int)), "Malloc Error d_adjList");
    checkCudaError(cudaMalloc((void**)&d_sizeAdj, n_vertices * sizeof(int)), "Malloc Error d_sizeAdj");
    checkCudaError(cudaMalloc((void**)&d_LCMSize, n_vertices * sizeof(int)), "Malloc Error d_sizeAdj");

    cudaMemcpy(d_adjList, adjList, totalSize * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sizeAdj, sizeAdj, n_vertices * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_LCMSize, LCMSize_Calc, n_vertices * sizeof(int), cudaMemcpyHostToDevice);

    dim3 DimGrid(ceil(n_vertices/numThreads), 1, 1);   
    if (n_vertices%numThreads) 
    {
        DimGrid.x++;
    }

    dim3 DimBlock(numThreads, 1, 1);
    int totLCMSize = 0;
    printf("Launching Size Kernel...\n");
    Get_LCMSize_Kernel<<<DimGrid,DimBlock>>>(d_adjList, d_sizeAdj, d_LCMSize, n_vertices);
    cudaThreadSynchronize();
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(), "Checking Last Error, Size Kernel Launch");
    cudaMemcpy(LCMSize_Calc, d_LCMSize, n_vertices * sizeof(int), cudaMemcpyDeviceToHost);
    
    for(int i = 0; i<n_vertices; i++)
    {
        totLCMSize += LCMSize_Calc[i];
        LCMSize[i] = LCMSize_Calc[i];
    }

    for(int i = 1; i<n_vertices; i++)
    {
        if(i>0)
            LCMSize[i] += LCMSize[i - 1];
    }
    printf("%d - %d\n", totalSize, totLCMSize);
    
    lcmMatrix = (int *) malloc(totLCMSize * sizeof(int));
    memset(lcmMatrix, 0, totLCMSize*sizeof(int));
    checkCudaError(cudaMalloc((void**)&d_lcmMatrix, totLCMSize * sizeof(int)), "Malloc Error d_lcmMatrix");
    cudaMemcpy(d_lcmMatrix, lcmMatrix, totLCMSize * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_LCMSize, LCMSize, n_vertices * sizeof(int), cudaMemcpyHostToDevice);
    printf("Launching LCM Kernel...\n");
    
    // LCM_Kernel<<<DimGrid,DimBlock>>>(d_adjList, d_sizeAdj, d_lcmMatrix, d_LCMSize, n_vertices);
    Get_LCM_Kernel<<<DimGrid,DimBlock>>>(d_adjList, d_sizeAdj, d_lcmMatrix, d_LCMSize, n_vertices);
    
	cudaThreadSynchronize();
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(), "Checking Last Error, Kernel Launch");
    printf("Copying to CPU Memory...\n");
    checkCudaError(cudaMemcpy(lcmMatrix, d_lcmMatrix, totLCMSize * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy Error d_lcmMatrix");
    // checkCudaError(cudaMemcpy(LCMSize, d_LCMSize, n_vertices * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy Error LCMSize");
    
    // cudaFree(d_lcmMatrix);
    // cudaFree(d_LCMSize);
    cudaFree(d_adjList);
    cudaFree(d_sizeAdj);
    free(sizeAdj);
    free(adjList);

    printf("Allocating Histogram...\n");
    int *histo, *d_histogram;
    histo = (int *) malloc(n_vertices * sizeof(int));
    memset(histo, 0, sizeof(int)*n_vertices);
    checkCudaError(cudaMalloc((void**)&d_histogram, n_vertices * sizeof(int)), "Malloc Error d_histogram");
    cudaMemcpy(d_histogram, histo, n_vertices * sizeof(int), cudaMemcpyHostToDevice);

    printf("Sorting LCM...\n");
    
	for(int i = 0; i< n_vertices; i++)
	{
		int offset = 0;
		if(i > 0)
		{
			offset = LCMSize[i - 1];
		}
		qsort(lcmMatrix + offset, LCMSize_Calc[i], sizeof(int), compare);
	}

    cudaMemcpy(d_lcmMatrix, lcmMatrix, totLCMSize * sizeof(int), cudaMemcpyHostToDevice);

	printf("Computing Histogram...\n");
	// return 0;
    printf("Launching Histogram Kernel...\n");
    
    LCM_Hist_Kernel<<<DimGrid,DimBlock>>>(d_lcmMatrix, d_LCMSize, d_histogram, n_vertices);
    cudaThreadSynchronize();
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(), "Checking Last Error, Kernel Launch");
    printf("Copying to CPU Memory...\n");
    checkCudaError(cudaMemcpy(histo, d_histogram, n_vertices * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy Error d_lcmMatrix");

    
    printf("Finished Histogram...\n");
    for(int i = 1; i < n_vertices; i++) {
        if ((long) (histo[i] / i) > 0)
            printf("%d    %d\n", i, (int) (histo[i] / i));
    }
    free(lcmMatrix);
    free(LCMSize_Calc);
    free(LCMSize);
    free(histo);
    cudaFree(d_histogram);
    cudaFree(d_lcmMatrix);
    cudaFree(d_LCMSize);
	return 0;
}

int main(int argc, char** argv)
{

    // TEST_PREP();
    int numThreads = 16;
    //checks arguments
    if(argc < 3) {

        printf("\nToo few arguments. Usage: ./%s graphFile all/out\n", argv[0]);
        return -1;
    }

    //graph direction out or all
    if(!strcmp(argv[2], "out"))
        OUTALL = IGRAPH_OUT;
    else
        OUTALL = IGRAPH_ALL;
    
    struct timeval stop, start;
    gettimeofday(&start, NULL);

    //opens graph file passed as 1st argument
    FILE *inputFile;
    inputFile = fopen(argv[1], "r");
    if(inputFile == NULL)
    {
        printf("Could not load input file...\n");
        return 1;
    }
    
    igraph_t graph;

    //builds graph from file
    igraph_read_graph_ncol(&graph, inputFile, NULL, true, IGRAPH_ADD_WEIGHTS_NO, IGRAPH_DIRECTED);

    // TEST_PREP();
    
    LCM_Kernel_Prep(graph, OUTALL, numThreads);
    
    // LCM_Test(graph, OUTALL, numThreads);

    gettimeofday(&stop, NULL);
    printf("took %2f\n", (stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f);
}

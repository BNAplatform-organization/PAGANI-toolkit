#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <cuda_runtime_api.h>

using namespace std;
#include "BC_GPU.h"
void CUBC(float *BC, int * row, int * col, int numVertices, int numEdges)
{
	int count[1];
	cudaError_t error;
	error = cudaGetDeviceCount(count); 
	if (error != cudaSuccess)
	{
		cerr<<"no CUDA device found."<<endl;
		return;
	}
	int device[10];
	cudaGetDevice(device); 
	cudaDeviceProp prop[1];
	int best, bestCount = 0;
	for (int i = 0; i < count[0]; i++)
	{
		cudaGetDeviceProperties(prop, device[i]);
		if (prop-> multiProcessorCount > bestCount)
		{
			bestCount = prop-> multiProcessorCount;
			best = i;
		}
	}
	cudaSetDevice(device[best]);
	cudaGetDeviceProperties(prop, device[best]);
	printf("GPU Device %d: \"%s\" with compute capability %d.%d\n", device[best], prop[0].name, prop[0].major, prop[0].minor);

	size_t GM_size = prop->totalGlobalMem;
	int maxgrid;
	//cout<<GM_size;
	//cout<<endl<<sizeof(bool);
	
	maxgrid = GM_size / (2*numEdges + 6 * numVertices) / sizeof(int);
	
	//cout <<endl<< maxgrid<<endl;
    #if (ARRAY)
	maxgrid = (GM_size - 2 * numEdges * sizeof(int)) / (6 * numVertices * sizeof(int) + numEdges)/1.5 ;
	#endif
	cout <<endl<< maxgrid<<endl;

	int grid = min(128, maxgrid);
	if (grid > 8)
		grid -= grid%8;
	else if (grid > 4)
		grid = 4;
	else if (grid > 2)
		grid = 2;
	else 
		grid = 1;
    #if (ARRAY)
	grid = max(grid, 12);
	#endif

	int thread = 512;
	cout<<"sizeGrid = "<<grid<<endl;

	// computing Betweenness
	cout<<"Start BC computing on GPU. n = "<<numVertices<<", m = "<<numEdges<<endl;
	clock_t gpu_time = clock();

	#if (NODE && ARRAY)
	Betweenness_GPU_node_array(row, col, numVertices, numEdges, BC, grid, thread);
	#endif

	#if (NODE && PRED)
	Betweenness_GPU_node_pred(row, col, numVertices, numEdges, BC, grid, thread);
	#endif

	#if (NODE && SUCC)
	Betweenness_GPU_node_succ(row, col, numVertices, numEdges, BC, grid, thread);
	#endif

	#if (EDGE && ARRAY)
	int * row_full = new int [numEdges];
	for (int i = 0; i < numVertices; i++)
	{
		for (int j = row[i]; j < row[i + 1]; j++)
		{
			row_full[j] = i;
		}
	}
	Betweenness_GPU_edge_array(row, row_full, col, numVertices, numEdges, BC, grid, thread);
	delete []row_full;
	#endif

	#if (EDGE && PRED)
	int * row_full = new int [numEdges];
	for (int i = 0; i < numVertices; i++)
	{
		for (int j = row[i]; j < row[i + 1]; j++)
		{
			row_full[j] = i;
		}
	}
	Betweenness_GPU_edge_pred(row, row_full, col, numVertices, numEdges, BC, grid, thread);
	delete []row_full;
	#endif

	#if (EDGE && SUCC)
	int * row_full = new int [numEdges];
	for (int i = 0; i < numVertices; i++)
	{
		for (int j = row[i]; j < row[i + 1]; j++)
		{
			row_full[j] = i;
		}
	}
	Betweenness_GPU_edge_succ(row, row_full, col, numVertices, numEdges, BC, grid, thread);
	delete []row_full;
	#endif

	gpu_time = clock() - gpu_time;
	cout<<gpu_time<<" ms for GPU BC"<<endl;
	cout<<(double)numEdges*numVertices/gpu_time*1000/1024/1024<<" MTEPS"<<endl;
}
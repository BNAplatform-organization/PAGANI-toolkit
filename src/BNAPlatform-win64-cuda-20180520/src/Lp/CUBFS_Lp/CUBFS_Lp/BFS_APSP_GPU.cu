#include <stdio.h>
#include <cuda_runtime.h>
#include "device_functions.h"
#include "BFS_APSP_GPU.cuh"
#define VIRTUAL_WARP 8

__global__ void APSP_BFS_node_kernel(int *r, int *c, int * dist, int numVertices, int numEdges, int offset_source) 
{
int offset_vertices = blockIdx.x * numVertices;
	int offset_edge = blockIdx.x * numEdges;

	for (int i = threadIdx.x; i < numVertices; i += blockDim.x)
	{
		dist[offset_vertices + i] = -1;
	}

	int edge_index = threadIdx.x % VIRTUAL_WARP;
	int vertice_index = threadIdx.x / VIRTUAL_WARP;
	int source = blockIdx.x + offset_source;
	if (source >= numVertices)
		return;
	__shared__ bool done;
	done = false;
	
	int level = 0;
	dist[offset_vertices + source] = level++;
	while (!done)
	{
		__syncthreads(); // attention: this sync is neccessary
		done = true;
		for (int current = vertice_index; current < numVertices; current += blockDim.x / VIRTUAL_WARP)
		{
			if (dist[offset_vertices + current] != level - 1)
				continue;
			for (int j = r[current] + edge_index; j < r[current + 1]; j += VIRTUAL_WARP)
			{
				int next = c[j];
				int read_dist = dist[offset_vertices + next];
				if (read_dist == -1)
				{
					dist[offset_vertices + next] = level;
					done = false;
				}
			}
		}
		level ++;
			__syncthreads();  	
	}							  	
}

void APSP_GPU(int * dist, int *r, int *c, int numVertices, int numEdges, int grid, int thread)
{
    int devID;
    cudaDeviceProp deviceProps;

    devID = findCudaDevice();

    // get number of SMs on this GPU
    checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
    //printf("CUDA device [%s] has %d Multi-Processors\n", deviceProps.name, deviceProps.multiProcessorCount);

	//int thread = 256;
	//int grid = 100;
    // allocate device memory
    int* d_r; 
	int* d_c;
	int* d_dist;

    checkCudaErrors( cudaMalloc( (void**) &d_r, sizeof(int) * (numVertices + 1)));
	checkCudaErrors( cudaMalloc( (void**) &d_c, sizeof(int) * numEdges));

	
    // copy host memory to device
	checkCudaErrors( cudaMemcpy( d_r, r, sizeof(int) * (numVertices + 1), cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy( d_c, c, sizeof(int) * numEdges, cudaMemcpyHostToDevice) );

    // allocate device memory for result
	checkCudaErrors( cudaMalloc( (void**) &d_dist, sizeof(int) * numVertices * grid));

	clock_t kernel_time = 0;
	clock_t transfer_time = 0;

    // execute the kernel
	for (int offset_source = 0; offset_source < numVertices; offset_source += grid)
	{
		clock_t time = clock();
		APSP_BFS_node_kernel<<<grid, thread>>>(d_r, d_c, d_dist, numVertices, numEdges, offset_source);
		// check if kernel execution generated and error
		getLastCudaError("Kernel execution failed");
		cudaThreadSynchronize();
		time = clock() - time;

		cout<<offset_source<<" done. Time = "<<time<<"ms."<<endl;
		kernel_time += time;

		time = clock();
		// copy result from device to host
		if(numVertices - offset_source > grid)			 
			checkCudaErrors(cudaMemcpy(dist + (long long)offset_source * numVertices, d_dist, sizeof(float) * numVertices * grid, cudaMemcpyDeviceToHost));
		else
			checkCudaErrors(cudaMemcpy(dist + (long long) offset_source * numVertices, d_dist, sizeof(float) * numVertices * (numVertices%grid), cudaMemcpyDeviceToHost));
		time = clock() - time;
		transfer_time += time;
	}

	cout<<"total kernel time: "<<kernel_time<<"ms."<<endl;
	cout<<"total transfering time: "<<transfer_time<<"ms."<<endl;

    // cleanup memory
    checkCudaErrors(cudaFree(d_r));
    checkCudaErrors(cudaFree(d_c));
	checkCudaErrors(cudaFree(d_dist));
    cudaDeviceReset();
}
#include "BC_GPU.cuh"

using namespace std;

__global__ void
betweenness_edge_array_kernel(float *BC, int * r, int * edge_begin, int *edge_end, int * dist, float * sigma, float * delta, bool * P, int numVertices, int numEdges, int offset_source) 
{
	int offset_vertices = blockIdx.x * numVertices;
	int offset_edge = blockIdx.x * numEdges;

	for (int i = threadIdx.x; i < numVertices; i += blockDim.x)
	{
		dist[offset_vertices + i] = -1;
		sigma[offset_vertices + i] = 0;
		delta[offset_vertices + i] = 0;
	}
	for (int i = threadIdx.x; i < numEdges; i += blockDim.x)
	{
		P[offset_edge + i] = false;
	}

	int source = blockIdx.x + offset_source;
	if (source >= numVertices)
		return;
	__shared__ bool done;
	done = false;
	
	int level = 0;
	dist[offset_vertices + source] = level++;
	sigma[offset_vertices + source] = 1; 
	while (!done)
	{
		__syncthreads(); // attention: this sync is neccessary
		done = true;
		for (int edge = threadIdx.x; edge < numEdges; edge += blockDim.x)
		{
			int current = edge_begin[edge];
			if (dist[offset_vertices + current] != level - 1)
				continue;

			int next = edge_end[edge];
			int read_dist = dist[offset_vertices + next];
			if (read_dist == -1)
			{
				dist[offset_vertices + next] = level;
				done = false;
			}

			if (read_dist < level && read_dist >= 0)
				continue;

			atomicAdd(sigma + offset_vertices + next, sigma[offset_vertices + current]); //atomic!
			P[offset_edge + edge] = true;
			
		}
		level ++;
		__syncthreads();
	}
	for (int i = level - 2; i >= 0; i--)
	{
		//// slower method due to atomic operations and bandwidth wastes. 
		//for (int edge = threadIdx.x; edge < numEdges; edge += blockDim.x)
		//{
		//	int current = edge_begin[edge];
		//	if (dist[offset_vertices + current] != i)
		//		continue;
		//	
		//	if (P[offset_edge + edge] != 1)
		//		continue;
		//	int next = edge_end[edge];
		//	atomicAdd(delta + offset_vertices + current, (double) sigma[offset_vertices + current] / sigma[offset_vertices + next]*(1 + delta[offset_vertices + next]));

		//}
		//__syncthreads();

		for (int current = threadIdx.x; current < numVertices; current += blockDim.x)
		{
			if (dist[offset_vertices + current] != i)
				continue;
			for (int j = r[current]; j < r[current + 1]; j += 1)
			{
				if (P[offset_edge + j] != true)
					continue;
				int next = edge_end[j];
				delta[offset_vertices + current] += (double) sigma[offset_vertices + current] / sigma[offset_vertices + next]*(1 + delta[offset_vertices + next]);
			}
		}
		__syncthreads();
	}

	for (int current = threadIdx.x; current < numVertices; current += blockDim.x)
	{
		if(current != source)
			atomicAdd(BC + current, delta[offset_vertices + current]);
	}
}

void Betweenness_GPU_edge_array(int * r, int *r_full, int *c, int numVertices, int numEdges, float *BC, int grid, int thread)
{
    int devID;
    cudaDeviceProp deviceProps;

    devID = findCudaDevice();

    // get number of SMs on this GPU
    checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
    printf("CUDA device [%s] has %d Multi-Processors\n", deviceProps.name, deviceProps.multiProcessorCount);

	//int thread = 256;
	//int grid = 100;
    // allocate device memory
	int* d_r;
	int* d_c;
	int* d_r_full;
	int* dist;
	float* sigma;
	float* delta;
	bool* P;

	checkCudaErrors( cudaMalloc( (void**) &d_r, sizeof(int) * numVertices));
	checkCudaErrors( cudaMalloc( (void**) &d_r_full, sizeof(int) * numEdges));
	checkCudaErrors( cudaMalloc( (void**) &d_c, sizeof(int) * numEdges));

	checkCudaErrors( cudaMalloc( (void**) &dist, sizeof(int) * numVertices * grid));
	checkCudaErrors( cudaMalloc( (void**) &sigma, sizeof(int) * numVertices * grid));
	checkCudaErrors( cudaMalloc( (void**) &delta, sizeof(int) * numVertices * grid));
	checkCudaErrors( cudaMalloc( (void**) &P, sizeof(bool) * numEdges * grid));
    // copy host memory to device
	checkCudaErrors( cudaMemcpy( d_r, r, sizeof(int) * numVertices, cudaMemcpyHostToDevice) );
    checkCudaErrors( cudaMemcpy( d_c, c, sizeof(int) * numEdges, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy( d_r_full, r_full, sizeof(int) * numEdges, cudaMemcpyHostToDevice) );

    // allocate device memory for result
    float* d_BC;
    checkCudaErrors( cudaMalloc( (void**) &d_BC, sizeof(float) * numVertices));
	checkCudaErrors( cudaMemset( d_BC, 0, sizeof(float) * numVertices));

	clock_t kernel_time = 0;

    // execute the kernel
	for (int offset_source = 0; offset_source < numVertices; offset_source += grid)
	{
		clock_t time = clock();
		betweenness_edge_array_kernel<<<grid, thread>>>(d_BC, d_r, d_r_full, d_c, dist, sigma, delta, P, numVertices, numEdges, offset_source);
		// check if kernel execution generated and error
		getLastCudaError("Kernel execution failed");
		cudaThreadSynchronize();
		time = clock() - time;

		cout<<offset_source<<" done. Time = "<<time<<"ms."<<endl;
		kernel_time += time;
	}

	cout<<"total kernel time: "<<kernel_time<<"ms."<<endl;

	// copy result from device to host
	checkCudaErrors(cudaMemcpy(BC, d_BC, sizeof(float) * numVertices, cudaMemcpyDeviceToHost));

    // cleanup memory
	checkCudaErrors(cudaFree(d_r));
    checkCudaErrors(cudaFree(d_r_full));
    checkCudaErrors(cudaFree(d_c));
	checkCudaErrors(cudaFree(d_BC));
	checkCudaErrors(cudaFree(dist));
    checkCudaErrors(cudaFree(sigma));
	checkCudaErrors(cudaFree(delta));
	checkCudaErrors(cudaFree(P));
    cudaDeviceReset();
}
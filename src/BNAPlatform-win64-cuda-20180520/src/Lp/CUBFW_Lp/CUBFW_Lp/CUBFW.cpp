# include <ctime>
# include <iostream>
using namespace std;

# include "CUBFW_ph1.h"
# include "CUBFW_ph2.h"
# include "CUBFW_ph3.h"
# include "BFW.h"
# include "mytimer.h"
# include "cuTranspose.h"
# include "cuda_runtime.h"
float * kernel_primary_block;
float * kernel_dst_ij;
float * kernel_src_ik;
float * kernel_src_kj;

void myprint_blocked(float * costmat, int numVertices, int block_size);

void cuAPSP(float *costmat, const int numVertices, const int block_size)
{
	int block_cnt = numVertices / block_size;
	cudaMalloc((void **)&kernel_primary_block, sizeof(float) * block_size * block_size);
	cudaMalloc((void **)&kernel_dst_ij, sizeof(float) * block_size * block_size);
	cudaMalloc((void **)&kernel_src_ik, sizeof(float) * block_size * block_size);
	cudaMalloc((void **)&kernel_src_kj, sizeof(float) * block_size * block_size);

	clock_t time_gpu = clock();

	for (int k_block = 0; k_block < block_cnt; k_block++)
	{
		cout<<"Round "<<k_block<<":"<<endl;
		// phase 1 gpu
		int block_row = k_block;
		int block_col = k_block;
		long long offset_ph1 = ((long long)block_row * block_cnt + block_col) * block_size * block_size;

		clock_t time_gpu_ph1 = clock();
		BFW_CUDA_small_Katz(costmat + offset_ph1, kernel_primary_block, block_size);
		time_gpu_ph1 = clock() - time_gpu_ph1;
		cout<<"Phase 1 time: "<<time_gpu_ph1<<" ms."<<endl;

		////phase 1 cpu
		//clock_t time_ph1_cpu = clock();
		//int block_row = k_block;
		//int block_col = k_block;
		//int offset_ph1 = (block_row * block_cnt + block_col) * block_size * block_size;

		//BFW_one_block_C(costmat + offset_ph1, costmat + offset_ph1, costmat + offset_ph1, block_size);
		//time_ph1_cpu = clock() - time_ph1_cpu;
		//cout<<"Phase1 cpu time: "<<time_ph1_cpu<<" ms."<<endl;

		//phase 2a gpu
		clock_t time_gpu_ph2 = clock();
		cudaMemcpy(kernel_primary_block, costmat + offset_ph1, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);
		for (int colrow = 0; colrow < block_cnt; colrow ++) 
		{
			int block_row = colrow;
			int block_col = k_block;
			long long offset_ph2 = ((long long)block_row * block_cnt + block_col) * block_size * block_size;
			cudaMemcpy(kernel_src_ik, costmat + offset_ph2, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);					
			BFW_CUDA_ph2(kernel_dst_ij, kernel_src_ik, kernel_primary_block, block_size);						
			cudaMemcpy(costmat + offset_ph2, kernel_dst_ij, sizeof(float) * block_size * block_size, cudaMemcpyDeviceToHost);
		}
		
		//phase 2b
		for (int colrow = 0; colrow < block_cnt; colrow ++) 
		{
			int block_row = k_block;
			int block_col = colrow;
			long long offset_ph2 = ((long long)block_row * block_cnt + block_col) * block_size * block_size;
			cudaMemcpy(kernel_src_kj, costmat + offset_ph2, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);
			BFW_CUDA_ph2(kernel_dst_ij, kernel_primary_block, kernel_src_kj, block_size);
			cudaMemcpy(costmat + offset_ph2, kernel_dst_ij, sizeof(float) * block_size * block_size, cudaMemcpyDeviceToHost);
		}
		time_gpu_ph2 = clock() - time_gpu_ph2;
		cout<<"Phase 2 time: "<<time_gpu_ph2<<" ms."<<endl;


		//// phase2 cpu
		//clock_t time_ph2_cpu = clock();
		//for (int colrow = 0; colrow < block_cnt; colrow ++) 
		//{
		//	if (colrow == k_block)
		//		continue;
		//	int block_row = k_block;
		//	int block_col = colrow;
		//	int offset_ph2 = (block_row * block_cnt + block_col) * block_size * block_size;
		//	BFW_one_block_C(costmat + offset_ph2, costmat + offset_ph1, costmat + offset_ph2, block_size);
		//}
		//for (int colrow = 0; colrow < block_cnt; colrow ++) 
		//{
		//	if (colrow == k_block)
		//		continue;
		//	int block_row = colrow;
		//	int block_col = k_block;
		//	int offset_ph2 = (block_row * block_cnt + block_col) * block_size * block_size;
		//	BFW_one_block_C(costmat + offset_ph2, costmat + offset_ph2, costmat + offset_ph1, block_size);
		//}
		//time_ph2_cpu = clock() - time_ph2_cpu;
		//cout<<"Phase2 cpu time: "<<time_ph2_cpu<<" ms."<<endl;
		
		//phase 3 gpu
		clock_t time_gpu_ph3 = clock();
		for (int block_row = 0; block_row < block_cnt; block_row ++) 
		{
			if (block_row == k_block) continue;
			int src_ik_row = block_row;
			int src_ik_col = k_block;
			long long offset_src_ik = ((long long)src_ik_row * block_cnt + src_ik_col) * block_size * block_size;
			cudaMemcpy(kernel_src_ik, costmat + offset_src_ik, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);	

			for (block_col = 0; block_col < block_cnt; block_col ++) 
			{
				if ( block_col == k_block ) continue;				
				int src_kj_row = k_block;
				int src_kj_col = block_col;
				long long offset_ph3 = ((long long)block_row * block_cnt + block_col) * block_size * block_size;			
				long long offset_src_kj = ((long long)src_kj_row * block_cnt + src_kj_col) * block_size * block_size;		
				
				cudaMemcpy(kernel_dst_ij, costmat + offset_ph3, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);					
				cudaMemcpy(kernel_src_kj, costmat + offset_src_kj, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);				
				BFW_CUDA_ph3(kernel_dst_ij, kernel_src_ik, kernel_src_kj, block_size);							
				cudaMemcpy(costmat + offset_ph3, kernel_dst_ij, sizeof(float) * block_size * block_size, cudaMemcpyDeviceToHost);
			}
		}
		time_gpu_ph3 = clock() - time_gpu_ph3;
		cout<<"Phase 3 time: "<<time_gpu_ph3<<" ms."<<endl;

		//// phase 3 cpu
		//clock_t time_cpu_ph3 = clock();
		//for (int block_row = 0; block_row < block_cnt; block_row ++) 
		//{
		//	if (block_row == k_block) continue;
		//	int src_ik_row = block_row;
		//	int src_ik_col = k_block;
		//	long long offset_src_ik = ((long long)src_ik_row * block_cnt + src_ik_col) * block_size * block_size;
		//	//cudaMemcpy(kernel_src_ik, costmat + offset_src_ik, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);	

		//	for (block_col = 0; block_col < block_cnt; block_col ++) 
		//	{
		//		if ( block_col == k_block ) continue;				
		//		int src_kj_row = k_block;
		//		int src_kj_col = block_col;
		//		long long offset_ph3 = ((long long)block_row * block_cnt + block_col) * block_size * block_size;			
		//		long long offset_src_kj = ((long long)src_kj_row * block_cnt + src_kj_col) * block_size * block_size;		
		//		BFW_one_block_C(costmat + offset_ph3, costmat + offset_src_ik, costmat + offset_src_kj, block_size);
		//		//cudaMemcpy(kernel_dst_ij, costmat + offset_ph3, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);					
		//		//cudaMemcpy(kernel_src_kj, costmat + offset_src_kj, sizeof(float) * block_size * block_size, cudaMemcpyHostToDevice);				
		//		//BFW_CUDA_ph3(kernel_dst_ij, kernel_src_ik, kernel_src_kj, block_size);							
		//		//cudaMemcpy(costmat + offset_ph3, kernel_dst_ij, sizeof(float) * block_size * block_size, cudaMemcpyDeviceToHost);
		//	}
		//}
		//time_cpu_ph3 = clock() - time_cpu_ph3;
		//cout<<"Phase 3 time: "<<time_cpu_ph3<<" ms."<<endl;

	}
	cudaThreadSynchronize();
	time_gpu = clock() - time_gpu;
	cout<<"¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï¡ï"<<endl;
	double n1 = numVertices / 1024.;
	double Gflop = 2 * n1 * n1 * n1;
	cout<<"Total GPU time (kernel + transfer): "<<time_gpu<<" ms. "<<endl;
	cudaFree(kernel_primary_block);
	cudaFree(kernel_dst_ij);
	cudaFree(kernel_src_ik);
	cudaFree(kernel_src_kj);

}

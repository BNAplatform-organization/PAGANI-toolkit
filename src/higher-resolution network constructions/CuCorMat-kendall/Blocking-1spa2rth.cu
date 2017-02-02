#include "cublas_v2.h"
#include "cuda_runtime.h"
#include "memory.h"
#include <iostream>
#include <ctime>
# include <fstream>
#include <vector>
#include <Windows.h>
#include<iomanip>
#include <thrust/device_vector.h>
#include"histogram.h"
#include <thrust/binary_search.h>
#include <sm_20_atomic_functions.h>
#include "help_func.cuh"
#include "pre_process.cuh"

//#define width  0.1//0.0001//best to be the multiples
#define width  0.000001//0.0001//best to be the multiples

using namespace std;
#pragma comment(lib,"cublas.lib")

extern __global__ void standardKernel(real__t* devCormat, int Batch_size, bool diagnoal);

void updataPointer(real__t** pointer, real__t* devBOLD, int L, int Batch_size, const int Transferblocknum)
{
	*pointer =  *pointer ==  ( devBOLD + L * Batch_size * (Transferblocknum-1) )  ? devBOLD : *pointer + L * Batch_size; 
}
void updateDevBOLD(int jj, real__t* devBOLD, real__t* BOLD_t, int L, int Batch_size, int* hosPoint, int* count, real__t* pointer, const int Num_Blocks)
{
	if (*count + 1 < Num_Blocks)
	{
		real__t* devJiAddr = pointer;
		if (*hosPoint == Num_Blocks - 1)
		{
			(*count)++;
			(*hosPoint) =  (*count);	
		}
		else
			(*hosPoint) ++;
	real__t* hosJiAddr = BOLD_t + (*hosPoint) * Batch_size * L;
	checkCudaErrors( cudaMemcpy(devJiAddr, hosJiAddr, sizeof(real__t) * L * Batch_size, cudaMemcpyHostToDevice) );
	}
}

/* 
function CorMat_spa2rth
Calculate the correlation threshold for 
each sparsity threshold by a histogram method.
*/
real__t CorMat_spa2rth_blocking(string OutCor, real__t * BOLD_t, int N,  int L, int  Batch_size,real__t *s_thresh, real__t* result, int num_spa, const int Transferblocknum)
{
	const int Num_Blocks = (N + Batch_size - 1) / Batch_size;
	decltype(N) N0 = Num_Blocks * Batch_size;
	size_t L2 = L * ( L - 1 ) / 2;
	unsigned long long *zeroAmount = new unsigned long long [num_spa];
	for (uint__t i = 0; i < num_spa; i++)
	{
		unsigned long long amount = N * s_thresh[i] * (N-1)  / 100.0 ;///2.0 ;
#ifdef myLiteDebug
		cout<<"amount:"<<amount<<endl;
#endif
		//amount += amount%2;
		zeroAmount[i] = N;
		zeroAmount[i] *= N ;
		zeroAmount[i] -= amount;
		//cout<<zeroAmount[i]<<endl;
	}
		
	
	long long subtractor = (long long)N0 * N0 - (long long)N * N;
	uint__t Num_Bins = 1.0 / width + 1; //take care!
#ifdef myLiteDebug
	cout<<"Num_Bins:"<<Num_Bins<<endl;
#endif
	uint__t position = 0;
	// transposing the BOLD signal
	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;               //optimize:bitmap+ for 1 -1 0; 
	real__t * devBOLD, * devCormat, * devBOLD_x, * devBOLD_y;
	size_t * tiecount; //perhaps _y _x tiecount should be integer type.
	real__t * pointer;
	int count = 0; //used to control pointer
	int hosPoint = Transferblocknum - 1;
	//uint__t *devhisto;
	checkCudaErrors ( cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * Batch_size * Transferblocknum) ) ;
	checkCudaErrors ( cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size) );
	checkCudaErrors ( cudaMalloc ((void**)&tiecount, sizeof(size_t) * N0) );
	checkCudaErrors ( cudaMalloc ((void**)&devBOLD_x, sizeof(real__t) * L2 * Batch_size) );
	checkCudaErrors ( cudaMalloc ((void**)&devBOLD_y, sizeof(real__t) * L2 * Batch_size) );
	checkCudaErrors( cudaMemcpy(devBOLD, BOLD_t, sizeof(real__t) * L * Batch_size * Transferblocknum, cudaMemcpyHostToDevice) );
	pointer = devBOLD;	
	stat = cublasCreate(&handle) ;
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
	
	cout<<"generating histogram..."<<endl;
	cout<<"block number per row: "<<Num_Blocks<<endl;
	thrust::device_vector<long long> histogram(Num_Bins,0);
	clock_t time;
	time = clock();
	//uint__t blocknum = ( Batch_size * Batch_size + thread_num -1 ) / thread_num;
	bool flag = true;
	size_t shareSize = sizeof(real__t) * L + sizeof(size_t) * L;
	dim3 Grid(Batch_size/thread_num2D, Batch_size/thread_num2D), Block(thread_num2D,thread_num2D);
	real__t* biPointer =  devBOLD_y;
	for (int ii = 0; ii < Num_Blocks; ii++)
	{
		if ( ii > 0 )
			flag = false;
		for (int jj = ii; jj < Num_Blocks; jj++)
		{
#ifdef myDebug
	cout<<"loop flag:"<<ii<<" "<<jj<<endl;
#endif
			biPointer = ii == jj ? devBOLD_y : devBOLD_x;
#ifdef figure
			int thread_num = L;
			int block_num = 180;
#endif
			pre_process <<<block_num, thread_num, shareSize>>>(biPointer, pointer, L, L2, Batch_size, tiecount + jj * Batch_size, flag);
	        checkCublasErrors( cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L2,  &alpha, devBOLD_y, L2, biPointer, L2, &beta, devCormat, Batch_size) ); //Is this really right?
			dividedByDenominatorAndStandardedKernelWith2DBlock<<<Grid, Block>>>(devCormat, Batch_size, L, tiecount + ii * Batch_size, tiecount + jj * Batch_size, ii==jj);
			updateDevBOLD(jj, devBOLD, BOLD_t, L, Batch_size, &hosPoint, &count, pointer, Num_Blocks);
			updataPointer(&pointer, devBOLD, L, Batch_size, Transferblocknum);
#ifdef myDebug
			//cout<<"loop flag:"<<ii<<" "<<jj<<endl;
			gpuOutput(devCormat, sizeof(real__t) * Batch_size * Batch_size, OutCor, true);
			if (jj==1)
			{
				thrust::device_vector<size_t> dt(tiecount,tiecount+129);
				print_vector("common!",dt);
			}
#endif
			if (ii==jj)
				dense_histogram(thrust::device_pointer_cast(devCormat),Batch_size * Batch_size, width, Num_Bins, histogram); //Somehow num_bins is generated wrongly in histogram.h files. So just transfer it directly.
			else
			    dense_histogram2(thrust::device_pointer_cast(devCormat),Batch_size * Batch_size, width, Num_Bins, histogram);//difference: Multiply temphisto by 2; Somehow num_bins is generated wrongly in histogram.h files. So just transfer it directly.
		}
	}
	
	substract(histogram, subtractor); //subtract extra 0 emerged by ending batch
	for (size_t i = 0; i < num_spa; i++)
	{
		//cout<<zeroAmount[i]<<endl;
		position = thrust::upper_bound(histogram.begin(),histogram.end(),zeroAmount[i]) - histogram.begin();//flag~ N*N
		result[i] = width * position + width/2.0;
#ifdef myLiteDebug
		cout<<"threshold:"<<result[i]<<endl;
		cout<<position<<endl;
		/*cout<<"histogram[position-1]:"<<histogram[position-1]<<endl;
		cout<<"histogram[position]:"<<histogram[position]<<endl;
		cout<<"histogram[position+1]:"<<histogram[position+1]<<endl;*/
#endif
	}
	//display and put out
	time = clock() - time;
	cout<<"histogram time: "<<time<<"ms"<<endl;
	
	checkCudaErrors ( cudaFree (devBOLD));
	checkCudaErrors ( cudaFree (devCormat));
	checkCudaErrors ( cudaFree (devBOLD_x));
	checkCudaErrors ( cudaFree (devBOLD_y));
	checkCudaErrors ( cudaFree (tiecount));

	/*thrust::adjacent_difference(histogram.begin(), histogram.end(), histogram.begin());
	print_vector("histogram:",histogram);*/
	thrust::device_vector<long long>().swap(histogram);
	
	delete[] zeroAmount;
	stat = cublasDestroy(handle);
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
		
}

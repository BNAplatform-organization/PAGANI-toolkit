#include "cublas_v2.h"
#include "cuda_runtime.h"
#include "memory.h"
#include <iostream>
#include <ctime>
# include <fstream>
#include <vector>
#include <Windows.h>
#include<iomanip>
#include"histogram.h"
#include <thrust/binary_search.h>

#include "help_func.cuh"
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/functional.h>
#include <cmath>

using namespace std;

#define ep  1e-6  //third question
//#define width  0.1//0.0001//best to be the multiples
#define width  0.000001//0.0001//best to be the multiples

#pragma comment(lib,"cublas.lib")
typedef float real__t;
typedef unsigned int uint__t;

const int thread_num = 1024; 
const int block_num = 30; 
bool myfunction (real__t i,real__t j) { return (i>j); }
#define TOM(byteValue) (byteValue/1024/1024)

__global__ void standardKernel(real__t* devCormat, int Batch_size, bool diagnoal)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = blockDim.x * gridDim.x;
	while(i<Batch_size*Batch_size)
	{
		if (!(devCormat[i]>0.0f && devCormat[i]<(1+ep)) ) 
		{
			devCormat[i] = 0;
		}
		if(diagnoal==true)
		{
			if (i%(Batch_size+1)==0)
		    {
			devCormat[i] = 0;
		    }
		}
		i += offset;
	}
}

/* 
function CorMat_spa2rth
Calculate the correlation threshold for 
each sparsity threshold by a histogram method.
*/
real__t CorMat_spa2rth(string OutCor, real__t * BOLD_t, int N,  int L, int  Batch_size, real__t *s_thresh, real__t* result, int num_spa)
{
	//cout<<"*s_thresh: "<<*s_thresh<<endl;
	//real__t * BOLD_t1, * BOLD_t2;// * tempout;
	const int Num_Blocks = (N + Batch_size - 1) / Batch_size;
	decltype(N) N0 = Num_Blocks * Batch_size;
	
	//amount += amount%2;
	//modify1 zeroAmount to an array
	//uint__t num_spa = sizeof(s_thresh) / sizeof(real__t);
	long long *zeroAmount = new long long [num_spa];
	for (uint__t i = 0; i < num_spa; i++)
	{
		long long amount = N * s_thresh[i] * (N-1)  / 100.0 ;///2.0 ;
		//amount += amount%2;
	
		zeroAmount[i] = N;
		zeroAmount[i] *= N ;//
		zeroAmount[i] -= amount;
		//cout<<zeroAmount[i]<<endl;
	/*	zeroAmount[i] = (long long)N * s_thresh[i] * (N-1)  / 100.0 ;
		zeroAmount[i] -= ((long long)N * N);
		zeroAmount[i] *= -1;*/
	}
		
	
	long long subtractor = (long long)N0 * N0 - (long long)N * N;
	//cout<<N0 * N0<<endl;
	//cout<<N * N<<endl;
	//cout<<subtractor<<endl;
	//uint__t invaccount = N - account;
	uint__t Num_Bins = 1.0 / width + 1; //take care!
	uint__t position = 0;

	// transposing the BOLD signal
	
	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;
	real__t * devBOLD, * devCormat;// * devCormatLower, * devCormatPacked;
	
	
	checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * N0)) ;
	checkCudaErrors (cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size));	
	//checkCudaErrors (cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size));	

	stat = cublasSetMatrix(N0, L, sizeof(real__t), BOLD_t, N0, devBOLD, N0);
	
	stat = cublasCreate(&handle) ;
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
	cout<<"generating histogram..."<<endl;
	cout<<"block number per row: "<<Num_Blocks<<endl;
	const float alpha = 1.0;
	const float beta = 0;	
	thrust::device_vector<long long> histogram(Num_Bins,0);
	clock_t time;
	time = clock();
	//uint__t blocknum = ( Batch_size * Batch_size + thread_num -1 ) / thread_num;
	for (int kk = 0, ii = 0; ii < Num_Blocks; ii++)
	{
		for (int jj = ii; jj < Num_Blocks; jj++)
		{
			//time = clock();
			stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L,  &alpha, devBOLD + ii * Batch_size * L, L, devBOLD + jj * Batch_size * L, L, &beta, devCormat, Batch_size);
			if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;
			//1clear NaN and negative number
			standardKernel<<<block_num,thread_num>>>(devCormat, Batch_size,ii==jj);
			//2.statistics;the results of non-diagnoal batch shall multiply by 2 accordingly
			//thrust::device_vector<uint__t> temphisto(Num_Bins,0); //dontforgettofree!
			thrust::device_ptr<real__t> dev_ptr(devCormat);
			
			/*if (ii==jj)
				dense_histogram(dev_ptr,Batch_size * Batch_size, width,temphisto, histogram);
			else
				dense_histogram2(dev_ptr,Batch_size * Batch_size, width,temphisto, histogram);*/

			//it seems that new approach reach higher speed.
			if (ii==jj)
				dense_histogram_new(dev_ptr,Batch_size * Batch_size, width, Num_Bins, histogram);
			else
				dense_histogram2_new(dev_ptr,Batch_size * Batch_size, width, Num_Bins, histogram);
			
			//thrust::device_vector<uint__t>().swap(temphisto);
			
			//thrust::raw_pointer_cast(dev_ptr);//have a try
			//time = clock() - time;
			//cout<<"thrust::histogram time: "<<time<<"ms"<<endl;
			//cout<<"loop flag "<<"ii: "<<ii<<" jj: "<<jj<<endl;
		}
		//cout<<"Fulfill the "<<ii+1<<"th block."<<endl;
	}

	substract(histogram, subtractor);
	
	for (size_t i = 0; i < num_spa; i++)
	{
		//cout<<zeroAmount[i]<<endl;
		position = thrust::upper_bound(histogram.begin(),histogram.end(),zeroAmount[i]) - histogram.begin();//flag~ N*N
		//cout<<position<<endl;
		result[i] = width * position + width/2.0;
#ifdef myLiteDebug
		cout<<"threshold:"<<result[i]<<endl;
		cout<<"histogram[position-1]:"<<histogram[position-1]<<endl;
		cout<<"histogram[position]:"<<histogram[position]<<endl;
		cout<<"histogram[position+1]:"<<histogram[position+1]<<endl;
#endif
	}

		//another interface is needed to return the subscript
		//display and put out
	time = clock() - time;
	cout<<"histogram time: "<<time<<"ms"<<endl;
	//checkCudaErrors (cudaFree(thrust::raw_pointer_cast(histogram.data())));
	thrust::device_vector<long long>().swap(histogram);
	checkCudaErrors ( cudaFree (devBOLD));
	checkCudaErrors ( cudaFree (devCormat));

	delete[] zeroAmount;
	stat = cublasDestroy(handle);
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
		
}

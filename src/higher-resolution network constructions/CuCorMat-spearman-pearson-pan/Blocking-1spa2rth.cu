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
extern bool myfunction (real__t i,real__t j);// { return (i>j); }
#define TOM(byteValue) (byteValue/1024/1024)

extern __global__ void standardKernel(real__t* devCormat, int Batch_size, bool diagnoal);


void updateOperationii(real__t** operationii, real__t* devBOLDii, int L, int Batch_size, int ii)
{
	*operationii = ii%2 ? devBOLDii + L * Batch_size : devBOLDii;
}
void updateDevBOLDii(int ii, real__t* devBOLDii, real__t* BOLD_t, int L, int Batch_size, const int Num_Blocks)
{
	if ((ii + 1)<Num_Blocks)
	{
		real__t* devIiAddr = ii%2 ? devBOLDii : devBOLDii + L * Batch_size; 
		real__t* hosIiAddr = BOLD_t + (ii + 1) * Batch_size * L;
		checkCudaErrors( cudaMemcpy(devIiAddr, hosIiAddr, sizeof(real__t) * L * Batch_size, cudaMemcpyHostToDevice) );
	}
}
void updatedevBOLDjj(int jj, real__t* devBOLDjj, real__t* BOLD_t, int L, int Batch_size, int* hosPoint, int* count, real__t* operationjj, const int Num_Blocks)
{
	if (*count + 1 < Num_Blocks)
	{
		real__t* devJiAddr = operationjj;
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
void updateOperationjj(real__t** operationjj, real__t* devBOLDjj, int L, int Batch_size, const int numjj)
{
	*operationjj =  *operationjj ==  ( devBOLDjj + L * Batch_size * (numjj-1) )  ? devBOLDjj : *operationjj + L * Batch_size; 
	
}


#ifdef myDebug
void gpuOutput( real__t* gpuAddr, unsigned int byteNo, string OutCor, bool nameFlag)
{
	ofstream fout;
	real__t* cpuAddr = new real__t[byteNo/sizeof(real__t)];
	checkCudaErrors ( cudaMemcpy(cpuAddr, gpuAddr, byteNo, cudaMemcpyDeviceToHost) );
	string filename;
	if(!nameFlag)
		filename = OutCor.append("compareW_wrong.matrix");
	else
		filename = OutCor.append("compareW_right.matrix");
	fout.open(filename.c_str(), ios::binary | ios::out);
	if (!fout)
	{
		cout<<"create outfile(gpu) unsuccessfully. error code:  "<<"fighting!"<<endl;		
		system("pause");
	}	
	fout.write((const char*)cpuAddr,byteNo);
	fout.close();
	delete[] cpuAddr;
}

#endif


/* 
function CorMat_spa2rth
Calculate the correlation threshold for 
each sparsity threshold by a histogram method.
*/
real__t CorMat_spa2rth_blocking(string OutCor, real__t * BOLD_t, int N,  int L, int  Batch_size,real__t *s_thresh, real__t* result, int num_spa, const int Transferblocknum)
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
	real__t * devBOLDii, * devBOLDjj; //only refer to head address; actually is a const.
	real__t * operationii, * operationjj;
	
	checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * Batch_size * Transferblocknum)) ;
	checkCudaErrors (cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size));	
	devBOLDii = devBOLD;
	devBOLDjj = devBOLD + L * Batch_size * 2;
	operationii = devBOLDii; 
	operationjj = devBOLDjj; 
	const int numjj = Transferblocknum - 2;
	int count = 0; //used to control devboldjj
	int hosPoint = numjj - 1;

	checkCudaErrors( cudaMemcpy(devBOLDii, BOLD_t, sizeof(real__t) * L * Batch_size, cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(devBOLDjj, BOLD_t, sizeof(real__t) * L * Batch_size * numjj, cudaMemcpyHostToDevice) );
	
	stat = cublasCreate(&handle) ;
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
	cout<<"generating histogram (blocked transmission)..."<<endl;
	cout<<"block number per row: "<<Num_Blocks<<endl;
	const float alpha = 1.0;
	const float beta = 0;	
	thrust::device_vector<long long> histogram(Num_Bins,0);
	clock_t time;
	time = clock();
	//uint__t blocknum = ( Batch_size * Batch_size + thread_num -1 ) / thread_num;
	for (int kk = 0, ii = 0; ii < Num_Blocks; ii++)
	{
		updateOperationii(&operationii, devBOLDii, L, Batch_size, ii);
		updateDevBOLDii(ii,devBOLDii, BOLD_t, L, Batch_size, Num_Blocks);
		
		cout<<"loop flag "<<"ii: "<<ii<<endl;
		for (int jj = ii; jj < Num_Blocks; jj++)
		{
			//time = clock();
			//checkCublasErrors( cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L,  &alpha, devBOLD + ii * Batch_size * L, L, devBOLD + jj * Batch_size * L, L, &beta, devCormat, Batch_size) );
			checkCublasErrors( cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L,  &alpha, operationii, L, operationjj, L, &beta, devCormat, Batch_size) );
			updatedevBOLDjj(jj,devBOLDjj, BOLD_t, L, Batch_size, &hosPoint, &count, operationjj, Num_Blocks);
			updateOperationjj(&operationjj, devBOLDjj, L, Batch_size, numjj);
			//1clear NaN and negative number
			standardKernel<<<block_num,thread_num>>>(devCormat, Batch_size,ii==jj);
			//2.statistics;the results of non-diagnoal batch shall multiply by 2 accordingly
			//thrust::device_vector<uint__t> temphisto(Num_Bins,0); //dontforgettofree!
			thrust::device_ptr<real__t> dev_ptr(devCormat);
			
			/*if (ii==jj)
				dense_histogram(dev_ptr,Batch_size * Batch_size, width,temphisto, histogram);
			else
				dense_histogram2(dev_ptr,Batch_size * Batch_size, width,temphisto, histogram);*/

#ifdef myDebug
		gpuOutput( devCormat, sizeof(real__t) * Batch_size * Batch_size, OutCor, false);
#endif
			//it seems that new approach reach higher speed.
			if (ii==jj)
				dense_histogram_new(dev_ptr,Batch_size * Batch_size, width, Num_Bins, histogram);
			else
				dense_histogram2_new(dev_ptr,Batch_size * Batch_size, width, Num_Bins, histogram);
			//thrust::device_vector<uint__t>().swap(temphisto);
			//thrust::raw_pointer_cast(dev_ptr);//have a try
			//time = clock() - time;
			//cout<<"thrust::histogram time: "<<time<<"ms"<<endl;
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
#ifdef myDebug
		cout<<"threshold:"<<result[i]<<endl;
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

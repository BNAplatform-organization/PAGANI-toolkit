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

using namespace std;

#define ep  1e-6  //third question
#define width  0.000001//0.0001//best to be the multiples 

#pragma comment(lib,"cublas.lib")
typedef float real__t;
typedef unsigned int uint__t;

const int thread_num = 256; //maybe redefinition
const int block_num = 48;     //(bat*bat+thrednum-1)/threadnum
const int blocksize = 1024*1024*48;
/*
bool IsNumber(double x)
{
	return (x == x);
}
 bool IsFiniteNumber(double x)
{
	return (x <= DBL_MAX && x >= -DBL_MAX);
}
*/
bool myfunction (real__t i,real__t j) { return (i>j); }
//void select(vector<real__t>::iterator A,long long n,long long k);
#define TOM(byteValue) (byteValue/1024/1024)
#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)
inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
	if(cudaSuccess != err)
	{
		fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString( err ) );
		//exit(-1);
		system("pause");
	}
}
__global__ void standardKernel(real__t* devCormat, int Batch_size, bool diagnoal)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = blockDim.x * gridDim.x;
	while(i<Batch_size*Batch_size) 
	{
		if ((!(devCormat[i] == devCormat[i]))||(!(devCormat[i] <= DBL_MAX && devCormat[i] >= -DBL_MAX))||devCormat[i]<0) 
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
void MatrixMultiplication(real__t * BOLD_t1, real__t * BOLD_t2,real__t * out,int Batch_size,int L);
real__t CorMat_spa2rth(string OutCor, real__t * BOLD, int N, int L, int Batch_size,real__t *s_thresh, clock_t * aggregate)
{
	real__t * BOLD_t1, * BOLD_t2, * tempout;
	const int Num_Blocks = (N + Batch_size - 1) / Batch_size;
	uint__t N0 = Num_Blocks * Batch_size;
	uint__t amount = N * (*s_thresh) * (N-1)  / 100.0 ;///2.0 ;
	//amount += amount%2;
	uint__t zeroAmount = N * N - amount;
	cout<<"nonzero numbers: "<<amount<<endl;
	//uint__t invaccount = N - account;
	uint__t Num_Bins = 1.0 / width + 1; //take care! 
	uint__t position = 0;

	// transposing the BOLD signal
	real__t * BOLD_t = new real__t [L * N0];
	tempout = new real__t[Batch_size * Batch_size];
	memset(BOLD_t, 0, sizeof(real__t) * L * N0);
	for (int i = 0; i < L; i ++)
		for (int j = 0; j < N; j++)
		{
			BOLD_t[j * L + i] = BOLD[i * N + j];
		}
		
		for (long i = L * N; i < L * N0; i++)
		{
			BOLD_t[i] = 0;
		}	
		
		// Normalize
		for (int i = 0; i < N; i++)
		{
			real__t * row = BOLD_t + i * L;
			double sum1 = 0, sum2 = 0;
			for (int l = 0; l < L; l++)
			{
				sum1 += row[l];
			}
			sum1 /= L;
			for (int l = 0; l < L; l++)
			{
				sum2 += (row[l] - sum1) * (row[l] - sum1);
			}
			sum2 = sqrt(sum2);
			for (int l = 0; l < L; l++)
			{
				row[l] = (row[l] - sum1) / sum2;;
			}
		}
		cudaError_t cudaStat;
		cublasStatus_t stat;
		cublasHandle_t handle;
		real__t * devBOLD, * devCormat;// * devCormatLower, * devCormatPacked;
		//uint__t *devhisto;
		cudaStat = cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * N0) ;
		//cudaError_t err = cudaGetLastError();
        //printf("%s\n",cudaGetErrorString(err));
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
		cudaStat = cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
		stat = cublasSetMatrix(N0, L, sizeof(real__t), BOLD_t, N0, devBOLD, N0);
		stat = cublasCreate(&handle) ;
		if (stat != CUBLAS_STATUS_SUCCESS)
			return stat;
		delete []BOLD_t;

		//是指GPU block的个数！
		cout<<"generating histogram..."<<endl;
		cout<<"block numbers: "<<Num_Blocks<<endl;
		const float alpha = 1.0;
		const float beta = 0;
		const float pbeta = -1.0;
		vector< vector<real__t> > bin;
		thrust::device_vector<uint__t> histogram(Num_Bins,0); 
		clock_t time;
		time = clock();
		//uint__t blocknum = ( Batch_size * Batch_size + thread_num -1 ) / thread_num;
		for (int kk = 0, ii = 0; ii < Num_Blocks; ii++)
		{
			for (int jj = ii; jj < Num_Blocks; jj++)
			{
				  
				//BOLD_t1 = BOLD_t + ii * Batch_size * L;
				//BOLD_t2 = BOLD_t + jj * Batch_size * L;
   			//	real__t *out = new real__t[Batch_size * Batch_size];

#ifdef CPUCormat
                MatrixMultiplication_s(BOLD_t1, BOLD_t2, out, Batch_size,L);//need modify as well.
#else            
				//time = clock();
				stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L,  &alpha, devBOLD + jj * Batch_size * L, L, devBOLD + ii * Batch_size * L, L, &beta, devCormat, Batch_size);
				if (stat != CUBLAS_STATUS_SUCCESS)
					return stat;
				//1clear NaN and negative number 
				standardKernel<<<block_num,thread_num>>>(devCormat, Batch_size,ii==jj);
				//2.statistics;the results of non-diagnoal batch shall multiply by 2 accordingly
				thrust::device_vector<uint__t> temphisto(Num_Bins,0); //dontforgettofree!
				thrust::device_ptr<real__t> dev_ptr(devCormat);
				
				if (ii==jj)
				{
					dense_histogram(dev_ptr,Batch_size * Batch_size, width,temphisto, histogram);
				}
				else
				{
					dense_histogram2(dev_ptr,Batch_size * Batch_size, width,temphisto, histogram);
				}
				thrust::device_vector<uint__t>().swap(temphisto);
				thrust::raw_pointer_cast(dev_ptr);//have a try
				//time = clock() - time;
				//cout<<"thrust::histogram time: "<<time<<"ms"<<endl;
				cout<<"loop flag "<<"ii: "<<ii<<" jj: "<<jj<<endl;
#endif
			}
			cout<<"Fulfill the "<<ii+1<<"th block."<<endl;
		}
		
		uint__t subtractor = N0 * N0 - N * N;
		substract(histogram, subtractor);
		position = thrust::upper_bound(histogram.begin(),histogram.end(),zeroAmount) - histogram.begin();
		real__t result = width * position + width/2.0; 
		//another interface is needed to return the subscript		
		    //display and put out 
			time = clock() - time;
			cout<<"histogram time: "<<time<<"ms"<<endl;
			*aggregate = time;
			//time_t nowTime;
			/*
			unsigned int FreeMem = 0;
			MEMORYSTATUS MemStat;
			MemStat.dwLength = sizeof(MEMORYSTATUS);
			GlobalMemoryStatus(&MemStat);
			FreeMem = TOM(MemStat.dwAvailPhys);
			cout << "bytes of physical memory: " << TOM(MemStat.dwTotalPhys) <<"M" <<endl;
			cout << "percent of memory in use: " << MemStat.dwMemoryLoad <<"%" <<endl;
			cout << "free physical memory bytes: " << TOM(MemStat.dwAvailPhys) <<"M" <<endl;
			*/
			cudaFree (devBOLD); 

			cudaFree (devCormat);
			stat = cublasDestroy(handle);
			if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;
		//	delete []BOLD_t;
			return result;
}
real__t interval(void)
{
	return width;
}
//void MatrixMultiplication_s(real__t * BOLD_t1, real__t * BOLD_t2,real__t * out,int Batch_size,int L)//do not announce
//{
//	long kk = 0;
//	for (int k = 0; k < Batch_size; k++)
//	{
//		for (int i = 0; i < Batch_size; i++)
//		{   
//			double sum3 = 0.0;
//			for (int j = 0; j < L; j++)
//			{
//				sum3 += 1.0*BOLD_t1[k*L+j] * BOLD_t2[i*L+j];
//			}
//			out[kk++] = sum3;
//		}
//	}
//	
//}
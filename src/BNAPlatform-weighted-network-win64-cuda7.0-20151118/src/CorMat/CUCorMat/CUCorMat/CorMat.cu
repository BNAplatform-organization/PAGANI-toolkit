#include "cublas_v2.h"
#include "cuda_runtime.h"
#include "memory.h"
#include <iostream>
#include <ctime>
using namespace std;

typedef float real__t;
typedef unsigned int uint__t;

const int thread_num = 256;
const int block_num = 48;
const int blocksize = 1024*1024*48;

void select(real__t *A,long long n,long long k);

int CorMat_gpu(real__t * Cormat, real__t * BOLD, int N, int L, int Batch_size)
{
	real__t * BOLD_t1, * BOLD_t2, * out, * tempout;
	int Num_Blocks = (N + Batch_size - 1) / Batch_size;
	uint__t N0 = Num_Blocks * Batch_size;

	// transposing the BOLD signal
	real__t * BOLD_t = new real__t [L * N0];
	tempout = new real__t[Batch_size * Batch_size];
	memset(BOLD_t, 0, sizeof(real__t) * L * N0);
	for (int i = 0; i < L; i ++)
		for (int j = 0; j < N; j++)
		{
			BOLD_t[j * L + i] = BOLD[i * N + j];
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

		//// column major in every block
		//real__t * BOLD_t_col = new real__t [L * N0];
		//for (int k = 0; k < Num_Blocks; k++)
		//{
		//	for (int i = 0; i < Batch_size; i ++)
		//		for (int j = 0; j < L; j++)
		//		{
		//			BOLD_t_col[k * Batch_size * L + j * Batch_size + i] = BOLD_t[k * Batch_size * L + i * L + j];
		//		}
		//}

		cudaError_t cudaStat;
		cublasStatus_t stat;
		cublasHandle_t handle;
		real__t * devBOLD, * devCormat;
//		stat = cublasAlloc(L*N0, sizeof(real__t), (void**)&devBOLD);
		cudaStat = cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * N0) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
//		stat = cublasAlloc(Batch_size * Batch_size, sizeof(real__t), (void**)&devCormat);		
		cudaStat = cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
		stat = cublasSetMatrix(N0, L, sizeof(real__t), BOLD_t, N0, devBOLD, N0);
//		cudaStat = cudaMemcpy(devBOLD, BOLD_t, sizeof(real__t) * L * N0, cudaMemcpyHostToDevice);
		stat = cublasCreate(&handle) ;
		if (stat != CUBLAS_STATUS_SUCCESS)
			return stat;


		cout<<Num_Blocks<<endl;
		const float alpha = 1.0;
		const float beta = 0;
		for (int kk = 0, ii = 0; ii < Num_Blocks; ii++)
			for (int jj = ii; jj < Num_Blocks; jj++)
			{
				BOLD_t1 = BOLD_t + ii * Batch_size * L;
				BOLD_t2 = BOLD_t + jj * Batch_size * L;
				out = Cormat + (long long) kk * Batch_size * Batch_size;
				kk++;
				stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L,  &alpha, devBOLD + jj * Batch_size * L, L, devBOLD + ii * Batch_size * L, L, &beta, devCormat, Batch_size);

				if (stat != CUBLAS_STATUS_SUCCESS)
					return stat;
				
				cudaStat = cudaMemcpy(out, devCormat, sizeof(real__t) * Batch_size * Batch_size, cudaMemcpyDeviceToHost);
				if (cudaStat != cudaSuccess) 
					return cudaStat;
				
		/*		for (int i = 0; i < Batch_size; i ++)
					for (int j = 0; j < Batch_size; j++)
					{
						out[j * Batch_size + i] = tempout[i * Batch_size + j];
					}*/
				/*float *testA = new float[Batch_size * L];
				cudaMemcpy(testA, devBOLD + ii * Batch_size * L, sizeof(real__t) * Batch_size * L, cudaMemcpyDeviceToHost);
				float *testB = new float[Batch_size * L];
				cudaMemcpy(testB, devBOLD + jj * Batch_size * L, sizeof(real__t) * Batch_size * L, cudaMemcpyDeviceToHost);
				cout<<"A"<<endl;
				for (int i = 0; i < Batch_size; i++)
				{
				for (int j = 0; j < L; j++)
				{
				cout<<testA[i * L + j]<<"\t";
				}
				cout<<endl;
				}
				getchar();
				cout<<"B"<<endl;

				for (int i = 0; i < Batch_size; i++)
				{
				for (int j = 0; j < L; j++)
				{
				cout<<testB[i * L + j]<<"\t";
				}
				cout<<endl;
				}
				getchar();
				cout<<"C"<<endl;
				for (int i = 0; i < Batch_size; i++)
				{
				for (int j = 0; j < Batch_size; j++)
				{
				cout<<out[i * Batch_size + j]<<"\t";
				}
				cout<<endl;
				}
				getchar();*/
				/*	double sum3;
				for (int k = 0, i = 0; i < Batch_size; i++)
					for (int j = 0; j < Batch_size; j++)
					{
						sum3 = 0;
						for (int l = 0; l < L; l++)
						{
							sum3 += BOLD_t1[i*L+l] * BOLD_t2[j*L+l];
						}
						out[k++] = sum3;
					}*/
			}
			cudaFree (devBOLD); 
			cudaFree (devCormat);
			stat = cublasDestroy(handle);
			if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;
			delete []BOLD_t;
			return 1;
}

__global__ void partition_kernel ( real__t x, real__t *A,  double *j_l, double *j_r,const long N)
{
	__shared__ int tmp_j_l[thread_num];
	__shared__ int tmp_j_r[thread_num];
	int tmp_l = 0;
	int tmp_r = 0;
	int offset = 0;
	const int threadid =blockIdx.x*blockDim.x + threadIdx.x;

	for(offset=threadid; offset<N; offset+=blockDim.x*gridDim.x)
	{		
		if(A[offset]>x)
			tmp_l++ ;
		else if(A[offset]>x)
			tmp_r++;		
	}
	tmp_j_l[threadIdx.x] = tmp_l;
	tmp_j_r[threadIdx.x] = tmp_r;

	syncthreads();

	for(offset=1; offset+threadIdx.x<thread_num; offset*=2)
	{
			if (threadIdx.x%(2*offset)==0)  tmp_j_l[threadIdx.x]+= tmp_j_l[threadIdx.x+offset] ;
			syncthreads();
	}	
	if(threadIdx.x==0)
		j_l[blockIdx.x]=(double) tmp_j_l[0];
	if(threadIdx.x==1)
		j_l[blockIdx.x]=(double) tmp_j_r[0];
}



/*int *partition_gpu()
{}*/
real__t select_GPU(real__t *Cormat, long long M1, long long k)
{
	long long offset;

	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;
	
	int index,i=0;
	//long long k_next = 0, M_next = 0;
	//int onethread_n = blocksize/block_num/thread_num;


	real__t *Cormat_block;
	cudaStat = cudaMalloc ((void**)&Cormat_block, sizeof(real__t) * blocksize) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
			
	int blocksize_num = (M1-1+blocksize)/blocksize;
	double *j_l;
	cudaMalloc ((void**)&j_l, sizeof(double) * block_num*blocksize_num) ;
	double *j_r;
	cudaMalloc ((void**)&j_r, sizeof(double) * block_num*blocksize_num) ;
	double hj_l = 0;
	double hj_r = 0; 
	
	
	real__t left  = 0.0;
	real__t right = 1.0;
	real__t x;
	stat = cublasCreate(&handle) ;
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
	
	i=0;
	clock_t time = clock();
	
	x = (left+right)/2.0;
	for (offset = 0; offset < M1; offset += blocksize)
	{
		int size = (M1-offset > blocksize? blocksize : M1-offset);
		cudaMemcpy(Cormat_block, Cormat+offset, sizeof(real__t) * size, cudaMemcpyHostToDevice);
		partition_kernel<<<block_num,thread_num>>>(x, Cormat_block,  j_l+block_num*i, j_r+block_num*i,(long) size);
		stat = cublasDasum(handle, block_num*blocksize_num, j_l, 1, &hj_l);
		if (stat != CUBLAS_STATUS_SUCCESS)
			return stat;
		stat = cublasDasum(handle, block_num*blocksize_num, j_r, 1, &hj_r);
		if (stat != CUBLAS_STATUS_SUCCESS)
			return stat;
		i++;
	}
	
	//cudaMemcpy(hj_r, j_r, sizeof(int)*block_num*((M1-1)/blocksize), cudaMemcpyDeviceToHost);
	if(hj_l + hj_r > M1) cout<<"partition error! hj_l = "<<hj_l<<" ; hj_r = "<<hj_r<<endl;
	else if(hj_l < k && hj_r < M1-k) return x;
	else if ( hj_l < k)  right = x;	
	else left = x;
	
	time = clock()-time;
	cout<<"first partition time = "<<time<<";  hj_l = "<<hj_l<<endl;

		/*cout<<"round "<<offset/blocksize<<endl;
		
		cout<<"bound for partition "<<x<<endl; 
		for (int i = 0; i < block_num; i++)
		{
			cout<<i<<" block :"<<hj_l[i]<<"  "<<hj_r[i]<<"  "<<hj_l[i]+hj_r[i]<<endl;
			if (hj_l[i]+hj_r[i] != thread_num * onethread_n)
			{
				cout<< "partition error! " <<i<<" block;  "<<hj_l[i]<<"  "<<hj_r[i]<<endl;
				system("pause");
			}
				int j = 0;
		}*/
				//	if (offset+i*onethread_n*thread_num>=M1)
		//		break;
		
	
	return (0);
	
}






long long find_max(real__t *Cormat, long long M1)
{
	long long offset;

	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;
	long long blocksize = 1024*1024*48*8;
	int index,i=0;

	real__t *Cormat_block;
	cudaStat = cudaMalloc ((void**)&Cormat_block, sizeof(real__t) * blocksize) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
	
	stat = cublasCreate(&handle) ;
		if (stat != CUBLAS_STATUS_SUCCESS)
			return stat;

	int segnum = (M1+blocksize)/blocksize;
	real__t *tmp = new real__t [segnum];
	
	long long *tmp_index = new long long [segnum] ;
	

	for (offset = 0; offset + blocksize < M1; offset += blocksize)
	{	
		cublasSetVector(blocksize, sizeof(real__t), Cormat+offset, 1, Cormat_block, 1);
		stat = cublasIsamax(handle, blocksize, Cormat_block, 1, &index);
		tmp[i] = *(Cormat+offset+index-1);
		//cout << tmp[i]<<endl;
		tmp_index[i++] = index+offset-1;
		//cout<< tmp_index[i-1]<<endl;
	}
	cublasSetVector(M1-offset, sizeof(real__t), Cormat+offset, 1, Cormat_block, 1);
	stat = cublasIsamax(handle, M1-offset, Cormat_block, 1, &index);
	tmp[i] = *(Cormat+offset+index);
	tmp_index[i] = index+offset;
	//cout << tmp[i]<<endl;
	//cout<< tmp_index[i]<<endl;

	//real__t *tmp_gpu;
	//cudaStat = cudaMalloc ((void**)&tmp_gpu, sizeof(real__t) * segnum) ;
	//if (cudaStat != CUBLAS_STATUS_SUCCESS) 
	//		return cudaStat;
	
	//cublasSetVector(segnum, sizeof(real__t), tmp, 1, tmp_gpu, 1);
	//stat = cublasIsamax(handle, segnum, tmp_gpu, 1, &index);
	real__t max_r = tmp[0];
	for (i = 1; i<segnum; i++)
		if (tmp[i] > max_r)
		{	max_r = tmp[i]; index = i;  }  

	return tmp_index[index];
}


	
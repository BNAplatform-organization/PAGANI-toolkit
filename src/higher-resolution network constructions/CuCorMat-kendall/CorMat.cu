#include "cublas_v2.h"
#include "cuda_runtime.h"
#include "memory.h"

#include <iostream>
#include <ctime>

using namespace std;

typedef float real__t;
typedef unsigned int uint__t;

#define WARP 32
const int thread_num = 256;
//const int Bv_size = 256;
const int block_num = 48;
const int blocksize = 1024*1024*48;

#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)
inline void __getLastCudaError(const char *errorMessage, const char *file, const int line )
{
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
			file, line, errorMessage, (int)err, cudaGetErrorString( err ) );
		exit(-1);
	}
}

__global__ void pre_process (real__t * devBOLD, real__t * BOLD_ori, int L, int L2, int Batch_size, real__t * tiecount, bool sumtieflag )
{
	__shared__ real__t BOLD_v[thread_num]; //wrong!supposed to load entire sequence,but this kernel can only handle sequence with 1024 length, alterstion needed later  
	__shared__ int tcount[thread_num];  //sole for each sequence.
	int tid = threadIdx.x; //each thead in distinct block
	int current_offs;
              //each block have a v   
	for(int v = blockIdx.x; v < Batch_size; v += gridDim.x)
	{
		long long offset_v_obj = v*L2; //??
	
		int tie_count = 0;
		real__t tmp = 0;
		
		if (tid<L)       //not load completely??
			BOLD_v[tid] = BOLD_ori[v*L+tid];//distinct block handle distinct quantity.
		else
			BOLD_v[tid] = 0; 
	
		syncthreads();
	
		for (int ii = tid/WARP; ii <= (L-2)/2; ii+=thread_num/WARP) //question: two loop? //elements NO. 1 warp for 1 elements
		{
			current_offs = (2*L-ii-1)*ii/2;//offset address fornula: (L- 1 + L - ii) * ii / 2.0 
			for (int j = tid%WARP; j < L-1-ii; j+=WARP) //number of pair of each elements
			{
				tmp = (real__t)(BOLD_v[j+ii+1]>BOLD_v[ii]) - (BOLD_v[j+ii+1]<BOLD_v[ii]); // greater than benchmark 1; less than benchmark -1; equal to benchmark:0;
				tie_count += real__t (tmp==0);
				devBOLD[offset_v_obj + current_offs + j] =  tmp;			
			}		
		}
	                 
		for (int ii = L-2-tid/WARP; ii >(L-2)/2; ii-=thread_num/WARP) 
		{
			current_offs = (2*L-1-ii)*ii/2;
			for (int j = tid%WARP; j < L-1-ii; j+=WARP)
			{
				tmp = (real__t) (BOLD_v[ii+j+1]>BOLD_v[ii]) - (BOLD_v[ii+j+1]<BOLD_v[ii]);
				tie_count += real__t (tmp==0);
				devBOLD[offset_v_obj + current_offs + j] = tmp;			
			}
		}
		syncthreads();

		if(sumtieflag) // two definition; in order to be immune to repetitive calculation.
		{		
			tcount[tid] = tie_count; 
			syncthreads();
			for (int i = thread_num/2; i > 0; i /= 2) //add together like tree
			{
				if (tid<i) tcount[tid] += tcount[tid + i];
				syncthreads();
			}
			if (tid==0)
				tiecount[v] = (real__t) tcount[tid];
		}		
		syncthreads();
	}
}

/*
__global__ void calctie(real__t *vec, real__t *tiecount, int L2, int Batch_size )
{
	__shared__ int tmp[thread_num];
	int tid = threadIdx.x;
	int tcount = 0;
	//int tidy = threadIdx.x/WARP;
	//int n_per_block = thread_num/WARP;
	for (int v = blockIdx.x; v<Batch_size; v+=gridDim.x )
	{
		tmp[threadIdx.x] = 0;
		for(int i = tid; i < L2; i+=thread_num)
		{
			tcount+=(vec[v*L2+i]==0);
		}
		tmp[tid] = tcount;
		syncthreads();

		for (int j = thread_num/2; j> 0; j = j/2)
		{
			if(tid <j) tmp[tid] += tmp[tid + j];
			syncthreads();
		}

		if (tid==0)
			tiecount[v] = (real__t) tmp[tid];
		syncthreads();
	}
}*/

int CorMat_gpu(real__t * Cormat, real__t * BOLD, int N, int L, int Batch_size,real__t * tie_count)
{
	//real__t *BOLD_t1, *BOLD_t2;
	real__t * out, * tempout;
	int L2 = L*(L-1)/2;
	int Num_Blocks = (N + Batch_size - 1) / Batch_size;
	uint__t N0 = Num_Blocks * Batch_size;
	

	// transposing the BOLD signal
	real__t * BOLD_t = new real__t [L * N0];
	//tempout = new real__t[Batch_size * Batch_size];
	memset(BOLD_t, 0, sizeof(real__t) * L * N0);
	for (int i = 0; i < L; i ++)
		for (int j = 0; j < N; j++)
		{
			BOLD_t[j * L + i] = BOLD[i * N + j];
		}

		// Normalize
	/*
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
	}*/

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
		real__t * devBOLD_x, * devBOLD_y, * devBOLD_ori, * devCormat, * tiecount;
//		stat = cublasAlloc(L*N0, sizeof(real__t), (void**)&devBOLD);
		cudaStat = cudaMalloc ((void**)&devBOLD_ori, sizeof(real__t) * L * N0) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
//		stat = cublasAlloc(Batch_size * Batch_size, sizeof(real__t), (void**)&devCormat);		
		cudaStat = cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
		stat = cublasSetMatrix(N0, L, sizeof(real__t), BOLD_t, N0, devBOLD_ori, N0);
//		cudaStat = cudaMemcpy(devBOLD_ori, BOLD_t, sizeof(real__t) * L * N0, cudaMemcpyHostToDevice);
		
		cudaStat = cudaMalloc ((void**)&tiecount, sizeof(real__t) * N0) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
		cudaStat = cudaMalloc ((void**)&devBOLD_x, sizeof(real__t) * L2 * Batch_size) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
		cudaStat = cudaMalloc ((void**)&devBOLD_y, sizeof(real__t) * L2 * Batch_size) ;
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;		

		/******************************************/
		/*
		clock_t time = clock();
		real__t * temp_host = new real__t [L2*Batch_size];
		long long c = 0;
		for (int i = 0; i < Batch_size; i++)
			for (int j = 0; j < L-1; j++)
				for (int k = j+1; k < L; k++)
					temp_host[c++] = (BOLD_t[i*L+j]<BOLD_t[i*L+k]) - (BOLD_t[i*L+j]>BOLD_t[i*L+k]);
		cout<<"check c:"<<c-L2*Batch_size<<endl;
		time = clock()-time;
		cout<<"CPU time:"<<time<<endl;
		//cudaMemset(devBOLD_x, 0, sizeof(real__t)*L2*Batch_size);
		//for (int i = 0; i < Batch_size/32; i++)
		time = clock();
		pre_process <<<64, thread_num>>>(devBOLD_x, devBOLD_ori, L, L2,Batch_size); 
		getLastCudaError("Kernel execution failed");
		real__t * temp_device = new real__t[L2*Batch_size];
		time = clock()-time;
		cout<<"GPU time:"<<time<<endl;
		cudaStat = cudaMemcpy(temp_device, devBOLD_x, sizeof(real__t) * L2 * Batch_size, cudaMemcpyDeviceToHost);
		if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;

		
		int check;
		for (c=0; c<L2*Batch_size; c++)
			if (temp_host[c]!=temp_device[c]) 
			{   cout<<c/L2<<"  "<<c%L2<<";  host:"<<temp_host[c]<< ";  device:"<<temp_device[c]<<endl;
				cin>>check;
			}

		//cout<<"check device host: "<<check<<endl;
		*/		
		/****************************************/

		stat = cublasCreate(&handle) ;
		if (stat != CUBLAS_STATUS_SUCCESS)
			return stat;

		const float alpha = 1.0;
		const float beta = 0;
		for (int kk = 0, ii = 0; ii < Num_Blocks; ii++)
		{                                                //iith dataWidth == batch_size
			pre_process <<<64, thread_num>>>(devBOLD_y, devBOLD_ori+ii * Batch_size * L, L, L2,Batch_size, tiecount+ii*Batch_size, true);
			out = Cormat + kk * Batch_size * Batch_size;
			kk++;

			stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L2,  &alpha, devBOLD_y, L2, devBOLD_y , L2, &beta, devCormat, Batch_size);
			if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;
			cudaStat = cudaMemcpy(out, devCormat, sizeof(real__t) * Batch_size * Batch_size, cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess) 
				return cudaStat;

			for (int jj = ii+1; jj < Num_Blocks; jj++)
			{
				//BOLD_t1 = BOLD_t + ii * Batch_size * L;
				//BOLD_t2 = BOLD_t + jj * Batch_size * L;
				out = Cormat + kk * Batch_size * Batch_size;
				kk++;
								
				pre_process <<<64, thread_num>>>(devBOLD_x, devBOLD_ori+jj * Batch_size * L, L, L2,Batch_size, tiecount, false);
				stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L2,  &alpha, devBOLD_x, L2, devBOLD_y , L2, &beta, devCormat, Batch_size);
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
		}
		cudaStat = cudaMemcpy(tie_count, tiecount, sizeof(real__t) * N, cudaMemcpyDeviceToHost);
		if (cudaStat != cudaSuccess) 
			return cudaStat;

		cudaFree (tiecount);
		cudaFree (devBOLD_x);
		cudaFree (devBOLD_y);
		cudaFree (devBOLD_ori);
		cudaFree (devCormat);
		stat = cublasDestroy(handle);
		if (stat != CUBLAS_STATUS_SUCCESS)
			return stat;
		delete []BOLD_t;
		return 1;
}




	
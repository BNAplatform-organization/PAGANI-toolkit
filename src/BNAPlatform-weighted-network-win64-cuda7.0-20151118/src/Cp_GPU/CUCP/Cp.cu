
#include "cuda_runtime.h"
#include "cusparse_v2.h"
#include "cublas_v2.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

const int n_thread = 256;
const int ny_thread = 16;
const int nx_thread = 16;

__global__ void cal_k_kernel(long N, int *d_R, float *d_CYC3 )
{		
	int i;
	const int blockid   = blockIdx.x;
	const int threadid = threadIdx.x;
	
	for(i = blockid*blockDim.x+threadid; i<N; i+=blockDim.x*gridDim.x )
	{
		float temp = (float) d_R[i+1];
		temp -= d_R[i];
		if (temp<2)
			d_CYC3[i] = 0;
		else 
			d_CYC3[i] = 1.0/(temp*(temp-1));
		//syncthreads();
	}
}

__global__ void pow3_vector(long N, float *d_vec )
{		
	int i;
	const int blockid   = blockIdx.x;
	const int threadid = threadIdx.x;
	
	for(i = blockid*blockDim.x+threadid; i<N; i+=blockDim.x*gridDim.x )
	{
		double temp = d_vec[i];
		temp = pow(temp,1.0/3);
		d_vec[i] = (float) temp;
		//syncthreads();
	}
}

__global__ void init_block (int *R , int*C, float *V, float *S,int N, int block_size)
{
	//__shared__ R_shared[ny_thread+1];
	int i = blockIdx.x*blockDim.y + threadIdx.y;
	int j = 0;
	//R_shared[threadIdx.y] = R[i]; 
	if (i<block_size)
	{
		float temp = 0;
		int startidx = R[i];
		int endidx = R[i+1];
		for (int k = startidx+threadIdx.x ; k < endidx; k += blockDim.x)
		{
			j = C[k];
			temp = V[k];
			S[i*N + j] = temp;
		}
	}
}

__global__ void calc_cci_block (int N, int size, int *R, int *C, float *V, float *S, float *cci)
{
	__shared__ float a[ny_thread][nx_thread];
	int i = blockIdx.x*blockDim.y + threadIdx.y;
	int j = 0;
	if (i<size)
	{
		float temp = 0;
		for (int k = R[i]+threadIdx.x ; k < R[i+1]; k += blockDim.x)
		{
			j = C[k];
			temp +=  V[k]*S[i*N + j];
		}
		a[threadIdx.y][threadIdx.x] = temp;
		syncthreads();

		for(j = ny_thread/2; j > 0; j/=2)
		{
			if (threadIdx.x < j)  a[threadIdx.y][threadIdx.x]+=a[threadIdx.y][threadIdx.x+j] ;
			syncthreads();
		}
		if(threadIdx.x==0)
			cci[i] = a[threadIdx.y][0];
	}
}



__global__ void dot_vv(long N, float alpha, float *d_a, float *d_b )
{		
	int i;
	const int blockid   = blockIdx.x;
	const int threadid = threadIdx.x;
	
	for(i = blockid*blockDim.x+threadid; i<N; i+=blockDim.x*gridDim.x )
	{
		float temp = d_a[i];
		temp = alpha * temp * d_b[i];
		d_a[i] = (float) temp;
		syncthreads();
	}
}

double Cp_2(int * C, int * R,float * V, float * Cp, int N)
{	

	int count[1];
	cudaError_t error;
	error = cudaGetDeviceCount(count); 
	if (error != cudaSuccess)
	{
		cerr<<"no CUDA device found."<<endl;
		return -1;
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
	

	int edge_num = R[N];
	cudaError_t cudaStat; //1,cudaStat2,cudaStat3,cudaStat4,cudaStat5,cudaStat6;
	cublasStatus_t stat_blas;
	cublasHandle_t handle_blas;
	cusparseStatus_t status;
	cusparseHandle_t handle=0;
	cusparseMatDescr_t descr=0;
	
	int * dev_C;
	int * dev_R;
	float *dev_V;
	//float *dev_Cp;
	float *dev_CYC3_INV;
	//float * dev_cyc3;
	float * dev_cci;
	
		
	cudaStat = cudaMalloc((void**)&dev_R, (N+1)*sizeof(int));
	if (cudaStat != cudaSuccess)
	{	cout << "dev_R malloc failed\n";
		return -1;
	}
	cudaStat = cudaMalloc((void**)&dev_C, R[N]*sizeof(int));
	if(cudaStat != cudaSuccess)
	{	cout << "dev_C malloc failed\n";
		return -1;
	}
	cudaStat = cudaMalloc((void**)&dev_V, R[N]*sizeof(int));
	if(cudaStat != cudaSuccess)
	{	cout << "dev_V malloc failed\n";
		return -1;
	}

	cudaStat = cudaMalloc((void**)&dev_CYC3_INV, N*sizeof(float));
	if(cudaStat != cudaSuccess)
	{	cout << "dev_CYC3_INV malloc failed\n";
		return -1;
	}

	/*cudaStat = cudaMalloc((void**)&dev_cyc3, N*sizeof(float));	
	if(cudaStat != cudaSuccess)
	{	cout << "dev_cyc3 alloc failed\n";
		return -1;
	}*/

	cudaStat = cudaMalloc((void**)&dev_cci, N*sizeof(float));
	if(cudaStat != cudaSuccess)
	{	cout << "dev_cci alloc failed\n";
		return -1;
	}

	cudaStat = cudaMemcpy(dev_R, R,(size_t)((N+1)*sizeof(R[0])),cudaMemcpyHostToDevice);
	if(cudaStat != cudaSuccess)
	{	cout << "dev_R memcpyh2d failed\n";
		return -1;
	}
	
	cudaStat = cudaMemcpy(dev_C, C,(size_t)(R[N]*sizeof(C[0])),cudaMemcpyHostToDevice);
	if(cudaStat != cudaSuccess)
	{	cout << "dev_C memcpyh2d failed\n";
		return -1;
	}

	cudaStat = cudaMemcpy(dev_V, V,(size_t)(R[N]*sizeof(V[0])),cudaMemcpyHostToDevice);
	if(cudaStat != cudaSuccess)
	{	cout << "dev_V memcpyh2d failed\n";
		return -1;
	}
	/*calc the inverse of CYC3*/

	

	cal_k_kernel<<<96,256>>>((long)N, dev_R, dev_CYC3_INV );

	stat_blas = cublasCreate(&handle_blas) ;
	if (stat_blas != CUBLAS_STATUS_SUCCESS)
	{	cout<<"blas failed\n";
		return stat_blas;
	}

	int Vm_idx;
	stat_blas = cublasIsamax( handle_blas, edge_num, dev_V, 1, &Vm_idx);
	if (stat_blas != CUBLAS_STATUS_SUCCESS)
	{	cout<<"blas failed\n";
		return stat_blas;
	}

	float alpha;
	cudaStat = cudaMemcpy(&alpha, dev_V+Vm_idx-1,(size_t) sizeof(float),cudaMemcpyDeviceToHost);
	if(cudaStat != cudaSuccess)
	{	cout << "V_max memcpyh2d failed\n";
		return -1;
	}
			
	cout<<"V_max = "<<alpha<<endl;
	
	alpha = 1.0/ (alpha);
	stat_blas = cublasSscal(handle_blas, edge_num, &alpha, dev_V, 1);
	if (stat_blas != CUBLAS_STATUS_SUCCESS)
	{	cout<<"scal failed\n";
		return stat_blas;
	}

	pow3_vector <<<96,256>>> (edge_num, dev_V);
	
	/* initialize cusparse library */
	status= cusparseCreate(&handle);
	if(status != CUSPARSE_STATUS_SUCCESS) 
	{	cout << "CUSPARSE Library initialization failed\n";
		return -1;
	}

	/* create and setup matrix descriptor */
	status= cusparseCreateMatDescr(&descr); 
	if(status != CUSPARSE_STATUS_SUCCESS) 
	{
		cout<<"Matrix descriptor initialization failed";
		return -1;
	} 
	cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO); 
		
	int block_size;
	block_size = GM_size/N/sizeof(float)/4;
	if (block_size > N)
		block_size = N;
	else if (block_size > 1024)
		block_size -= block_size%1024; 
	else if (block_size>512)
		block_size = 512;
	else if (block_size>32)
		block_size = block_size%32;
	else if (block_size>8)
		block_size = 8;
	else if (block_size>2)
		block_size = 2;
	else block_size = 1;
	
	block_size = 96*32;

	float alpha1 = 1.0;
	float beta = 0;
	float * dev_S;
	cudaStat = cudaMalloc((void**)&dev_S, sizeof(float)*N*block_size);
	if(cudaStat != cudaSuccess)
	{	cout << "dev_S alloc failed\n";
		return -1;
	}

	float * dev_SS;
	cudaStat = cudaMalloc((void**)&dev_SS, sizeof(float)*N*block_size);
	if(cudaStat != cudaSuccess)
	{	cout << "dev_SS alloc failed\n";
		return -1;
	}

	dim3 threadnum(ny_thread,nx_thread);
	int blocknum;
	
	int number_blocksize = (N + block_size - 1)/block_size;
	cout<<"block_size : "<<block_size<<"\nnumber of blocks : "<<number_blocksize<<endl;

	for (int i = 0; i < number_blocksize; i++ )
	{
		cudaMemset(dev_S, 0, sizeof(float)*N*block_size);
		cudaMemset(dev_SS, 0, sizeof(float)*N*block_size);
		
		//cout<<block_size;
		int size = block_size;
		if (i == number_blocksize - 1) size = N - i*block_size;
		blocknum = (size+threadnum.y-1)/threadnum.y;
		init_block<<<blocknum,threadnum>>> (dev_R + i*block_size, dev_C, dev_V, dev_S, N, size);
		
		//float * tmp_test = new float[N*block_size];
		//cudaMemcpy(tmp_test , dev_S, sizeof(float)*N*block_size, cudaMemcpyDeviceToHost);
		//cout<<endl<<endl;
		

		status = cusparseScsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
			N, size, N, R[N], &alpha1,
			descr, dev_V, dev_R, dev_C,
			dev_S, N, &beta, dev_SS, N);
		if(status != CUSPARSE_STATUS_SUCCESS) 
		{
			cout<<"Matrix mm failed";
			return -1;
		} 
		

		/*float * tmp_test1 = new float[N*block_size];
		cudaMemcpy(tmp_test1 , dev_SS, sizeof(float)*N*block_size, cudaMemcpyDeviceToHost);
		for (int k = 0; k < 100; k++)
			cout<<tmp_test1[k]<<endl;
		cout<<endl;
		for (int k = 4900; k < 5000; k++)
			cout<<tmp_test1[k]<<endl;
		delete[]tmp_test1;
		*/

		calc_cci_block<<<blocknum,threadnum>>>(N, size, dev_R + i*block_size, dev_C, dev_V, dev_SS, dev_cci + i*block_size);
		
		/*float * tmp_test2 = new float[size];
		cudaMemcpy(tmp_test2 , dev_cci + i*block_size, sizeof(float)*size, cudaMemcpyDeviceToHost);
		cout<<endl;
		for (int k = 0; k < 50; k++)
			cout<<tmp_test2[k]<<endl;
		cout<<endl;
		*/
		/*
		int baseS, nnzS;
		// nnzTotalDevHostPtr points to host memory
		int *dev_SR;
		int *dev_SC;
		float *dev_SV;

		cusparseMatDescr_t descrS = 0;
		status= cusparseCreateMatDescr(&descrS); 
		if(status != CUSPARSE_STATUS_SUCCESS) 
		{
		cout<<"Matrix S descriptor initialization failed";
		return -1;
		} 
		cusparseSetMatType(descrS,CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(descrS,CUSPARSE_INDEX_BASE_ZERO); 

		int *nnzTotalDevHostPtr = &nnzS;
		cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
		cudaMalloc((void**)&dev_SR, sizeof(int)*(N+1));

		status = cusparseXcsrgemmNnz(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, N, 
		descr, edge_num, dev_R, dev_C,
		descr, edge_num, dev_R, dev_C,
		descrS, dev_SR, nnzTotalDevHostPtr );
		if(status != CUSPARSE_STATUS_SUCCESS) 
		{
			cout<<"Matrix Xmm failed";
			return -1;
		} 

		if(NULL != nnzTotalDevHostPtr){
			nnzS = *nnzTotalDevHostPtr;
		}
		else{
			cudaMemcpy(&nnzS , dev_SR+N, sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(&baseS, dev_SR , sizeof(int), cudaMemcpyDeviceToHost);
			nnzS -= baseS;
		}
		cout<<nnzS<<endl;
		cudaStat = cudaMalloc((void**)&dev_SC, sizeof(int)*nnzS);
		if(cudaStat != cudaSuccess)
		{	cout << "dev_SC alloc failed\n";
			return -1;
		}
		cudaStat = cudaMalloc((void**)&dev_SV , sizeof(float)*nnzS);
		if(cudaStat != cudaSuccess)
		{	cout << "dev_SV alloc failed\n";
			return -1;
		}
		status = cusparseScsrgemm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, N,
		descr, edge_num,dev_V, dev_R, dev_C,
		descr, edge_num,dev_V, dev_R, dev_C,
		descrS, dev_SV, dev_SR, dev_SC);
		if(status != CUSPARSE_STATUS_SUCCESS) 
		{
			cout<<"Matrix Smm failed";
			return -1;
		} 
		
		int *SR = new int[N+1];
		cudaStat = cudaMemcpy(SR , dev_SR, (N+1)*sizeof(int), cudaMemcpyDeviceToHost);
		if(cudaStat != cudaSuccess)
		{	cout << "dev_SR memcpyd2h failed\n";
			return -1;
		}
		*/
	

	    /*
		for (int ii = i*block_size; ii < i*block_size+size; ii++)
		{	
			/*
			int idx = SR[i];
			int nnz = SR[i+1]-SR[i];
			cudaStat = cudaMemset(dev_cyc3, 0, N*sizeof(float));
			if(cudaStat != cudaSuccess)
			{	cout << "dev_cyc3 memset failed\n";
				return -1;
			}
			status = cusparseSsctr(handle,nnz, dev_SV+idx, dev_SC+idx, dev_cyc3, CUSPARSE_INDEX_BASE_ZERO);
			if(status != CUSPARSE_STATUS_SUCCESS) 
			{
				cout<<"Vector scatter failed";
				return -1;
			} 
			*/
			/*	
			int idx = R[ii];
			int nnz = R[ii+1]-R[ii];
			cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_DEVICE);
			status = cusparseSdoti(handle, nnz, dev_V+idx, dev_C+idx, dev_SS+(ii-i*block_size)*N, dev_cci+ii, CUSPARSE_INDEX_BASE_ZERO);
			if(status != CUSPARSE_STATUS_SUCCESS) 
			{
				cout<<"Vector doti failed";
				return -1;
			} 
			//cout<<endl<<tmp_test[i];
		}
		*/
	}
	
	
	/*float * tmp_test = new float[N];
	cudaMemcpy(tmp_test , dev_cci, sizeof(float)*N, cudaMemcpyDeviceToHost);
	for (int k = 0; k < 100; k++)
		cout<<tmp_test[k]<<endl;
	*/
	
	cudaFree(dev_S);
	cudaFree(dev_SS);
	dot_vv<<<96,256>>> (N, 1.0, dev_cci, dev_CYC3_INV);
	cudaStat = cudaMemcpy(Cp , dev_cci, N*sizeof(float), cudaMemcpyDeviceToHost);
	if(cudaStat != cudaSuccess)
	{	cout << "dev_cci memcpyd2h failed\n";
		return -1;
	}

	float mean_Cp = 0;
	cublasSasum( handle_blas, N, dev_cci, 1, &mean_Cp);
	
	mean_Cp = mean_Cp/N;
	
	//for (int k = 0; k < 50; k++)
		

	cout<<endl<<mean_Cp<<endl;
	
	//delete []SR;

	status = cusparseDestroy(handle);
	handle = 0;
	if(status != CUSPARSE_STATUS_SUCCESS) {
		cout<<"CUSPARSE Library release of resources failed\n";
		return -1;
	}

	stat_blas = cublasDestroy(handle_blas);
	if (stat_blas != CUBLAS_STATUS_SUCCESS)
		return -1;
	
	cudaFree (dev_C);
	cudaFree (dev_R);
	cudaFree (dev_V);
	cudaFree (dev_cci);
	cudaFree (dev_CYC3_INV);

	//cudaFree (dev_SR);
	//cudaFree (dev_SC);
	//cudaFree (dev_SV);
	//cudaFree (dev_cyc3);
	

	return ((double) mean_Cp);

}
		
	
	

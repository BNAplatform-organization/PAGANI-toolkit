# include <iostream>
# include "mytimer.h"
# include "FormBlock.h"
using namespace std;

# define SIZE_ROW_PER_THREAD_KATZ 6
# define SIZE_COL_PER_THREAD_KATZ 6
const int Katz_Block_Size = SIZE_ROW_PER_THREAD_KATZ * 16;

# define SHARED_FETCH_SIZE 16
__global__ void cuBFW_small_Katz_ph123(float * dst_ij, int k_block, int phase);
__global__ void cuBFW_small_Katz_ph2(float * dst_ij, int k_block);
__global__ void cuBFW_small_Katz_ph3(float * dst_ij, int k_block);

#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)
inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
	if(cudaSuccess != err)
	{
		fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString( err ) );
		exit(-1);        
	}
}

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

void BFW_CUDA_small_Katz(float *costmat, float * kernel_primary_block, int dim)
{
	cudaFuncSetCacheConfig(cuBFW_small_Katz_ph123, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(cuBFW_small_Katz_ph2, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(cuBFW_small_Katz_ph3, cudaFuncCachePreferShared);

	float * costmat_katz = new float[dim * dim];
	FormBlock(costmat_katz, costmat, dim, Katz_Block_Size);
	checkCudaErrors(cudaMemcpy(kernel_primary_block, costmat_katz, sizeof(float) * dim * dim, cudaMemcpyHostToDevice));

	int block_cnt = dim / Katz_Block_Size;

	dim3 dimBlock(Katz_Block_Size / SIZE_COL_PER_THREAD_KATZ, Katz_Block_Size / SIZE_ROW_PER_THREAD_KATZ); 
	dim3 dimGrid(block_cnt, block_cnt);

	STimer tmr;
	TimerInit(&tmr);
	TimerStart(&tmr);	
	for (int k_block = 0; k_block < block_cnt; k_block++)
	{
		cuBFW_small_Katz_ph123<<<dimGrid, dimBlock>>>(kernel_primary_block, k_block, 1);
		cuBFW_small_Katz_ph2<<<dimGrid, dimBlock>>>(kernel_primary_block, k_block);
		cuBFW_small_Katz_ph3<<<dimGrid, dimBlock>>>(kernel_primary_block, k_block);
	}
	getLastCudaError("Kernel execution failed");
/*	cudaThreadSynchronize();
	TimerStop(&tmr);

	double n1 = Katz_Block_Size / 1024.;
	double Gflop = 2*n1*n1*n1 * (block_cnt)* (block_cnt) * block_cnt;
	cout<<"Phase 1 time: "<<TimerGetRuntime(&tmr)*1000<<"ms, Achieving: "<< Gflop / TimerGetRuntime(&tmr)<<"Gflop/s"<<endl;	
*/
	checkCudaErrors(cudaMemcpy(costmat_katz, kernel_primary_block, sizeof(float) * dim * dim, cudaMemcpyDeviceToHost));
	DeFormBlock(costmat, costmat_katz, dim, Katz_Block_Size);
	delete []costmat_katz;
}

# define SHARED_FETCH_SIZE 16

__global__ void cuBFW_small_Katz_ph3(float * dst_ij, int k_block)
{
	if (blockIdx.y == k_block || blockIdx.x == k_block)
		return;
	float * src_ik = dst_ij + (blockIdx.y * gridDim.x + k_block) * Katz_Block_Size * Katz_Block_Size;
	float * src_kj = dst_ij + (k_block * gridDim.x + blockIdx.x) * Katz_Block_Size * Katz_Block_Size;
	dst_ij += (blockIdx.y * gridDim.x + blockIdx.x) * Katz_Block_Size * Katz_Block_Size;
	dst_ij += threadIdx.y * SIZE_ROW_PER_THREAD_KATZ * Katz_Block_Size + threadIdx.x * SIZE_COL_PER_THREAD_KATZ;
	int ID = threadIdx.y * blockDim.x + threadIdx.x;

	__shared__ float shared_src_kj1[Katz_Block_Size * SHARED_FETCH_SIZE]; 
	__shared__ float shared_src_ik1[Katz_Block_Size * SHARED_FETCH_SIZE]; 

	float * shared_src_ik_calc = shared_src_ik1 + threadIdx.y * SIZE_ROW_PER_THREAD_KATZ * SHARED_FETCH_SIZE;
	float * shared_src_kj_calc = shared_src_kj1 + threadIdx.x * SIZE_COL_PER_THREAD_KATZ;

	float reg_dst_ij_00, reg_dst_ij_01, reg_dst_ij_02, reg_dst_ij_03, reg_dst_ij_04, reg_dst_ij_05;
	float reg_dst_ij_10, reg_dst_ij_11, reg_dst_ij_12, reg_dst_ij_13, reg_dst_ij_14, reg_dst_ij_15;
	float reg_dst_ij_20, reg_dst_ij_21, reg_dst_ij_22, reg_dst_ij_23, reg_dst_ij_24, reg_dst_ij_25;
	float reg_dst_ij_30, reg_dst_ij_31, reg_dst_ij_32, reg_dst_ij_33, reg_dst_ij_34, reg_dst_ij_35;
	float reg_dst_ij_40, reg_dst_ij_41, reg_dst_ij_42, reg_dst_ij_43, reg_dst_ij_44, reg_dst_ij_45;
	float reg_dst_ij_50, reg_dst_ij_51, reg_dst_ij_52, reg_dst_ij_53, reg_dst_ij_54, reg_dst_ij_55;

	float reg_src_ik_0_, reg_src_ik_1_, reg_src_ik_2_, reg_src_ik_3_, reg_src_ik_4_, reg_src_ik_5_;
	float reg_src_kj__0, reg_src_kj__1, reg_src_kj__2, reg_src_kj__3, reg_src_kj__4, reg_src_kj__5;

	reg_dst_ij_00 = dst_ij[0 * Katz_Block_Size + 0];
	reg_dst_ij_01 = dst_ij[0 * Katz_Block_Size + 1];
	reg_dst_ij_02 = dst_ij[0 * Katz_Block_Size + 2];
	reg_dst_ij_03 = dst_ij[0 * Katz_Block_Size + 3];
	reg_dst_ij_04 = dst_ij[0 * Katz_Block_Size + 4];
	reg_dst_ij_05 = dst_ij[0 * Katz_Block_Size + 5];

	reg_dst_ij_10 = dst_ij[1 * Katz_Block_Size + 0];
	reg_dst_ij_11 = dst_ij[1 * Katz_Block_Size + 1];
	reg_dst_ij_12 = dst_ij[1 * Katz_Block_Size + 2];
	reg_dst_ij_13 = dst_ij[1 * Katz_Block_Size + 3];
	reg_dst_ij_14 = dst_ij[1 * Katz_Block_Size + 4];
	reg_dst_ij_15 = dst_ij[1 * Katz_Block_Size + 5];

	reg_dst_ij_20 = dst_ij[2 * Katz_Block_Size + 0];
	reg_dst_ij_21 = dst_ij[2 * Katz_Block_Size + 1];
	reg_dst_ij_22 = dst_ij[2 * Katz_Block_Size + 2];
	reg_dst_ij_23 = dst_ij[2 * Katz_Block_Size + 3];
	reg_dst_ij_24 = dst_ij[2 * Katz_Block_Size + 4];
	reg_dst_ij_25 = dst_ij[2 * Katz_Block_Size + 5];

	reg_dst_ij_30 = dst_ij[3 * Katz_Block_Size + 0];
	reg_dst_ij_31 = dst_ij[3 * Katz_Block_Size + 1];
	reg_dst_ij_32 = dst_ij[3 * Katz_Block_Size + 2];
	reg_dst_ij_33 = dst_ij[3 * Katz_Block_Size + 3];
	reg_dst_ij_34 = dst_ij[3 * Katz_Block_Size + 4];
	reg_dst_ij_35 = dst_ij[3 * Katz_Block_Size + 5];

	reg_dst_ij_40 = dst_ij[4 * Katz_Block_Size + 0];
	reg_dst_ij_41 = dst_ij[4 * Katz_Block_Size + 1];
	reg_dst_ij_42 = dst_ij[4 * Katz_Block_Size + 2];
	reg_dst_ij_43 = dst_ij[4 * Katz_Block_Size + 3];
	reg_dst_ij_44 = dst_ij[4 * Katz_Block_Size + 4];
	reg_dst_ij_45 = dst_ij[4 * Katz_Block_Size + 5];

	reg_dst_ij_50 = dst_ij[5 * Katz_Block_Size + 0];
	reg_dst_ij_51 = dst_ij[5 * Katz_Block_Size + 1];
	reg_dst_ij_52 = dst_ij[5 * Katz_Block_Size + 2];
	reg_dst_ij_53 = dst_ij[5 * Katz_Block_Size + 3];
	reg_dst_ij_54 = dst_ij[5 * Katz_Block_Size + 4];
	reg_dst_ij_55 = dst_ij[5 * Katz_Block_Size + 5];

//	# pragma unroll 2
	for (int kk = 0; kk < Katz_Block_Size; kk += SHARED_FETCH_SIZE)
	{
		__syncthreads();
		shared_src_ik1[ID + 0 * 256] = src_ik[kk + (threadIdx.y + 0 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 1 * 256] = src_ik[kk + (threadIdx.y + 1 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 2 * 256] = src_ik[kk + (threadIdx.y + 2 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 3 * 256] = src_ik[kk + (threadIdx.y + 3 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 4 * 256] = src_ik[kk + (threadIdx.y + 4 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 5 * 256] = src_ik[kk + (threadIdx.y + 5 * blockDim.y) * Katz_Block_Size + threadIdx.x];

		shared_src_kj1[ID + 0 * 256] = src_kj[kk * Katz_Block_Size + ID + 0 * 256];
		shared_src_kj1[ID + 1 * 256] = src_kj[kk * Katz_Block_Size + ID + 1 * 256];
		shared_src_kj1[ID + 2 * 256] = src_kj[kk * Katz_Block_Size + ID + 2 * 256];
		shared_src_kj1[ID + 3 * 256] = src_kj[kk * Katz_Block_Size + ID + 3 * 256];
		shared_src_kj1[ID + 4 * 256] = src_kj[kk * Katz_Block_Size + ID + 4 * 256];
		shared_src_kj1[ID + 5 * 256] = src_kj[kk * Katz_Block_Size + ID + 5 * 256];
		__syncthreads();

		//# pragma unroll 2
		for (int k = 0; k < SHARED_FETCH_SIZE; k++)
		{		
			reg_src_ik_0_ = shared_src_ik_calc[0 * SHARED_FETCH_SIZE + k];
			reg_src_ik_1_ = shared_src_ik_calc[1 * SHARED_FETCH_SIZE + k];
			reg_src_ik_2_ = shared_src_ik_calc[2 * SHARED_FETCH_SIZE + k];
			reg_src_ik_3_ = shared_src_ik_calc[3 * SHARED_FETCH_SIZE + k];
			reg_src_ik_4_ = shared_src_ik_calc[4 * SHARED_FETCH_SIZE + k];
			reg_src_ik_5_ = shared_src_ik_calc[5 * SHARED_FETCH_SIZE + k];
									 
			reg_src_kj__0 = shared_src_kj_calc[0 + k * Katz_Block_Size];
			reg_src_kj__1 = shared_src_kj_calc[1 + k * Katz_Block_Size];
			reg_src_kj__2 = shared_src_kj_calc[2 + k * Katz_Block_Size];
			reg_src_kj__3 = shared_src_kj_calc[3 + k * Katz_Block_Size];
			reg_src_kj__4 = shared_src_kj_calc[4 + k * Katz_Block_Size];
			reg_src_kj__5 = shared_src_kj_calc[5 + k * Katz_Block_Size];

			reg_dst_ij_00 = fminf(reg_dst_ij_00, reg_src_ik_0_ + reg_src_kj__0);
			reg_dst_ij_10 = fminf(reg_dst_ij_10, reg_src_ik_1_ + reg_src_kj__0);
			reg_dst_ij_20 = fminf(reg_dst_ij_20, reg_src_ik_2_ + reg_src_kj__0);
			reg_dst_ij_30 = fminf(reg_dst_ij_30, reg_src_ik_3_ + reg_src_kj__0);
			reg_dst_ij_40 = fminf(reg_dst_ij_40, reg_src_ik_4_ + reg_src_kj__0);
			reg_dst_ij_50 = fminf(reg_dst_ij_50, reg_src_ik_5_ + reg_src_kj__0);
																
			reg_dst_ij_01 = fminf(reg_dst_ij_01, reg_src_ik_0_ + reg_src_kj__1);
			reg_dst_ij_11 = fminf(reg_dst_ij_11, reg_src_ik_1_ + reg_src_kj__1);
			reg_dst_ij_21 = fminf(reg_dst_ij_21, reg_src_ik_2_ + reg_src_kj__1);
			reg_dst_ij_31 = fminf(reg_dst_ij_31, reg_src_ik_3_ + reg_src_kj__1);
			reg_dst_ij_41 = fminf(reg_dst_ij_41, reg_src_ik_4_ + reg_src_kj__1);
			reg_dst_ij_51 = fminf(reg_dst_ij_51, reg_src_ik_5_ + reg_src_kj__1);
																			
			reg_dst_ij_02 = fminf(reg_dst_ij_02, reg_src_ik_0_ + reg_src_kj__2);
			reg_dst_ij_12 = fminf(reg_dst_ij_12, reg_src_ik_1_ + reg_src_kj__2);
			reg_dst_ij_22 = fminf(reg_dst_ij_22, reg_src_ik_2_ + reg_src_kj__2);
			reg_dst_ij_32 = fminf(reg_dst_ij_32, reg_src_ik_3_ + reg_src_kj__2);
			reg_dst_ij_42 = fminf(reg_dst_ij_42, reg_src_ik_4_ + reg_src_kj__2);
			reg_dst_ij_52 = fminf(reg_dst_ij_52, reg_src_ik_5_ + reg_src_kj__2);
																			
			reg_dst_ij_03 = fminf(reg_dst_ij_03, reg_src_ik_0_ + reg_src_kj__3);
			reg_dst_ij_13 = fminf(reg_dst_ij_13, reg_src_ik_1_ + reg_src_kj__3);
			reg_dst_ij_23 = fminf(reg_dst_ij_23, reg_src_ik_2_ + reg_src_kj__3);
			reg_dst_ij_33 = fminf(reg_dst_ij_33, reg_src_ik_3_ + reg_src_kj__3);
			reg_dst_ij_43 = fminf(reg_dst_ij_43, reg_src_ik_4_ + reg_src_kj__3);
			reg_dst_ij_53 = fminf(reg_dst_ij_53, reg_src_ik_5_ + reg_src_kj__3);
																			
			reg_dst_ij_04 = fminf(reg_dst_ij_04, reg_src_ik_0_ + reg_src_kj__4);
			reg_dst_ij_14 = fminf(reg_dst_ij_14, reg_src_ik_1_ + reg_src_kj__4);
			reg_dst_ij_24 = fminf(reg_dst_ij_24, reg_src_ik_2_ + reg_src_kj__4);
			reg_dst_ij_34 = fminf(reg_dst_ij_34, reg_src_ik_3_ + reg_src_kj__4);
			reg_dst_ij_44 = fminf(reg_dst_ij_44, reg_src_ik_4_ + reg_src_kj__4);
			reg_dst_ij_54 = fminf(reg_dst_ij_54, reg_src_ik_5_ + reg_src_kj__4);
														   				
			reg_dst_ij_05 = fminf(reg_dst_ij_05, reg_src_ik_0_ + reg_src_kj__5);
			reg_dst_ij_15 = fminf(reg_dst_ij_15, reg_src_ik_1_ + reg_src_kj__5);
			reg_dst_ij_25 = fminf(reg_dst_ij_25, reg_src_ik_2_ + reg_src_kj__5);
			reg_dst_ij_35 = fminf(reg_dst_ij_35, reg_src_ik_3_ + reg_src_kj__5);
			reg_dst_ij_45 = fminf(reg_dst_ij_45, reg_src_ik_4_ + reg_src_kj__5);
			reg_dst_ij_55 = fminf(reg_dst_ij_55, reg_src_ik_5_ + reg_src_kj__5);
		}
	}

	dst_ij[0 * Katz_Block_Size + 0] = reg_dst_ij_00;
	dst_ij[0 * Katz_Block_Size + 1] = reg_dst_ij_01;
	dst_ij[0 * Katz_Block_Size + 2] = reg_dst_ij_02;
	dst_ij[0 * Katz_Block_Size + 3] = reg_dst_ij_03;
	dst_ij[0 * Katz_Block_Size + 4] = reg_dst_ij_04;
	dst_ij[0 * Katz_Block_Size + 5] = reg_dst_ij_05;

	dst_ij[1 * Katz_Block_Size + 0] = reg_dst_ij_10;
	dst_ij[1 * Katz_Block_Size + 1] = reg_dst_ij_11;
	dst_ij[1 * Katz_Block_Size + 2] = reg_dst_ij_12;
	dst_ij[1 * Katz_Block_Size + 3] = reg_dst_ij_13;
	dst_ij[1 * Katz_Block_Size + 4] = reg_dst_ij_14;
	dst_ij[1 * Katz_Block_Size + 5] = reg_dst_ij_15;

	dst_ij[2 * Katz_Block_Size + 0] = reg_dst_ij_20;
	dst_ij[2 * Katz_Block_Size + 1] = reg_dst_ij_21;
	dst_ij[2 * Katz_Block_Size + 2] = reg_dst_ij_22;
	dst_ij[2 * Katz_Block_Size + 3] = reg_dst_ij_23;
	dst_ij[2 * Katz_Block_Size + 4] = reg_dst_ij_24;
	dst_ij[2 * Katz_Block_Size + 5] = reg_dst_ij_25;

	dst_ij[3 * Katz_Block_Size + 0] = reg_dst_ij_30;
	dst_ij[3 * Katz_Block_Size + 1] = reg_dst_ij_31;
	dst_ij[3 * Katz_Block_Size + 2] = reg_dst_ij_32;
	dst_ij[3 * Katz_Block_Size + 3] = reg_dst_ij_33;
	dst_ij[3 * Katz_Block_Size + 4] = reg_dst_ij_34;
	dst_ij[3 * Katz_Block_Size + 5] = reg_dst_ij_35;

	dst_ij[4 * Katz_Block_Size + 0] = reg_dst_ij_40;
	dst_ij[4 * Katz_Block_Size + 1] = reg_dst_ij_41;
	dst_ij[4 * Katz_Block_Size + 2] = reg_dst_ij_42;
	dst_ij[4 * Katz_Block_Size + 3] = reg_dst_ij_43;
	dst_ij[4 * Katz_Block_Size + 4] = reg_dst_ij_44;
	dst_ij[4 * Katz_Block_Size + 5] = reg_dst_ij_45;

	dst_ij[5 * Katz_Block_Size + 0] = reg_dst_ij_50;
	dst_ij[5 * Katz_Block_Size + 1] = reg_dst_ij_51;
	dst_ij[5 * Katz_Block_Size + 2] = reg_dst_ij_52;
	dst_ij[5 * Katz_Block_Size + 3] = reg_dst_ij_53;
	dst_ij[5 * Katz_Block_Size + 4] = reg_dst_ij_54;
	dst_ij[5 * Katz_Block_Size + 5] = reg_dst_ij_55;
}

__global__ void cuBFW_small_Katz_ph2(float * dst_ij, int k_block)
{
	if ((blockIdx.y != k_block && blockIdx.x != k_block) || (blockIdx.y == k_block && blockIdx.x == k_block))
		return;
	float * src_ik = dst_ij + (blockIdx.y * gridDim.x + k_block) * Katz_Block_Size * Katz_Block_Size;
	float * src_kj = dst_ij + (k_block * gridDim.x + blockIdx.x) * Katz_Block_Size * Katz_Block_Size;
	dst_ij += (blockIdx.y * gridDim.x + blockIdx.x) * Katz_Block_Size * Katz_Block_Size;
	dst_ij += threadIdx.y * SIZE_ROW_PER_THREAD_KATZ * Katz_Block_Size + threadIdx.x * SIZE_COL_PER_THREAD_KATZ;
	int ID = threadIdx.y * blockDim.x + threadIdx.x;

	__shared__ float shared_src_kj1[Katz_Block_Size * SHARED_FETCH_SIZE]; 
	__shared__ float shared_src_ik1[Katz_Block_Size * SHARED_FETCH_SIZE]; 

	float * shared_src_ik_calc = shared_src_ik1 + threadIdx.y * SIZE_ROW_PER_THREAD_KATZ * SHARED_FETCH_SIZE;
	float * shared_src_kj_calc = shared_src_kj1 + threadIdx.x * SIZE_COL_PER_THREAD_KATZ;

	float reg_dst_ij_00, reg_dst_ij_01, reg_dst_ij_02, reg_dst_ij_03, reg_dst_ij_04, reg_dst_ij_05;
	float reg_dst_ij_10, reg_dst_ij_11, reg_dst_ij_12, reg_dst_ij_13, reg_dst_ij_14, reg_dst_ij_15;
	float reg_dst_ij_20, reg_dst_ij_21, reg_dst_ij_22, reg_dst_ij_23, reg_dst_ij_24, reg_dst_ij_25;
	float reg_dst_ij_30, reg_dst_ij_31, reg_dst_ij_32, reg_dst_ij_33, reg_dst_ij_34, reg_dst_ij_35;
	float reg_dst_ij_40, reg_dst_ij_41, reg_dst_ij_42, reg_dst_ij_43, reg_dst_ij_44, reg_dst_ij_45;
	float reg_dst_ij_50, reg_dst_ij_51, reg_dst_ij_52, reg_dst_ij_53, reg_dst_ij_54, reg_dst_ij_55;

	float reg_src_ik_0_, reg_src_ik_1_, reg_src_ik_2_, reg_src_ik_3_, reg_src_ik_4_, reg_src_ik_5_;
	float reg_src_kj__0, reg_src_kj__1, reg_src_kj__2, reg_src_kj__3, reg_src_kj__4, reg_src_kj__5;

	reg_dst_ij_00 = 1e10;//dst_ij[0 * Katz_Block_Size + 0];
	reg_dst_ij_01 = 1e10;//dst_ij[0 * Katz_Block_Size + 1];
	reg_dst_ij_02 = 1e10;//dst_ij[0 * Katz_Block_Size + 2];
	reg_dst_ij_03 = 1e10;//dst_ij[0 * Katz_Block_Size + 3];
	reg_dst_ij_04 = 1e10;//dst_ij[0 * Katz_Block_Size + 4];
	reg_dst_ij_05 = 1e10;//dst_ij[0 * Katz_Block_Size + 5];

	reg_dst_ij_10 = 1e10;//dst_ij[1 * Katz_Block_Size + 0];
	reg_dst_ij_11 = 1e10;//dst_ij[1 * Katz_Block_Size + 1];
	reg_dst_ij_12 = 1e10;//dst_ij[1 * Katz_Block_Size + 2];
	reg_dst_ij_13 = 1e10;//dst_ij[1 * Katz_Block_Size + 3];
	reg_dst_ij_14 = 1e10;//dst_ij[1 * Katz_Block_Size + 4];
	reg_dst_ij_15 = 1e10;//dst_ij[1 * Katz_Block_Size + 5];

	reg_dst_ij_20 = 1e10;//dst_ij[2 * Katz_Block_Size + 0];
	reg_dst_ij_21 = 1e10;//dst_ij[2 * Katz_Block_Size + 1];
	reg_dst_ij_22 = 1e10;//dst_ij[2 * Katz_Block_Size + 2];
	reg_dst_ij_23 = 1e10;//dst_ij[2 * Katz_Block_Size + 3];
	reg_dst_ij_24 = 1e10;//dst_ij[2 * Katz_Block_Size + 4];
	reg_dst_ij_25 = 1e10;//dst_ij[2 * Katz_Block_Size + 5];

	reg_dst_ij_30 = 1e10;//dst_ij[3 * Katz_Block_Size + 0];
	reg_dst_ij_31 = 1e10;//dst_ij[3 * Katz_Block_Size + 1];
	reg_dst_ij_32 = 1e10;//dst_ij[3 * Katz_Block_Size + 2];
	reg_dst_ij_33 = 1e10;//dst_ij[3 * Katz_Block_Size + 3];
	reg_dst_ij_34 = 1e10;//dst_ij[3 * Katz_Block_Size + 4];
	reg_dst_ij_35 = 1e10;//dst_ij[3 * Katz_Block_Size + 5];

	reg_dst_ij_40 = 1e10;//dst_ij[4 * Katz_Block_Size + 0];
	reg_dst_ij_41 = 1e10;//dst_ij[4 * Katz_Block_Size + 1];
	reg_dst_ij_42 = 1e10;//dst_ij[4 * Katz_Block_Size + 2];
	reg_dst_ij_43 = 1e10;//dst_ij[4 * Katz_Block_Size + 3];
	reg_dst_ij_44 = 1e10;//dst_ij[4 * Katz_Block_Size + 4];
	reg_dst_ij_45 = 1e10;//dst_ij[4 * Katz_Block_Size + 5];

	reg_dst_ij_50 = 1e10;//dst_ij[5 * Katz_Block_Size + 0];
	reg_dst_ij_51 = 1e10;//dst_ij[5 * Katz_Block_Size + 1];
	reg_dst_ij_52 = 1e10;//dst_ij[5 * Katz_Block_Size + 2];
	reg_dst_ij_53 = 1e10;//dst_ij[5 * Katz_Block_Size + 3];
	reg_dst_ij_54 = 1e10;//dst_ij[5 * Katz_Block_Size + 4];
	reg_dst_ij_55 = 1e10;//dst_ij[5 * Katz_Block_Size + 5];

//	# pragma unroll 2
	for (int kk = 0; kk < Katz_Block_Size; kk += SHARED_FETCH_SIZE)
	{
		__syncthreads();
		shared_src_ik1[ID + 0 * 256] = src_ik[kk + (threadIdx.y + 0 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 1 * 256] = src_ik[kk + (threadIdx.y + 1 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 2 * 256] = src_ik[kk + (threadIdx.y + 2 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 3 * 256] = src_ik[kk + (threadIdx.y + 3 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 4 * 256] = src_ik[kk + (threadIdx.y + 4 * blockDim.y) * Katz_Block_Size + threadIdx.x];
		shared_src_ik1[ID + 5 * 256] = src_ik[kk + (threadIdx.y + 5 * blockDim.y) * Katz_Block_Size + threadIdx.x];

		shared_src_kj1[ID + 0 * 256] = src_kj[kk * Katz_Block_Size + ID + 0 * 256];
		shared_src_kj1[ID + 1 * 256] = src_kj[kk * Katz_Block_Size + ID + 1 * 256];
		shared_src_kj1[ID + 2 * 256] = src_kj[kk * Katz_Block_Size + ID + 2 * 256];
		shared_src_kj1[ID + 3 * 256] = src_kj[kk * Katz_Block_Size + ID + 3 * 256];
		shared_src_kj1[ID + 4 * 256] = src_kj[kk * Katz_Block_Size + ID + 4 * 256];
		shared_src_kj1[ID + 5 * 256] = src_kj[kk * Katz_Block_Size + ID + 5 * 256];
		__syncthreads();

		//# pragma unroll 2
		for (int k = 0; k < SHARED_FETCH_SIZE; k++)
		{		
			reg_src_ik_0_ = shared_src_ik_calc[0 * SHARED_FETCH_SIZE + k];
			reg_src_ik_1_ = shared_src_ik_calc[1 * SHARED_FETCH_SIZE + k];
			reg_src_ik_2_ = shared_src_ik_calc[2 * SHARED_FETCH_SIZE + k];
			reg_src_ik_3_ = shared_src_ik_calc[3 * SHARED_FETCH_SIZE + k];
			reg_src_ik_4_ = shared_src_ik_calc[4 * SHARED_FETCH_SIZE + k];
			reg_src_ik_5_ = shared_src_ik_calc[5 * SHARED_FETCH_SIZE + k];
									 
			reg_src_kj__0 = shared_src_kj_calc[0 + k * Katz_Block_Size];
			reg_src_kj__1 = shared_src_kj_calc[1 + k * Katz_Block_Size];
			reg_src_kj__2 = shared_src_kj_calc[2 + k * Katz_Block_Size];
			reg_src_kj__3 = shared_src_kj_calc[3 + k * Katz_Block_Size];
			reg_src_kj__4 = shared_src_kj_calc[4 + k * Katz_Block_Size];
			reg_src_kj__5 = shared_src_kj_calc[5 + k * Katz_Block_Size];

			reg_dst_ij_00 = fminf(reg_dst_ij_00, reg_src_ik_0_ + reg_src_kj__0);
			reg_dst_ij_10 = fminf(reg_dst_ij_10, reg_src_ik_1_ + reg_src_kj__0);
			reg_dst_ij_20 = fminf(reg_dst_ij_20, reg_src_ik_2_ + reg_src_kj__0);
			reg_dst_ij_30 = fminf(reg_dst_ij_30, reg_src_ik_3_ + reg_src_kj__0);
			reg_dst_ij_40 = fminf(reg_dst_ij_40, reg_src_ik_4_ + reg_src_kj__0);
			reg_dst_ij_50 = fminf(reg_dst_ij_50, reg_src_ik_5_ + reg_src_kj__0);
																
			reg_dst_ij_01 = fminf(reg_dst_ij_01, reg_src_ik_0_ + reg_src_kj__1);
			reg_dst_ij_11 = fminf(reg_dst_ij_11, reg_src_ik_1_ + reg_src_kj__1);
			reg_dst_ij_21 = fminf(reg_dst_ij_21, reg_src_ik_2_ + reg_src_kj__1);
			reg_dst_ij_31 = fminf(reg_dst_ij_31, reg_src_ik_3_ + reg_src_kj__1);
			reg_dst_ij_41 = fminf(reg_dst_ij_41, reg_src_ik_4_ + reg_src_kj__1);
			reg_dst_ij_51 = fminf(reg_dst_ij_51, reg_src_ik_5_ + reg_src_kj__1);
																			
			reg_dst_ij_02 = fminf(reg_dst_ij_02, reg_src_ik_0_ + reg_src_kj__2);
			reg_dst_ij_12 = fminf(reg_dst_ij_12, reg_src_ik_1_ + reg_src_kj__2);
			reg_dst_ij_22 = fminf(reg_dst_ij_22, reg_src_ik_2_ + reg_src_kj__2);
			reg_dst_ij_32 = fminf(reg_dst_ij_32, reg_src_ik_3_ + reg_src_kj__2);
			reg_dst_ij_42 = fminf(reg_dst_ij_42, reg_src_ik_4_ + reg_src_kj__2);
			reg_dst_ij_52 = fminf(reg_dst_ij_52, reg_src_ik_5_ + reg_src_kj__2);
																			
			reg_dst_ij_03 = fminf(reg_dst_ij_03, reg_src_ik_0_ + reg_src_kj__3);
			reg_dst_ij_13 = fminf(reg_dst_ij_13, reg_src_ik_1_ + reg_src_kj__3);
			reg_dst_ij_23 = fminf(reg_dst_ij_23, reg_src_ik_2_ + reg_src_kj__3);
			reg_dst_ij_33 = fminf(reg_dst_ij_33, reg_src_ik_3_ + reg_src_kj__3);
			reg_dst_ij_43 = fminf(reg_dst_ij_43, reg_src_ik_4_ + reg_src_kj__3);
			reg_dst_ij_53 = fminf(reg_dst_ij_53, reg_src_ik_5_ + reg_src_kj__3);
																			
			reg_dst_ij_04 = fminf(reg_dst_ij_04, reg_src_ik_0_ + reg_src_kj__4);
			reg_dst_ij_14 = fminf(reg_dst_ij_14, reg_src_ik_1_ + reg_src_kj__4);
			reg_dst_ij_24 = fminf(reg_dst_ij_24, reg_src_ik_2_ + reg_src_kj__4);
			reg_dst_ij_34 = fminf(reg_dst_ij_34, reg_src_ik_3_ + reg_src_kj__4);
			reg_dst_ij_44 = fminf(reg_dst_ij_44, reg_src_ik_4_ + reg_src_kj__4);
			reg_dst_ij_54 = fminf(reg_dst_ij_54, reg_src_ik_5_ + reg_src_kj__4);
														   				
			reg_dst_ij_05 = fminf(reg_dst_ij_05, reg_src_ik_0_ + reg_src_kj__5);
			reg_dst_ij_15 = fminf(reg_dst_ij_15, reg_src_ik_1_ + reg_src_kj__5);
			reg_dst_ij_25 = fminf(reg_dst_ij_25, reg_src_ik_2_ + reg_src_kj__5);
			reg_dst_ij_35 = fminf(reg_dst_ij_35, reg_src_ik_3_ + reg_src_kj__5);
			reg_dst_ij_45 = fminf(reg_dst_ij_45, reg_src_ik_4_ + reg_src_kj__5);
			reg_dst_ij_55 = fminf(reg_dst_ij_55, reg_src_ik_5_ + reg_src_kj__5);
		}
	}

	dst_ij[0 * Katz_Block_Size + 0] = reg_dst_ij_00;
	dst_ij[0 * Katz_Block_Size + 1] = reg_dst_ij_01;
	dst_ij[0 * Katz_Block_Size + 2] = reg_dst_ij_02;
	dst_ij[0 * Katz_Block_Size + 3] = reg_dst_ij_03;
	dst_ij[0 * Katz_Block_Size + 4] = reg_dst_ij_04;
	dst_ij[0 * Katz_Block_Size + 5] = reg_dst_ij_05;

	dst_ij[1 * Katz_Block_Size + 0] = reg_dst_ij_10;
	dst_ij[1 * Katz_Block_Size + 1] = reg_dst_ij_11;
	dst_ij[1 * Katz_Block_Size + 2] = reg_dst_ij_12;
	dst_ij[1 * Katz_Block_Size + 3] = reg_dst_ij_13;
	dst_ij[1 * Katz_Block_Size + 4] = reg_dst_ij_14;
	dst_ij[1 * Katz_Block_Size + 5] = reg_dst_ij_15;

	dst_ij[2 * Katz_Block_Size + 0] = reg_dst_ij_20;
	dst_ij[2 * Katz_Block_Size + 1] = reg_dst_ij_21;
	dst_ij[2 * Katz_Block_Size + 2] = reg_dst_ij_22;
	dst_ij[2 * Katz_Block_Size + 3] = reg_dst_ij_23;
	dst_ij[2 * Katz_Block_Size + 4] = reg_dst_ij_24;
	dst_ij[2 * Katz_Block_Size + 5] = reg_dst_ij_25;

	dst_ij[3 * Katz_Block_Size + 0] = reg_dst_ij_30;
	dst_ij[3 * Katz_Block_Size + 1] = reg_dst_ij_31;
	dst_ij[3 * Katz_Block_Size + 2] = reg_dst_ij_32;
	dst_ij[3 * Katz_Block_Size + 3] = reg_dst_ij_33;
	dst_ij[3 * Katz_Block_Size + 4] = reg_dst_ij_34;
	dst_ij[3 * Katz_Block_Size + 5] = reg_dst_ij_35;

	dst_ij[4 * Katz_Block_Size + 0] = reg_dst_ij_40;
	dst_ij[4 * Katz_Block_Size + 1] = reg_dst_ij_41;
	dst_ij[4 * Katz_Block_Size + 2] = reg_dst_ij_42;
	dst_ij[4 * Katz_Block_Size + 3] = reg_dst_ij_43;
	dst_ij[4 * Katz_Block_Size + 4] = reg_dst_ij_44;
	dst_ij[4 * Katz_Block_Size + 5] = reg_dst_ij_45;

	dst_ij[5 * Katz_Block_Size + 0] = reg_dst_ij_50;
	dst_ij[5 * Katz_Block_Size + 1] = reg_dst_ij_51;
	dst_ij[5 * Katz_Block_Size + 2] = reg_dst_ij_52;
	dst_ij[5 * Katz_Block_Size + 3] = reg_dst_ij_53;
	dst_ij[5 * Katz_Block_Size + 4] = reg_dst_ij_54;
	dst_ij[5 * Katz_Block_Size + 5] = reg_dst_ij_55;
}

// no shared memory usage
__global__ void cuBFW_small_Katz_ph123(float * dst_ij, int k_block, int phase)
{
	if (phase == 3 && (blockIdx.y == k_block || blockIdx.x == k_block))
	{
		return;
	}
	if (phase == 2 && ((blockIdx.y != k_block && blockIdx.x != k_block) || (blockIdx.y == k_block && blockIdx.x == k_block)))
	{
		return;
	}
	if (phase == 1 && (blockIdx.y != k_block || blockIdx.x != k_block))
	{
		return;
	}

	float * src_ik = dst_ij + (blockIdx.y * gridDim.x + k_block) * Katz_Block_Size * Katz_Block_Size;
	float * src_kj = dst_ij + (k_block * gridDim.x + blockIdx.x) * Katz_Block_Size * Katz_Block_Size;
	dst_ij += (blockIdx.y * gridDim.x + blockIdx.x) * Katz_Block_Size * Katz_Block_Size;
	dst_ij += threadIdx.y * SIZE_ROW_PER_THREAD_KATZ * Katz_Block_Size + threadIdx.x * SIZE_COL_PER_THREAD_KATZ;
	src_kj += threadIdx.x * SIZE_COL_PER_THREAD_KATZ;
	src_ik += threadIdx.y * SIZE_ROW_PER_THREAD_KATZ * Katz_Block_Size;	

	float reg_dst_ij_00, reg_dst_ij_01, reg_dst_ij_02, reg_dst_ij_03, reg_dst_ij_04, reg_dst_ij_05;
	float reg_dst_ij_10, reg_dst_ij_11, reg_dst_ij_12, reg_dst_ij_13, reg_dst_ij_14, reg_dst_ij_15;
	float reg_dst_ij_20, reg_dst_ij_21, reg_dst_ij_22, reg_dst_ij_23, reg_dst_ij_24, reg_dst_ij_25;
	float reg_dst_ij_30, reg_dst_ij_31, reg_dst_ij_32, reg_dst_ij_33, reg_dst_ij_34, reg_dst_ij_35;
	float reg_dst_ij_40, reg_dst_ij_41, reg_dst_ij_42, reg_dst_ij_43, reg_dst_ij_44, reg_dst_ij_45;
	float reg_dst_ij_50, reg_dst_ij_51, reg_dst_ij_52, reg_dst_ij_53, reg_dst_ij_54, reg_dst_ij_55;

	float reg_src_ik_0_, reg_src_ik_1_, reg_src_ik_2_, reg_src_ik_3_, reg_src_ik_4_, reg_src_ik_5_;
	float reg_src_kj__0, reg_src_kj__1, reg_src_kj__2, reg_src_kj__3, reg_src_kj__4, reg_src_kj__5;

	reg_dst_ij_00 = dst_ij[0 * Katz_Block_Size + 0];
	reg_dst_ij_01 = dst_ij[0 * Katz_Block_Size + 1];
	reg_dst_ij_02 = dst_ij[0 * Katz_Block_Size + 2];
	reg_dst_ij_03 = dst_ij[0 * Katz_Block_Size + 3];
	reg_dst_ij_04 = dst_ij[0 * Katz_Block_Size + 4];
	reg_dst_ij_05 = dst_ij[0 * Katz_Block_Size + 5];

	reg_dst_ij_10 = dst_ij[1 * Katz_Block_Size + 0];
	reg_dst_ij_11 = dst_ij[1 * Katz_Block_Size + 1];
	reg_dst_ij_12 = dst_ij[1 * Katz_Block_Size + 2];
	reg_dst_ij_13 = dst_ij[1 * Katz_Block_Size + 3];
	reg_dst_ij_14 = dst_ij[1 * Katz_Block_Size + 4];
	reg_dst_ij_15 = dst_ij[1 * Katz_Block_Size + 5];

	reg_dst_ij_20 = dst_ij[2 * Katz_Block_Size + 0];
	reg_dst_ij_21 = dst_ij[2 * Katz_Block_Size + 1];
	reg_dst_ij_22 = dst_ij[2 * Katz_Block_Size + 2];
	reg_dst_ij_23 = dst_ij[2 * Katz_Block_Size + 3];
	reg_dst_ij_24 = dst_ij[2 * Katz_Block_Size + 4];
	reg_dst_ij_25 = dst_ij[2 * Katz_Block_Size + 5];

	reg_dst_ij_30 = dst_ij[3 * Katz_Block_Size + 0];
	reg_dst_ij_31 = dst_ij[3 * Katz_Block_Size + 1];
	reg_dst_ij_32 = dst_ij[3 * Katz_Block_Size + 2];
	reg_dst_ij_33 = dst_ij[3 * Katz_Block_Size + 3];
	reg_dst_ij_34 = dst_ij[3 * Katz_Block_Size + 4];
	reg_dst_ij_35 = dst_ij[3 * Katz_Block_Size + 5];

	reg_dst_ij_40 = dst_ij[4 * Katz_Block_Size + 0];
	reg_dst_ij_41 = dst_ij[4 * Katz_Block_Size + 1];
	reg_dst_ij_42 = dst_ij[4 * Katz_Block_Size + 2];
	reg_dst_ij_43 = dst_ij[4 * Katz_Block_Size + 3];
	reg_dst_ij_44 = dst_ij[4 * Katz_Block_Size + 4];
	reg_dst_ij_45 = dst_ij[4 * Katz_Block_Size + 5];

	reg_dst_ij_50 = dst_ij[5 * Katz_Block_Size + 0];
	reg_dst_ij_51 = dst_ij[5 * Katz_Block_Size + 1];
	reg_dst_ij_52 = dst_ij[5 * Katz_Block_Size + 2];
	reg_dst_ij_53 = dst_ij[5 * Katz_Block_Size + 3];
	reg_dst_ij_54 = dst_ij[5 * Katz_Block_Size + 4];
	reg_dst_ij_55 = dst_ij[5 * Katz_Block_Size + 5];

	k_block *= Katz_Block_Size;
	for (int k = 0; k < blockDim.x * SIZE_COL_PER_THREAD_KATZ; k++)
	{
		reg_src_ik_0_ = src_ik[0 * Katz_Block_Size + k];
		reg_src_ik_1_ = src_ik[1 * Katz_Block_Size + k];
		reg_src_ik_2_ = src_ik[2 * Katz_Block_Size + k];
		reg_src_ik_3_ = src_ik[3 * Katz_Block_Size + k];
		reg_src_ik_4_ = src_ik[4 * Katz_Block_Size + k];
		reg_src_ik_5_ = src_ik[5 * Katz_Block_Size + k];

		reg_src_kj__0 = src_kj[k * Katz_Block_Size + 0];
		reg_src_kj__1 = src_kj[k * Katz_Block_Size + 1];
		reg_src_kj__2 = src_kj[k * Katz_Block_Size + 2];
		reg_src_kj__3 = src_kj[k * Katz_Block_Size + 3];
		reg_src_kj__4 = src_kj[k * Katz_Block_Size + 4];
		reg_src_kj__5 = src_kj[k * Katz_Block_Size + 5];
		
		reg_dst_ij_00 = fminf(reg_dst_ij_00, reg_src_ik_0_ + reg_src_kj__0);
		reg_dst_ij_01 = fminf(reg_dst_ij_01, reg_src_ik_0_ + reg_src_kj__1);
		reg_dst_ij_02 = fminf(reg_dst_ij_02, reg_src_ik_0_ + reg_src_kj__2);
		reg_dst_ij_03 = fminf(reg_dst_ij_03, reg_src_ik_0_ + reg_src_kj__3);
		reg_dst_ij_04 = fminf(reg_dst_ij_04, reg_src_ik_0_ + reg_src_kj__4);
		reg_dst_ij_05 = fminf(reg_dst_ij_05, reg_src_ik_0_ + reg_src_kj__5);

		reg_dst_ij_10 = fminf(reg_dst_ij_10, reg_src_ik_1_ + reg_src_kj__0);
		reg_dst_ij_11 = fminf(reg_dst_ij_11, reg_src_ik_1_ + reg_src_kj__1);
		reg_dst_ij_12 = fminf(reg_dst_ij_12, reg_src_ik_1_ + reg_src_kj__2);
		reg_dst_ij_13 = fminf(reg_dst_ij_13, reg_src_ik_1_ + reg_src_kj__3);
		reg_dst_ij_14 = fminf(reg_dst_ij_14, reg_src_ik_1_ + reg_src_kj__4);
		reg_dst_ij_15 = fminf(reg_dst_ij_15, reg_src_ik_1_ + reg_src_kj__5);

		reg_dst_ij_20 = fminf(reg_dst_ij_20, reg_src_ik_2_ + reg_src_kj__0);
		reg_dst_ij_21 = fminf(reg_dst_ij_21, reg_src_ik_2_ + reg_src_kj__1);
		reg_dst_ij_22 = fminf(reg_dst_ij_22, reg_src_ik_2_ + reg_src_kj__2);
		reg_dst_ij_23 = fminf(reg_dst_ij_23, reg_src_ik_2_ + reg_src_kj__3);
		reg_dst_ij_24 = fminf(reg_dst_ij_24, reg_src_ik_2_ + reg_src_kj__4);
		reg_dst_ij_25 = fminf(reg_dst_ij_25, reg_src_ik_2_ + reg_src_kj__5);

		reg_dst_ij_30 = fminf(reg_dst_ij_30, reg_src_ik_3_ + reg_src_kj__0);
		reg_dst_ij_31 = fminf(reg_dst_ij_31, reg_src_ik_3_ + reg_src_kj__1);
		reg_dst_ij_32 = fminf(reg_dst_ij_32, reg_src_ik_3_ + reg_src_kj__2);
		reg_dst_ij_33 = fminf(reg_dst_ij_33, reg_src_ik_3_ + reg_src_kj__3);
		reg_dst_ij_34 = fminf(reg_dst_ij_34, reg_src_ik_3_ + reg_src_kj__4);
		reg_dst_ij_35 = fminf(reg_dst_ij_35, reg_src_ik_3_ + reg_src_kj__5);

		reg_dst_ij_40 = fminf(reg_dst_ij_40, reg_src_ik_4_ + reg_src_kj__0);
		reg_dst_ij_41 = fminf(reg_dst_ij_41, reg_src_ik_4_ + reg_src_kj__1);
		reg_dst_ij_42 = fminf(reg_dst_ij_42, reg_src_ik_4_ + reg_src_kj__2);
		reg_dst_ij_43 = fminf(reg_dst_ij_43, reg_src_ik_4_ + reg_src_kj__3);
		reg_dst_ij_44 = fminf(reg_dst_ij_44, reg_src_ik_4_ + reg_src_kj__4);
		reg_dst_ij_45 = fminf(reg_dst_ij_45, reg_src_ik_4_ + reg_src_kj__5);

		reg_dst_ij_50 = fminf(reg_dst_ij_50, reg_src_ik_5_ + reg_src_kj__0);
		reg_dst_ij_51 = fminf(reg_dst_ij_51, reg_src_ik_5_ + reg_src_kj__1);
		reg_dst_ij_52 = fminf(reg_dst_ij_52, reg_src_ik_5_ + reg_src_kj__2);
		reg_dst_ij_53 = fminf(reg_dst_ij_53, reg_src_ik_5_ + reg_src_kj__3);
		reg_dst_ij_54 = fminf(reg_dst_ij_54, reg_src_ik_5_ + reg_src_kj__4);
		reg_dst_ij_55 = fminf(reg_dst_ij_55, reg_src_ik_5_ + reg_src_kj__5);

		if (blockIdx.x * blockDim.x + threadIdx.x == (k_block + k + 1) / SIZE_COL_PER_THREAD_KATZ || blockIdx.y * blockDim.y + threadIdx.y == (k_block + k + 1) / SIZE_ROW_PER_THREAD_KATZ)
		//if (0)
		{
			dst_ij[0 * Katz_Block_Size + 0] = reg_dst_ij_00;
			dst_ij[0 * Katz_Block_Size + 1] = reg_dst_ij_01;
			dst_ij[0 * Katz_Block_Size + 2] = reg_dst_ij_02;
			dst_ij[0 * Katz_Block_Size + 3] = reg_dst_ij_03;
			dst_ij[0 * Katz_Block_Size + 4] = reg_dst_ij_04;
			dst_ij[0 * Katz_Block_Size + 5] = reg_dst_ij_05;

			dst_ij[1 * Katz_Block_Size + 0] = reg_dst_ij_10;
			dst_ij[1 * Katz_Block_Size + 1] = reg_dst_ij_11;
			dst_ij[1 * Katz_Block_Size + 2] = reg_dst_ij_12;
			dst_ij[1 * Katz_Block_Size + 3] = reg_dst_ij_13;
			dst_ij[1 * Katz_Block_Size + 4] = reg_dst_ij_14;
			dst_ij[1 * Katz_Block_Size + 5] = reg_dst_ij_15;

			dst_ij[2 * Katz_Block_Size + 0] = reg_dst_ij_20;
			dst_ij[2 * Katz_Block_Size + 1] = reg_dst_ij_21;
			dst_ij[2 * Katz_Block_Size + 2] = reg_dst_ij_22;
			dst_ij[2 * Katz_Block_Size + 3] = reg_dst_ij_23;
			dst_ij[2 * Katz_Block_Size + 4] = reg_dst_ij_24;
			dst_ij[2 * Katz_Block_Size + 5] = reg_dst_ij_25;

			dst_ij[3 * Katz_Block_Size + 0] = reg_dst_ij_30;
			dst_ij[3 * Katz_Block_Size + 1] = reg_dst_ij_31;
			dst_ij[3 * Katz_Block_Size + 2] = reg_dst_ij_32;
			dst_ij[3 * Katz_Block_Size + 3] = reg_dst_ij_33;
			dst_ij[3 * Katz_Block_Size + 4] = reg_dst_ij_34;
			dst_ij[3 * Katz_Block_Size + 5] = reg_dst_ij_35;

			dst_ij[4 * Katz_Block_Size + 0] = reg_dst_ij_40;
			dst_ij[4 * Katz_Block_Size + 1] = reg_dst_ij_41;
			dst_ij[4 * Katz_Block_Size + 2] = reg_dst_ij_42;
			dst_ij[4 * Katz_Block_Size + 3] = reg_dst_ij_43;
			dst_ij[4 * Katz_Block_Size + 4] = reg_dst_ij_44;
			dst_ij[4 * Katz_Block_Size + 5] = reg_dst_ij_45;

			dst_ij[5 * Katz_Block_Size + 0] = reg_dst_ij_50;
			dst_ij[5 * Katz_Block_Size + 1] = reg_dst_ij_51;
			dst_ij[5 * Katz_Block_Size + 2] = reg_dst_ij_52;
			dst_ij[5 * Katz_Block_Size + 3] = reg_dst_ij_53;
			dst_ij[5 * Katz_Block_Size + 4] = reg_dst_ij_54;
			dst_ij[5 * Katz_Block_Size + 5] = reg_dst_ij_55;
		}
		__syncthreads();	
	}

	dst_ij[0 * Katz_Block_Size + 0] = reg_dst_ij_00;
	dst_ij[0 * Katz_Block_Size + 1] = reg_dst_ij_01;
	dst_ij[0 * Katz_Block_Size + 2] = reg_dst_ij_02;
	dst_ij[0 * Katz_Block_Size + 3] = reg_dst_ij_03;
	dst_ij[0 * Katz_Block_Size + 4] = reg_dst_ij_04;
	dst_ij[0 * Katz_Block_Size + 5] = reg_dst_ij_05;

	dst_ij[1 * Katz_Block_Size + 0] = reg_dst_ij_10;
	dst_ij[1 * Katz_Block_Size + 1] = reg_dst_ij_11;
	dst_ij[1 * Katz_Block_Size + 2] = reg_dst_ij_12;
	dst_ij[1 * Katz_Block_Size + 3] = reg_dst_ij_13;
	dst_ij[1 * Katz_Block_Size + 4] = reg_dst_ij_14;
	dst_ij[1 * Katz_Block_Size + 5] = reg_dst_ij_15;

	dst_ij[2 * Katz_Block_Size + 0] = reg_dst_ij_20;
	dst_ij[2 * Katz_Block_Size + 1] = reg_dst_ij_21;
	dst_ij[2 * Katz_Block_Size + 2] = reg_dst_ij_22;
	dst_ij[2 * Katz_Block_Size + 3] = reg_dst_ij_23;
	dst_ij[2 * Katz_Block_Size + 4] = reg_dst_ij_24;
	dst_ij[2 * Katz_Block_Size + 5] = reg_dst_ij_25;

	dst_ij[3 * Katz_Block_Size + 0] = reg_dst_ij_30;
	dst_ij[3 * Katz_Block_Size + 1] = reg_dst_ij_31;
	dst_ij[3 * Katz_Block_Size + 2] = reg_dst_ij_32;
	dst_ij[3 * Katz_Block_Size + 3] = reg_dst_ij_33;
	dst_ij[3 * Katz_Block_Size + 4] = reg_dst_ij_34;
	dst_ij[3 * Katz_Block_Size + 5] = reg_dst_ij_35;

	dst_ij[4 * Katz_Block_Size + 0] = reg_dst_ij_40;
	dst_ij[4 * Katz_Block_Size + 1] = reg_dst_ij_41;
	dst_ij[4 * Katz_Block_Size + 2] = reg_dst_ij_42;
	dst_ij[4 * Katz_Block_Size + 3] = reg_dst_ij_43;
	dst_ij[4 * Katz_Block_Size + 4] = reg_dst_ij_44;
	dst_ij[4 * Katz_Block_Size + 5] = reg_dst_ij_45;

	dst_ij[5 * Katz_Block_Size + 0] = reg_dst_ij_50;
	dst_ij[5 * Katz_Block_Size + 1] = reg_dst_ij_51;
	dst_ij[5 * Katz_Block_Size + 2] = reg_dst_ij_52;
	dst_ij[5 * Katz_Block_Size + 3] = reg_dst_ij_53;
	dst_ij[5 * Katz_Block_Size + 4] = reg_dst_ij_54;
	dst_ij[5 * Katz_Block_Size + 5] = reg_dst_ij_55;
}
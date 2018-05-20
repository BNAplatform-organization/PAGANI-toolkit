# include <iostream>
# include "mytimer.h"
using namespace std;

# define SIZE_ROW_PER_THREAD_RX3 6
# define SIZE_COL_PER_THREAD_RX3 6
const int SHARED_BLOCK_SIZE = SIZE_ROW_PER_THREAD_RX3 * 16;

__global__ void cuBFW_ph2(float * dst_ij, float * src_ik, float * src_kj);

void BFW_CUDA_ph2(float * dst_ij, float * src_ik, float * src_kj, int block_size)
{
	cudaFuncSetCacheConfig(cuBFW_ph2, cudaFuncCachePreferShared);
	dim3 dimBlock(16, 16); 
	dim3 dimGrid(block_size / SHARED_BLOCK_SIZE, block_size / SHARED_BLOCK_SIZE);

	STimer tmr;
	TimerInit(&tmr);
	TimerStart(&tmr);	

	cuBFW_ph2<<<dimGrid, dimBlock>>>(dst_ij, src_ik, src_kj);

/*	cudaThreadSynchronize();
	TimerStop(&tmr);
	double n1 = block_size / 1024.;
	double Gflop = 2*n1*n1*n1;
	cout<<"Phase 2 time: "<<TimerGetRuntime(&tmr)*1000<<"ms, Achieving: "<< Gflop / TimerGetRuntime(&tmr)<<"Gflop/s"<<endl;	
	*/
}

# define SHARED_FETCH_SIZE 16

__global__ void cuBFW_ph2(float * dst_ij, float * src_ik, float * src_kj)
{
	
	int DIM = SIZE_COL_PER_THREAD_RX3 * blockDim.x * gridDim.x;
	dst_ij += blockIdx.y * SHARED_BLOCK_SIZE * DIM + blockIdx.x * SHARED_BLOCK_SIZE;
	dst_ij += threadIdx.y * SIZE_ROW_PER_THREAD_RX3 * DIM + threadIdx.x * SIZE_COL_PER_THREAD_RX3;
	src_ik += blockIdx.y * SHARED_BLOCK_SIZE * DIM;
	src_kj += blockIdx.x * SHARED_BLOCK_SIZE;

	__shared__ float shared_src_kj1[SHARED_BLOCK_SIZE * SHARED_FETCH_SIZE]; 
	__shared__ float shared_src_ik1[SHARED_BLOCK_SIZE * SHARED_FETCH_SIZE]; 

	float reg_dst_ij_00, reg_dst_ij_01, reg_dst_ij_02, reg_dst_ij_03, reg_dst_ij_04, reg_dst_ij_05;
	float reg_dst_ij_10, reg_dst_ij_11, reg_dst_ij_12, reg_dst_ij_13, reg_dst_ij_14, reg_dst_ij_15;
	float reg_dst_ij_20, reg_dst_ij_21, reg_dst_ij_22, reg_dst_ij_23, reg_dst_ij_24, reg_dst_ij_25;
	float reg_dst_ij_30, reg_dst_ij_31, reg_dst_ij_32, reg_dst_ij_33, reg_dst_ij_34, reg_dst_ij_35;
	float reg_dst_ij_40, reg_dst_ij_41, reg_dst_ij_42, reg_dst_ij_43, reg_dst_ij_44, reg_dst_ij_45;
	float reg_dst_ij_50, reg_dst_ij_51, reg_dst_ij_52, reg_dst_ij_53, reg_dst_ij_54, reg_dst_ij_55;

	float reg_src_ik_0_, reg_src_ik_1_, reg_src_ik_2_, reg_src_ik_3_, reg_src_ik_4_, reg_src_ik_5_;
	float reg_src_kj__0, reg_src_kj__1, reg_src_kj__2, reg_src_kj__3, reg_src_kj__4, reg_src_kj__5;

	reg_dst_ij_00 = 1e10;//dst_ij[0 * DIM + 0];
	reg_dst_ij_01 = 1e10;//dst_ij[0 * DIM + 1];
	reg_dst_ij_02 = 1e10;//dst_ij[0 * DIM + 2];
	reg_dst_ij_03 = 1e10;//dst_ij[0 * DIM + 3];
	reg_dst_ij_04 = 1e10;//dst_ij[0 * DIM + 4];
	reg_dst_ij_05 = 1e10;//dst_ij[0 * DIM + 5];

	reg_dst_ij_10 = 1e10;//dst_ij[1 * DIM + 0];
	reg_dst_ij_11 = 1e10;//dst_ij[1 * DIM + 1];
	reg_dst_ij_12 = 1e10;//dst_ij[1 * DIM + 2];
	reg_dst_ij_13 = 1e10;//dst_ij[1 * DIM + 3];
	reg_dst_ij_14 = 1e10;//dst_ij[1 * DIM + 4];
	reg_dst_ij_15 = 1e10;//dst_ij[1 * DIM + 5];

	reg_dst_ij_20 = 1e10;//dst_ij[2 * DIM + 0];
	reg_dst_ij_21 = 1e10;//dst_ij[2 * DIM + 1];
	reg_dst_ij_22 = 1e10;//dst_ij[2 * DIM + 2];
	reg_dst_ij_23 = 1e10;//dst_ij[2 * DIM + 3];
	reg_dst_ij_24 = 1e10;//dst_ij[2 * DIM + 4];
	reg_dst_ij_25 = 1e10;//dst_ij[2 * DIM + 5];

	reg_dst_ij_30 = 1e10;//dst_ij[3 * DIM + 0];
	reg_dst_ij_31 = 1e10;//dst_ij[3 * DIM + 1];
	reg_dst_ij_32 = 1e10;//dst_ij[3 * DIM + 2];
	reg_dst_ij_33 = 1e10;//dst_ij[3 * DIM + 3];
	reg_dst_ij_34 = 1e10;//dst_ij[3 * DIM + 4];
	reg_dst_ij_35 = 1e10;//dst_ij[3 * DIM + 5];

	reg_dst_ij_40 = 1e10;//dst_ij[4 * DIM + 0];
	reg_dst_ij_41 = 1e10;//dst_ij[4 * DIM + 1];
	reg_dst_ij_42 = 1e10;//dst_ij[4 * DIM + 2];
	reg_dst_ij_43 = 1e10;//dst_ij[4 * DIM + 3];
	reg_dst_ij_44 = 1e10;//dst_ij[4 * DIM + 4];
	reg_dst_ij_45 = 1e10;//dst_ij[4 * DIM + 5];

	reg_dst_ij_50 = 1e10;//dst_ij[5 * DIM + 0];
	reg_dst_ij_51 = 1e10;//dst_ij[5 * DIM + 1];
	reg_dst_ij_52 = 1e10;//dst_ij[5 * DIM + 2];
	reg_dst_ij_53 = 1e10;//dst_ij[5 * DIM + 3];
	reg_dst_ij_54 = 1e10;//dst_ij[5 * DIM + 4];
	reg_dst_ij_55 = 1e10;//dst_ij[5 * DIM + 5];

	//int ID = threadIdx.y * blockDim.x + threadIdx.x;
	float * shared_src_ik_calc = shared_src_ik1 + threadIdx.y * SIZE_ROW_PER_THREAD_RX3 * SHARED_FETCH_SIZE;
	float * shared_src_kj_calc = shared_src_kj1 + threadIdx.x * SIZE_COL_PER_THREAD_RX3;

//	# pragma unroll 2
	for (int kk = 0; kk < DIM; kk += SHARED_FETCH_SIZE)
	{
		__syncthreads();
		shared_src_ik1[threadIdx.y * SHARED_FETCH_SIZE + threadIdx.x + 0 * 256] = src_ik[kk + threadIdx.y * DIM + threadIdx.x + 0 * 16 * DIM];
		shared_src_ik1[threadIdx.y * SHARED_FETCH_SIZE + threadIdx.x + 1 * 256] = src_ik[kk + threadIdx.y * DIM + threadIdx.x + 1 * 16 * DIM];
		shared_src_ik1[threadIdx.y * SHARED_FETCH_SIZE + threadIdx.x + 2 * 256] = src_ik[kk + threadIdx.y * DIM + threadIdx.x + 2 * 16 * DIM];
		shared_src_ik1[threadIdx.y * SHARED_FETCH_SIZE + threadIdx.x + 3 * 256] = src_ik[kk + threadIdx.y * DIM + threadIdx.x + 3 * 16 * DIM];
		shared_src_ik1[threadIdx.y * SHARED_FETCH_SIZE + threadIdx.x + 4 * 256] = src_ik[kk + threadIdx.y * DIM + threadIdx.x + 4 * 16 * DIM];
		shared_src_ik1[threadIdx.y * SHARED_FETCH_SIZE + threadIdx.x + 5 * 256] = src_ik[kk + threadIdx.y * DIM + threadIdx.x + 5 * 16 * DIM];

		shared_src_kj1[threadIdx.y * SHARED_BLOCK_SIZE + threadIdx.x + 0 * 16] = src_kj[kk * DIM + threadIdx.y * DIM + threadIdx.x + 0 * 16];
		shared_src_kj1[threadIdx.y * SHARED_BLOCK_SIZE + threadIdx.x + 1 * 16] = src_kj[kk * DIM + threadIdx.y * DIM + threadIdx.x + 1 * 16];
		shared_src_kj1[threadIdx.y * SHARED_BLOCK_SIZE + threadIdx.x + 2 * 16] = src_kj[kk * DIM + threadIdx.y * DIM + threadIdx.x + 2 * 16];
		shared_src_kj1[threadIdx.y * SHARED_BLOCK_SIZE + threadIdx.x + 3 * 16] = src_kj[kk * DIM + threadIdx.y * DIM + threadIdx.x + 3 * 16];
		shared_src_kj1[threadIdx.y * SHARED_BLOCK_SIZE + threadIdx.x + 4 * 16] = src_kj[kk * DIM + threadIdx.y * DIM + threadIdx.x + 4 * 16];
		shared_src_kj1[threadIdx.y * SHARED_BLOCK_SIZE + threadIdx.x + 5 * 16] = src_kj[kk * DIM + threadIdx.y * DIM + threadIdx.x + 5 * 16];

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
									 
			reg_src_kj__0 = shared_src_kj_calc[0 + k * SHARED_BLOCK_SIZE];
			reg_src_kj__1 = shared_src_kj_calc[1 + k * SHARED_BLOCK_SIZE];
			reg_src_kj__2 = shared_src_kj_calc[2 + k * SHARED_BLOCK_SIZE];
			reg_src_kj__3 = shared_src_kj_calc[3 + k * SHARED_BLOCK_SIZE];
			reg_src_kj__4 = shared_src_kj_calc[4 + k * SHARED_BLOCK_SIZE];
			reg_src_kj__5 = shared_src_kj_calc[5 + k * SHARED_BLOCK_SIZE];

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

	dst_ij[0 * DIM + 0] = reg_dst_ij_00;
	dst_ij[0 * DIM + 1] = reg_dst_ij_01;
	dst_ij[0 * DIM + 2] = reg_dst_ij_02;
	dst_ij[0 * DIM + 3] = reg_dst_ij_03;
	dst_ij[0 * DIM + 4] = reg_dst_ij_04;
	dst_ij[0 * DIM + 5] = reg_dst_ij_05;

	dst_ij[1 * DIM + 0] = reg_dst_ij_10;
	dst_ij[1 * DIM + 1] = reg_dst_ij_11;
	dst_ij[1 * DIM + 2] = reg_dst_ij_12;
	dst_ij[1 * DIM + 3] = reg_dst_ij_13;
	dst_ij[1 * DIM + 4] = reg_dst_ij_14;
	dst_ij[1 * DIM + 5] = reg_dst_ij_15;

	dst_ij[2 * DIM + 0] = reg_dst_ij_20;
	dst_ij[2 * DIM + 1] = reg_dst_ij_21;
	dst_ij[2 * DIM + 2] = reg_dst_ij_22;
	dst_ij[2 * DIM + 3] = reg_dst_ij_23;
	dst_ij[2 * DIM + 4] = reg_dst_ij_24;
	dst_ij[2 * DIM + 5] = reg_dst_ij_25;

	dst_ij[3 * DIM + 0] = reg_dst_ij_30;
	dst_ij[3 * DIM + 1] = reg_dst_ij_31;
	dst_ij[3 * DIM + 2] = reg_dst_ij_32;
	dst_ij[3 * DIM + 3] = reg_dst_ij_33;
	dst_ij[3 * DIM + 4] = reg_dst_ij_34;
	dst_ij[3 * DIM + 5] = reg_dst_ij_35;

	dst_ij[4 * DIM + 0] = reg_dst_ij_40;
	dst_ij[4 * DIM + 1] = reg_dst_ij_41;
	dst_ij[4 * DIM + 2] = reg_dst_ij_42;
	dst_ij[4 * DIM + 3] = reg_dst_ij_43;
	dst_ij[4 * DIM + 4] = reg_dst_ij_44;
	dst_ij[4 * DIM + 5] = reg_dst_ij_45;

	dst_ij[5 * DIM + 0] = reg_dst_ij_50;
	dst_ij[5 * DIM + 1] = reg_dst_ij_51;
	dst_ij[5 * DIM + 2] = reg_dst_ij_52;
	dst_ij[5 * DIM + 3] = reg_dst_ij_53;
	dst_ij[5 * DIM + 4] = reg_dst_ij_54;
	dst_ij[5 * DIM + 5] = reg_dst_ij_55;

}
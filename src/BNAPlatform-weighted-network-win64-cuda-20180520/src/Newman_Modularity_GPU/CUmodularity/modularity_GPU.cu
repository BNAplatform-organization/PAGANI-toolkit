#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <memory.h>
#include <fstream>
#include <cstring>
#include "dirent.h" 
#include "device_functions.h"
#include "modularity_GPU.cuh"
#include <cmath>
#include <time.h> 
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/fill.h>
#include <thrust/extrema.h>
#include <math_constants.h>
//#include <cublas_v2.h>
//#include <cusparse.h>
using namespace std;
//
//void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);
//
//#define RANDOM_V0
//
//extern long long N, Ntemp;
//extern double * v, *vv; 
//extern double * v0, * verr; 
//extern	double * sumBG;
//extern long long seed;
//
typedef  unsigned int u_int;
// set some parameters

extern const int MAX_ITER;			// The maximum iteration times in the power method
extern const double BETA_Adjust;		// An optional parameter for quicker convergence. Its effect is uncertain
extern const double Epsilon;	// If |x - x0| < Epsilon, quit iteraion 
extern const double LAMBDA;		// if labmda > LAMBDA, initiate the division
extern const double DQ_MIN;
extern const int MIN_GROUP;			// The minimum nodes of an allowed module 

extern ofstream fout;	

const int threadnumx = 16;
const int threadnumy = 16;
const int  threadnum = 256;
const int blocknum    = 96;
//int* d_AD_init;
//int* d_AD;
//int* d_r; 
//int* d_c;
//int* d_orir;
//double* d_u;
//double* d_u0;
//double* d_uu;
//double * d_sumBG;
//double * d_norm;
////bool* d_G;
//double * temp_result;
//double * d_vector;
//double * d_vector1;
//double * d_vector2;
//double *d_k;
//double *d_orik;
//
//void Partition(int * R, int * C, int * Result);
//bool Sub_Partition(int * OriR, int * OriC, int * R, int * C, int M, long long innerM, int * Result, int * Max_Result,int * AD);
//double Lead_Vector(int * OriR, int * R, int * C, int M, double * sumBG1, double beta, int *AD, double *v, double *vv);
//

cusparseHandle_t s_handle=0;
cusparseMatDescr_t s_descr=0;
cublasHandle_t handle;

double Lead_Vector(int N, int * R, int * C, float *V, double *K,  double * sumBG, double m, double beta, double *u);
double calculate_dQ(signed char *S, int N, int * R, int * C, float *V, double *K, double *sumBG, double m);
double fine_tune_S(double dQ, int N, int * R, int * C, float * V, double * K, double m, signed char *S);
double qmax(int N, double *Qit, char *indSub, int *imax);
template <class Type> double VectorNorm(Type * x, int N);

__global__ void init_du(int N, double *d_u)
{
	int tid = blockDim.x*blockIdx.x+threadIdx.x;
	int i;
	for (i = tid; i < N; i+=blockDim.x*gridDim.x)
		d_u[i] = (1.0 * i) / (N-1);
}

__global__ void sign_eigv(int N, double *d_eigv)
{
	int tid = blockDim.x*blockIdx.x+threadIdx.x;
	int i;
	double temp = 0;
	for (i = tid; i<N; i+=blockDim.x*gridDim.x)
	{	
		temp = ((d_eigv[i] >= 0) ? 1.0 : -1.0);
		syncthreads();
		d_eigv[i] = temp;
	}
}

//
//__global__ void cal_k_kernel(long N, int *AD,int *d_OriR, double *d_K )
//{		
//	int offset;
//	const int blockid   = blockIdx.x;
//	const int threadid = threadIdx.x;
//	
//	for(offset=(blockid/2)*threadnum*2+threadid*2+blockid%2; offset<N; offset+=blockDim.x*gridDim.x )
//		d_K[offset]= (double)(d_OriR[AD[offset]+1]-d_OriR[AD[offset]]);
//}
//
//__global__ void sum_kj_kernel(long N, double *d_K ,double *sum_kj)
//{	__shared__ int sum[threadnum];
//	int temp=0;
//	int offset;
//   const int threadid =blockIdx.x*blockDim.x + threadIdx.x;
//	
//	//for(offset=(blockid/2)*threadnum*2+threadid*2+blockid%2; offset<N; offset+=blockDim.x*gridDim.x )
//		// temp+= dG[offset]*(d_R[offset+1]-d_R[offset]);
//	//for (offset=blockid*threadnum; offset+threadid<N; offset+=blockDim.x*gridDim.x ) 
//	 	//  temp +=dG[offset+threadid]*(d_R[offset + threadid+1]-d_R[offset + threadid]);	 		
//	
//	for(offset=threadid; offset<N; offset+=blockDim.x*gridDim.x)
//	{
//		temp+=(int)d_K[offset] ;
//	}
//	sum[threadIdx.x]=temp;
//	syncthreads();
//	
//	for(offset=1;offset+threadIdx.x<threadnum;offset*=2){
//			if (threadIdx.x%(2*offset)==0)  sum[threadIdx.x]+=sum[threadIdx.x+offset] ;
//			syncthreads();
//	}
//	if(threadIdx.x==0)
//	  sum_kj[blockIdx.x]=(double) sum[threadIdx.x];
//}
//
//__global__ void sum_kernel(int size, double *data, double scale)		   //将一些前面操作中用树形结构相加每个block的结果求和。
//{	
//	__shared__ double sum[blocknum];
//	int offset;
//	sum[threadIdx.x]=data[threadIdx.x];
//	syncthreads();
//	
//	for(offset=1;offset+threadIdx.x<blocknum; offset*=2){
//		if (threadIdx.x%(2*offset)==0)  sum[threadIdx.x]+=sum[threadIdx.x+offset] ;
//		syncthreads();
//	}
//	if(threadIdx.x==0)
//	    data[0]=sum[0]*scale; 
//}
//
//__global__ void spmv_one_thread(long N, long M, double *result, int *R, int *C, double *vv, double *dk, double vk, double *d_sum, double beta, double *v0)	   //计算Ai*v0,   每16个threads计算一行
//{
//	__shared__ int R_shared[threadnum+1];
//
//	double temp1=0;
//	int offset;  
//	
//	for (offset=blockIdx.x*threadnum+threadIdx.x;offset<N; offset+=gridDim.x*blockDim.x)
//	{
//		R_shared[threadIdx.x+1] = R[offset+1];
//		if(threadIdx.x==0) R_shared[threadIdx.x]=R[offset] ;
//		syncthreads();	 
//		temp1 =0 ;
//		for(int i=R_shared[threadIdx.x]; i<R_shared[threadIdx.x+1];i++)
//		{		
//			temp1+=vv[C[i]];			
//		}		
//		temp1-=vk/(2*M)*dk[offset]+(d_sum[offset]-beta)*v0[offset];
//		result[offset]=temp1 ;
//			//syncthreads();
//	}		
//}
//
__global__ void vvdot_subtr (long N, double *result_vector, double *vx1, double *vx2, double *vy)    //计算向量加法
 {
	 const int blockid   = blockIdx.x;
	 const int threadid  = threadIdx.x;
	 int offset;  	 
	 for(offset=threadid+blockid*blockDim.x; offset < N; offset+=blockDim.x*gridDim.x)
		 result_vector[offset] = vx1[offset] * vx2[offset] - vy[offset];
 }

__global__ void calculate_du (int N, double alpha, double beta, double *d_u, double *d_K, double *d_sumBG, double *d_u0)
{
	int tid = blockDim.x*blockIdx.x+threadIdx.x;
	int i;
	double temp;
	// d_u = d_u - alpha*K - (sumBG+beta).*u0
	for (i = tid; i<N; i+=blockDim.x*gridDim.x)
	{
		temp = d_u[i];
		syncthreads();
		temp -= alpha * d_K[i];
		syncthreads();
		temp -= (d_sumBG[i]+beta)*d_u0[i];		
		syncthreads();
		
		d_u[i] = temp ;
	}
}

__global__ void calculate_Qit (int N, double Q, double *d_Qit, double *d_temp_vector, double *d_K, double *d_S, double m, char *d_indSub)
{
	int tid = blockDim.x*blockIdx.x+threadIdx.x;
	double temp = 0, ki = 0, si = 0;
	char f = 0;
	// 
	for (int i = tid; i<N; i+=blockDim.x*gridDim.x)
	{
		f = d_indSub[i];
		if (f)
		{
			temp = d_temp_vector[i];
			//syncthreads();
			ki = d_K[i];
			//syncthreads();
			si = d_S[i];
			temp = temp*si +  ki*ki/m;
			//syncthreads();
			temp = Q - 4 * temp ;
			d_Qit[i] = (temp >= 0 ? temp : 0);
			//d_Qit[i] = temp;
		}
		else 
			d_Qit[i] = 0;

	}
}

__global__ void update_array(long imax, double *d_Sit)
{
	d_Sit[imax] = d_Sit[imax];
	//indSub[imax] = CUDART_NAN_F; //at risk
}
// 
// __global__ void sumBG_kernel (double *sumBG, double *d_orik, double *d_k, int *AD,long N, double innerMd2M)    //计算向量加法
// {
//	 //sumBG[i] = R[i+1] - R[i] - (OriR[AD[i]+1] - OriR[AD[i]]) * (double)innerM / 2 / M;
//	 const int blockid   = blockIdx.x;
//	 const int threadid  = threadIdx.x;
//	 int offset;  	 
//	 for(offset=threadid+blockid*threadnum; offset<N; offset+=threadnum*gridDim.x)
//	 {
//		sumBG[offset] = d_k[offset] - (d_orik[AD[offset]]) * innerMd2M;
//	 }
// }
// 
// __global__ void cal_VV_kernel ( long Ntemp, long N, int *AD, double *vv, double *v)
// {
//	 const int blockid   = blockIdx.x;
//	 const int threadid  = threadIdx.x;
//	 int offset;
//	 for(offset=threadid+blockid*threadnum; offset<Ntemp; offset+=threadnum*gridDim.x)
//		 vv[AD[offset]] = v[offset];
// }
//
// __global__ void calc_vector (long N, double *result, double *dk, double vk2m, double *d_sum, double beta, double *v0)	   //计算 ( v_k/2/M * (R[i+1] - R[i])+(sumBG1[i] - beta) * v0[i])，
// {		
//	 const int blockid   = blockIdx.x;
//	 const int threadid = threadIdx.x;
//	 int offset;
//	 double temp=0;
//	 
//	 for(offset=threadnum*blockid+threadid; offset<N; offset+=threadnum*gridDim.x)
//	 {	
//		  temp=vk2m*dk[offset]+(d_sum[offset]-beta)*v0[offset];
//		  result[offset]=temp;
//	 }
//  }
//
// /*__global__ void  Norm2_ph1(long N, double *norm,  double *v, bool *dG)	 //求向量的二范数，配合sum_kernel 得到最终结果
// {
//	 __shared__  double temp[threadnum];
//	 const int blockid   = blockIdx.x;
//	 const int threadid = threadIdx.x;
//	 double temp1=0;
//	 int offset;
//
//	 for(offset=blockid*threadnum+threadid; offset<N; offset+=threadnum*gridDim.x)
//		 temp1+=dG[offset]? v[offset]*v[offset] : 0;
//	 temp[threadid]=temp1;
//	 syncthreads();
//
//	 for(offset=1;offset+threadid<threadnum;offset*=2){
//		 if (threadid%(2*offset)==0)  temp[threadid]+=temp[threadid+offset] ;
//		 syncthreads();
//	 }
//	 if (threadid==0)
//		 norm[blockid] = temp[threadid];  	 
// }
// __global__ void Norm2_ph2(long N,  double norm, double *v)						  //将向量归一化，若为零向量，则每个元素除以1，保持不变。
// {	 	 
//	 for (int offset = threadIdx.x+blockIdx.x*threadnum; offset<N ; offset+= threadnum*gridDim.x  )
//		 v[offset]/=(norm? (norm) : 1);  
// } */
//
//
//
//
//
//
//
///* 
//This function returns the norm of the input vector x[G].
//G is the logic subscriber and N is the matrix dimension.
//*/
//
//
void cudaDevice_check()
{
	int devID;
	cudaDeviceProp deviceProps;
	devID = findCudaDevice();
	// get number of SMs on this GPU
	checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
}

double Lead_Vector_GPU(int N, int nnz, int *d_R, int *d_C, double *d_V, double *d_K, double *d_sumBG, double m, double beta, double *d_u)
{
	int i = 0, j = 0;
	
	// Initialize d_u. Two methods are optional. Define RANDOM_V0 if you want to use random starting vector
	init_du<<<blocknum,threadnum>>>(N, d_u);
//#ifdef RANDOM_V0
//	srand(time(0));
//	//srand(2016);
//	u[0] = 1;                        //normalized;
//	for (i = 1; i < N; i++)
//		u[i] = rand()*1.0f/RAND_MAX;
//	
//#else
//	for (i = 0; i < N; i++)
//		u[i] = 1.0f*i/(N-1);
//#endif
	
	double *d_u0;
	checkCudaErrors( cudaMalloc( (void**) &d_u0, sizeof(double) * (N)));
	//double err1 = 1, err2 = 1;
	double err = 1;
	u_int ITER = 0;
	double uNorm = 0, uNorm0 = 0;
	double u_k;
	double temp1;
	double alpha = 1.0;
	double cublasbeta = 0;
	double *K = new double [N];
	double *sumBG = new double [N];
	checkCudaErrors( cudaMemcpy(K, d_K, sizeof(double)*N, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(sumBG, d_sumBG, sizeof(double)*N, cudaMemcpyDeviceToHost));
	
	//double *u = new double [N];
	//double *u0 = new double [N];

	while (err > Epsilon && ITER < MAX_ITER)
	{
		//for (i = 0; i < N; i++)
		//	v0[i] =  v[i];
		//checkCudaErrors( cudaMemcpy(d_u, u,sizeof(double)*N, cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(d_u0,d_u,sizeof(double)*N, cudaMemcpyDeviceToDevice));
		
		// The dot product of u_k = k'*u, 
		u_k = 0;
		checkCublasErrors( cublasDdot (handle, N, d_u0, 1, d_K, 1, &u_k), "Ddot err in  u*k"); 
		//for (i = 0; i < N; i++)
		//	u_k += u0[i] * (K[i]);
		//cout<<"u_k = "<<u_k<<endl;
		
		// Do the matrix-vector multiplication
		alpha = 1.0;
		checkCusparseErrors( cusparseDcsrmv(s_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nnz, &alpha, s_descr, d_V, d_R, d_C, d_u0, &cublasbeta, d_u),
		"Dcsrmv err in calculate lead_vector");
		
		//checkCublasErrors (cublasDaxpy(handle, N, &alpha, d_K, 1, d_u, 1), "Daxpy err in calculate lead_vector");
		
		
		//checkCudaErrors( cudaMemcpy(u,d_u,sizeof(double)*N, cudaMemcpyDeviceToHost));
		alpha = u_k/m;
		calculate_du<<<blocknum,threadnum>>>(N, alpha, beta, d_u, d_K, d_sumBG, d_u0);
		
		
		/*for (i = 0; i < N; i++)
		{
			temp1 = 0;
			for (j = R[i]; j < R[i+1]; j++)
				temp1 += V[j] * u0[C[j]];
			temp1 -= u_k / m * (K[i]) + (sumBG[i] + beta) * u0[i];
			u[i] = temp1;
		}*/
		
		//ofstream debugfile;
		//debugfile.open("debug_u",ios::binary|ios::out);
		//debugfile.write((char *)u, N*sizeof(double));
		//debugfile.close();
		
		int idx;
		//double *x = new double [N];
		//checkCudaErrors( cudaMemcpy(d_u, u, sizeof(double)*N, cudaMemcpyHostToDevice));
		checkCublasErrors(cublasIdamax(handle, N, d_u, 1, &idx),"amax err in calculating lead_vector!\n");
		checkCudaErrors( cudaMemcpy(&uNorm,d_u+idx-1,sizeof(double), cudaMemcpyDeviceToHost));
		//uNorm = VectorNorm<double>(x, N);
		
		// Decide whether converge, using infinity norm
		err = fabs(uNorm-uNorm0);
		//cout<<uNorm<<" - "<<uNorm0<<" residual: "<<err<<endl;
		uNorm0 = uNorm;
		alpha = 1.0/uNorm;
		checkCublasErrors( cublasDscal(handle, N, &alpha, d_u, 1),"Dscal err!\n");
		//for (i = 0; i < N; i++)
		//	u[i] = u[i] / uNorm ;
				
		ITER++;
	}
	//delete []u0;
	//delete []u;
	cout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<err<<'\t';
	fout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<err<<'\t';
	cudaFree(d_u0);
	//delete []u0;
	//u0 = NULL;
	// return the eigenvalue
	//u_int max_index = 0;
	//u[0] = G[0] ? v[0] : 0;
	//for (i = 0; i < N; i++)
	//	if (fabs(u[i]) > fabs(u[max_index]))
	//		max_index = i;
	//cout<<uNorm<<endl;
	return (uNorm);
}




//	long long i = 0, j = 0;
//	/*double *k = new double [Ntemp];
//	for (int p=0; p<Ntemp; p++)
//	{
//		k[p] = p;
//	}*/
//	// Initialize v. Two methods are optional. Define RANDOM_V0 if you want to use random starting vector
//#ifdef RANDOM_V0
//	//srand(time(0));
//	srand(seed);
//	for (i = 0; i < Ntemp; i++){
//		v[i] = AD[i];
//		//if(i<20) cout<<v[i]<<endl;
//	}	 
//	for (i = 0; i < N; i++){
//		vv[i] = 0;
//		//if(i<20) cout<<v[i]<<endl;
//	}
//#else
//	for (i = 0; i < N && !G[i]; i++)
//#endif
//	cudaError_t cudaStat;
//	cublasStatus_t stat;
//	cublasHandle_t handle;
//	stat = cublasCreate(&handle) ;
//	checkCudaErrors( cudaMemset (d_u, 0, sizeof(double) * (Ntemp)));
//	checkCudaErrors( cudaMemcpy( d_u, v, sizeof(double) * Ntemp , cudaMemcpyHostToDevice) );
//	checkCudaErrors( cudaMemcpy( d_uu, vv, sizeof(double) * N , cudaMemcpyHostToDevice) );
//	if (stat != CUBLAS_STATUS_SUCCESS)
//		return stat;
//
//	double err1 = 1, err2 = 1;
//	int ITER = 0;
//	double vNorm = 0;
//	double temp2= -1;
//	double temp1=0;
//	//double *norm_v= new double ;
//	//double *check_du= new double [N];
//    double v_k;
//	
//	//int blocknum_spmv = N*HALF_WARP/threadnum+1;
//
//	//spmv_kernel<<<blocknum_spmv, threadnum>>>((long) N, d_r, d_c, d_G, d_u, d_u0);
//	dim3 blocknum_spmv ( Ntemp/threadnumy+(Ntemp%threadnumy?1:0) );
//	dim3 threadn(threadnumx,threadnumy);
//	//cal_k_kernel<<<6,threadnum>>>((long) Ntemp, d_AD, d_orir, d_k);
//
//	while (err1 > Epsilon &&  err2 > Epsilon && ITER < MAX_ITER)
//	{	  		
//		
//	   cublasDcopy(handle, (int) Ntemp, d_u, 1 ,d_u0, 1 );
//	   //这里需要计算vv！
//	   cal_VV_kernel<<<blocknum,threadnum>>>((long) Ntemp, (long) N, d_AD, d_uu, d_u);
//	   
//	   cublasDdot (handle, Ntemp, d_u0, 1, d_k, 1, &v_k);
//	   //calc_vector<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_k, v_k/(2*M), d_sumBG, beta, d_u0);
//	  /* checkCudaErrors( cudaMemcpy( vvector, d_vector, sizeof(double) * (Ntemp), cudaMemcpyDeviceToHost) );
//	   for (int ii = 0;ii < Ntemp;ii++)
//	   {
//		   cout<<"vector["<<ii<<"] = "<<vvector[ii]<<endl;
//	   }*/
//  		spmv_one_thread<<<blocknum , threadnum>>>((long)Ntemp, (long) M, d_u, d_r, d_c, d_uu, d_k, v_k, d_sumBG, beta, d_u0) ;
//		/*checkCudaErrors( cudaMemcpy( v, d_u, sizeof(double) * Ntemp , cudaMemcpyDeviceToHost) ); 
//	  for (int ii = 0;ii < Ntemp;ii++)
//	   {
//		   cout<<"v["<<ii<<"] = "<<v[ii]<<endl;
//	   }*/
//		//cublasDaxpy(handle, (int) Ntemp, &temp2, d_vector, 1, d_u, 1);
//		//checkCudaErrors( cudaMemcpy( v, d_u, sizeof(double) * (Ntemp), cudaMemcpyDeviceToHost) );
//	   /*for(int ii=0;ii<Ntemp;ii++)
//	   {
//			if (v[ii]!=0)
//				cout<<"v["<<ii<<"] = "<<v[ii]<<endl;
//	   }*/
//	   
//	    cublasDnrm2(handle, Ntemp, d_u, 1, &vNorm);
//		temp1=1/vNorm;
//		cublasDscal (handle, (int) Ntemp, &temp1, d_u, 1);   //Normalize v, v[i] = v[i]/vNorm
//	
//		vvplus<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_u, 1.0, d_u0, -1.0);
//	    cublasDnrm2(handle, Ntemp, d_vector, 1, &err1);
//		vvplus<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_u, 1.0, d_u0, 1.0);
//		cublasDnrm2(handle, Ntemp, d_vector, 1, &err2);
//				 
//		ITER++;
//	}	 
//	//system("pause");
//	cout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';
//	fout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';
//	
//	checkCudaErrors( cudaMemcpy( v, d_u, sizeof(double) * Ntemp , cudaMemcpyDeviceToHost) ); 
//	checkCudaErrors( cudaMemcpy(v0,d_u0, sizeof(double) * Ntemp,  cudaMemcpyDeviceToHost) );
//	cublasDestroy(handle);
//	long long max_index = 0;
//	for (i = 0; i < Ntemp; i++)
//		if (fabs(v[i]) > fabs(v[max_index]))
//			max_index = i;
//	return (v[max_index] * v0[max_index] > 0) ? vNorm: -vNorm;
//}
//
//
//

/* 
This funtion calculates dQ for each round to decide whether to split the current submodule.
It returns a double variable dQ.
*/
double calculate_dQ_GPU(double *d_S, int N, int nnz, int *d_R, int *d_C, double *d_V, double *d_K, double *d_sumBG, double m)
{
	double dQ = 0;
	//double * x = new double [N];
	int i,j;
	double temp;
	double k_s = 0;
		
	// x = S'*Bsub*S = S'*(Bsub*S) = S'*[(Asparse-k*k'/m)*S]

	checkCublasErrors( cublasDdot (handle, N, d_S, 1, d_K, 1, &k_s), "Ddot err in  calculate dQ function"); 
	/*for (i = 0; i < N; i++)
			k_s += S[i] * K[i];*/
	
	double alpha = 1.0, beta = 0;
	double *d_temp_vector;
	checkCudaErrors( cudaMalloc( (void**) &d_temp_vector, sizeof(double) * N ));
	
	checkCusparseErrors(cusparseDcsrmv(s_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nnz, &alpha, s_descr, d_V, d_R, d_C, d_S, &beta, d_temp_vector),
		"Dcsrmv err in calculate dQ function");
	
	
	/*double *temp_vector = new double [N];
	checkCudaErrors( cudaMemcpy( temp_vector, d_temp_vector, sizeof(double) * N, cudaMemcpyDeviceToHost ));
	for (i = 0; i < N; i++)
	{
		if (i<10)
			cout<<temp_vector[i]<<endl;
		dQ += temp_vector[i];
	}*/

	vvdot_subtr <<<blocknum,threadnum>>> (N, d_temp_vector, d_S, d_temp_vector, d_sumBG); //d_S.*temp_vector - d_sumBG
	
	thrust::device_ptr<double> d_tpvector(d_temp_vector);
	dQ = thrust::reduce(d_tpvector, d_tpvector+N, (double) 0, thrust::plus<double>());
	//checkCublasErrors ( cublasDasum(handle, N, d_temp_vector, 1, &dQ), "cublas sum dQ error!"); //asum is the sum of absolute value;
	
	
	/*for (i = 0; i < N; i++)
		{
			temp = 0;
			for (j = R[i]; j < R[i+1]; j++)
				temp += V[j] * S[C[j]];			
			dQ += S[i] * temp - sumBG[i];
		}*/
	
	dQ -= k_s*k_s/m;

	cudaFree(d_temp_vector);

	//delete []x;
	return (dQ);
}

double fine_tune_S_GPU(double dQ, int N, int nnz, int * d_R, int * d_C, double * d_V, double * d_K, double m, double *d_S)
{
	
	
	//signed char *Sit = new signed char [N];
	double *d_Sit;
	checkCudaErrors( cudaMalloc( (void**) &d_Sit, sizeof(double) * (N)));
	checkCudaErrors( cudaMemcpy( d_Sit, d_S, sizeof(double) * N, cudaMemcpyDeviceToDevice)) ;
	thrust::device_ptr<double> dev_Sit(d_Sit);
	thrust::device_ptr<double> dev_S(d_S);
	//memcpy(Sit,S,sizeof(signed char)*N);
	//signed char *S_Si = new signed char [N];
	
	double *d_Qit; // = new double [N];
	checkCudaErrors( cudaMalloc( (void**) &d_Qit, sizeof(double) * (N)));
	checkCudaErrors(cudaMemset(d_Qit,0,sizeof(double)*N));
	thrust::device_ptr<double> dev_Qit(d_Qit);

	char *d_indSub;
	checkCudaErrors( cudaMalloc( (void**) &d_indSub, sizeof(char) * (N)));
	thrust::device_ptr<char> dev_indSub(d_indSub);
	thrust::fill(dev_indSub,dev_indSub+N,TRUE);
	//fill(indSub,indSub+N,1); 

	
	double k_s = 0;
	int i = 0,j = 0;
	int imax = 0;
	int ITER = 0;
	double Q = dQ;
	double Qmax = dQ;
	double alpha = 0;
	//bool flag = TRUE;
	double *d_temp_vector;
	checkCudaErrors( cudaMalloc( (void**) &d_temp_vector, sizeof(double) * N ));

	//double *Qit = new double [N];
	//double *Sit = new double [N];
	//bool copy_S_flag = false;
	while (ITER<N)
	{
		//////////////////////////////////////////////////////////
		//calculate Qit[]=dQ - 4 *(Sit').*( Bsub-diag(Bsub) )*Sit
		checkCublasErrors( cublasDdot (handle, N, d_Sit, 1, d_K, 1, &k_s), "Ddot err in fine_tune_S function");
		alpha = 1.0;
		double beta = 0;
		checkCusparseErrors(cusparseDcsrmv(s_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nnz, &alpha, s_descr, d_V, d_R, d_C, d_Sit, &beta, d_temp_vector),
		"Dcsrmv err in fine_tune_S function");
		alpha = -1.0*k_s/m;
		//cublasDaxpy(cublasHandle_t handle, int n, const double *alpha, const double *x, int incx, double *y, int incy)


		checkCublasErrors(cublasDaxpy(handle, N, &alpha, d_K, 1, d_temp_vector, 1),"Daxpy err in fine_tune_S function");
		
		calculate_Qit<<<blocknum,threadnum>>>(N, Qmax, d_Qit, d_temp_vector, d_K, d_Sit, m, d_indSub);
		
		//checkCudaErrors( cudaMemcpy( Qit, d_Qit, sizeof(double) * N, cudaMemcpyDeviceToHost)) ;
		
		//for (i = 0; i < N; i++)
		//{
		//	if (indSub[i] == 0)
		//		continue;

			//for (j = 0; j < N; j++)
			//	S_Si[j] = Sit[j]*Sit[i];

			//Sit[k] = -Sit[k];
			//S_Si[i] = 0;
			//k_s = 0;
			
					

			//temp = 0;
			//for (j = R[i]; j < R[i+1]; j++)
			//	temp += S_Si[C[j]] * V[j];			
			////	Qit[k] += S[i] * temp + K[i]*K[i]/m;
			////}
			//for (j = 0; j < N; j++)
			//	temp -= S_Si[j]*K[i]*K[j]/m;
			//Qit[i] = Qmax - 4*temp;
			
			//cout<<Qit[k]<<'\t';
			////////////////////////////////////////////////////
			//Sit[k] = -Sit[k];
		//}

		//for (i = 0; i < N; i++)
		//	Qit[i] *= indSub[i];
		
		//Qmax = qmax(N, Qit, indSub, &imax);
		
		cublasIdamax(handle, N, d_Qit, 1, &imax); 
		imax = imax-1;
		Qmax = dev_Qit[imax];
		
		/*thrust::device_ptr<double> max_ptr = thrust::max_element(dev_Qit, dev_Qit + N);
		imax = &max_ptr[0] - &dev_Qit[0];
		Qmax = max_ptr[0];*/
		
		//cout<<"Qmax = "<<Qmax<<"\t imax = "<<imax<<endl;
		
		if (dev_indSub[imax])
			dev_indSub[imax] = 0;
		else
			break;
			//cout<<"indSub break flag!"<<endl;
		
		//if (Qmax > Q)
		//{
		//	Q = Qmax;
		//	copy_S_flag = true;
		//	//continue;						
		//}
		//else if (copy_S_flag)
		//{
		//	copy_S_flag = false;
		//	checkCudaErrors( cudaMemcpy(d_S, d_Sit, sizeof(double)*N, cudaMemcpyDeviceToDevice));			
		//}
				
		dev_Sit[imax] = -1.0*dev_Sit[imax];
		
		if (Qmax > Q)
		{
			Q = Qmax;
			dev_S[imax] = dev_Sit[imax];
		}
		else
			break;

		ITER++;
	}
	cout<<"\nFine tune ITER: "<<ITER<<endl;
	//if (copy_S_flag)
	//	checkCudaErrors( cudaMemcpy(d_S, d_Sit, sizeof(double)*N, cudaMemcpyDeviceToDevice));

	cudaFree(d_temp_vector);
	cudaFree(d_Qit);
	cudaFree(d_Sit);
	cudaFree(d_indSub);
	//delete []Sit;
	//delete []Qit;
	//delete []indSub;
	cout<<"fine tune Q: "<<Q<<endl;
	return Q;
}

//bool Sub_Partition_GPU(int N, int * ind, int * R, int * C, float *V, double *K, double *sumBG, double m, int * Result, int * Num_module)
//{
//	int devID;
//	cudaDeviceProp deviceProps;
//	devID = findCudaDevice();
//	// get number of SMs on this GPU
//	checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
//	
//	checkCublasErrors(cublasCreate(&handle), "Create cublas handle err!\n");
//	checkCusparseErrors( cusparseCreate(&s_handle), "Create cusparse handle err!\n" );
//	checkCusparseErrors(cusparseCreateMatDescr(&s_descr), "CreateMatDescr err!\n"); 
//	cusparseSetMatType(s_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
//    cusparseSetMatIndexBase(s_descr,CUSPARSE_INDEX_BASE_ZERO); 
//
//	int i,j;
//	int nnz = R[N]-R[0];
//	double *Vd = new double [nnz];
//	for (i = 0; i < nnz; i++)
//		Vd[i] = (double) V[i];
//	
//	int *d_R, *d_C;
//	double *d_V;
//	double *d_K, *d_sumBG;
//
//	checkCudaErrors( cudaMalloc( (void**) &d_R, sizeof(int) * (N + 1)));
//	checkCudaErrors( cudaMalloc( (void**) &d_C, sizeof(int) * (R[N]-R[0])));
//	checkCudaErrors( cudaMalloc( (void**) &d_V, sizeof(double) * (N)));
//	checkCudaErrors( cudaMalloc( (void**) &d_sumBG, sizeof(double) * (N)));
//	checkCudaErrors( cudaMalloc( (void**) &d_K, sizeof(double) * (N)));
//	
//	checkCudaErrors( cudaMemcpy( d_R, R, sizeof(int) * (N + 1), cudaMemcpyHostToDevice)) ;
//	checkCudaErrors( cudaMemcpy( d_C, C, sizeof(int) * R[N], cudaMemcpyHostToDevice)) ;
//	checkCudaErrors( cudaMemcpy( d_V, Vd, sizeof(double) * R[N], cudaMemcpyHostToDevice)) ;
//	checkCudaErrors( cudaMemcpy( d_K, K, sizeof(double) * N, cudaMemcpyHostToDevice)) ;
//	checkCudaErrors( cudaMemcpy( d_sumBG, sumBG, sizeof(double) * N, cudaMemcpyHostToDevice)) ;
//
//	double * d_eigv;
//	checkCudaErrors( cudaMalloc( (void**) &d_eigv, sizeof(double) * N));
//	
//	double lambda = 0;
//	lambda = Lead_Vector_GPU(N, nnz, d_R, d_C, d_V, d_K,  d_sumBG, m, BETA_Adjust, d_eigv);
//	//lambda = Lead_Vector(N, R, C, V, K,  sumBG, m, BETA_Adjust, eigv);
//	lambda += BETA_Adjust;
//	// If lambda < 0, calucate the leading eigenvalue for  B - lambda * I
//	if (lambda < 0)
//		lambda += Lead_Vector_GPU(N, nnz, d_R, d_C, d_V, d_K,  d_sumBG, m, lambda, d_eigv);
//		//lambda += Lead_Vector(N, R, C, V, K,  sumBG, m, lambda, eigv);
//	cout<<"Eigen Value: "<<lambda<<'\t';
//	fout<<"Eigen Value: "<<lambda<<'\t';
//
//	// Decide whether this round of partition is successful 
//	int subN = 0, subP = 0;
//	
//	/*for (i = 0; i < N; i++)
//	{
//		subP += (eigv[i] >= 0);  		
//		subN += (eigv[i] < 0);	 		
//		eigv[i] = ((eigv[i] >= 0) ? 1 : -1);
//		S[i] = ((eigv[i] >= 0) ? 1 : -1);
//	}*/
//	
//	
//	
//
//	sign_eigv<<<blocknum,threadnum>>>(N,d_eigv);
//	double S_sum = 0;
//	//stat = cublasDasum(handle, N, d_eigv, 1, &S_sum);  //asum sum of absolute value;
//	//if (stat != CUBLAS_STATUS_SUCCESS)
//	//	return stat;
//	
//	subP = (N + S_sum)/2;
//	subN = (N - S_sum)/2;
//	cout<<"subP: "<<subP<<endl;
//	cout<<"subN: "<<subN<<endl;
//	if (subP+subN != N)
//	{
//		cout<<"S allocation error!";
//		system("pause");
//	}
//
//	//calculate dQ;
//	double dQ = 0;
//
//	dQ = calculate_dQ_GPU(d_eigv, N, nnz, d_R, d_C, d_V, d_K, d_sumBG, m);
//	cout<< "dQ = "<<dQ<<endl;
//	if (dQ>DQ_MIN)  //fine tune results 
//	{
//		 dQ = fine_tune_S_GPU(dQ, N, nnz, d_R, d_C, d_V, d_K, m, d_eigv);	
//	}
//	subP = 0;
//	subN = 0;
//	
//	/*for (i = 0; i < N; i++)
//	{
//		subP += (S[i] >= 0);  		
//		subN += (S[i] < 0);	 		
//	}*/
//
//	checkCublasErrors( cublasDasum(handle, N, d_eigv, 1, &S_sum), "sum S err in sub_partition!");
//		
//	subP = (N + S_sum)/2;
//	subN = (N - S_sum)/2;
//
//	dQ = calculate_dQ_GPU(d_eigv, N, nnz, d_R, d_C, d_V, d_K, d_sumBG, m);
//	cout<< "after fine tune, dQ = "<<dQ<<endl;
//	
//	bool isSplit = (dQ > DQ_MIN && subP > MIN_GROUP && subN > MIN_GROUP);
//
//	cout<<"Divide?: "<<isSplit<<'\t';
//	fout<<"Divide?: "<<isSplit<<'\t';
//	// If not divided, return; otherwise update Result and Max_Result
//	double *S = new double [N];
//	if (isSplit)
//	{
//		checkCudaErrors( cudaMemcpy( S, d_eigv, sizeof(double) * N, cudaMemcpyDeviceToHost)) ;
//		(*Num_module) += 1;
//		if (subP>subN)
//			for (i = 0; i < N; i++)
//				Result[ind[i]] = ((S[i] >= 0 ) ? Result[ind[i]] : (*Num_module));
//		else
//			for (i = 0; i < N; i++)
//				Result[ind[i]] = ((S[i] < 0 ) ? Result[ind[i]] : (*Num_module));
//	}
//	
//	delete []S;
//	cublasDestroy(handle);
//
//
//	delete []Vd;
//	checkCudaErrors(cudaFree(d_eigv));
//	
//	checkCudaErrors(cudaFree(d_R));
//	checkCudaErrors(cudaFree(d_C));
//	checkCudaErrors(cudaFree(d_V));
//	checkCudaErrors(cudaFree(d_sumBG));
//	checkCudaErrors(cudaFree(d_K));
//	return isSplit;
//}

bool Sub_Partition_GPU_test(int N, int * ind, int * R, int * C, float *V, double *K, double *sumBG, double m, int * Result, int * Num_module)
{
		
	checkCublasErrors(cublasCreate(&handle), "Create cublas handle err!\n");
	checkCusparseErrors( cusparseCreate(&s_handle), "Create cusparse handle err in calculate dQ function" );
	checkCusparseErrors(cusparseCreateMatDescr(&s_descr), "CreateMatDescr err in calculate dQ function"); 
	cusparseSetMatType(s_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(s_descr,CUSPARSE_INDEX_BASE_ZERO); 

	int i,j;
	int nnz = R[N]-R[0];
	double *Vd = new double [nnz];
	for (i = 0; i < nnz; i++)
		Vd[i] = (double) V[i];
	
	int *d_R, *d_C;
	double *d_V;
	double *d_K, *d_sumBG;

	checkCudaErrors( cudaMalloc( (void**) &d_R, sizeof(int) * (N + 1)));
	checkCudaErrors( cudaMalloc( (void**) &d_C, sizeof(int) * nnz));
	checkCudaErrors( cudaMalloc( (void**) &d_V, sizeof(double) * nnz));
	checkCudaErrors( cudaMalloc( (void**) &d_sumBG, sizeof(double) * (N)));
	checkCudaErrors( cudaMalloc( (void**) &d_K, sizeof(double) * (N)));
	
	checkCudaErrors( cudaMemcpy( d_R, R, sizeof(int) * (N + 1), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( d_C, C, sizeof(int) * nnz, cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( d_V, Vd, sizeof(double) * nnz, cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( d_K, K, sizeof(double) * N, cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( d_sumBG, sumBG, sizeof(double) * N, cudaMemcpyHostToDevice)) ;

	double * d_eigv;
	checkCudaErrors( cudaMalloc( (void**) &d_eigv, sizeof(double) * N));
	
	double * eigv = new double [N];
	//signed char *S = new signed char [N];
	double lambda = 0;
	lambda = Lead_Vector_GPU(N, nnz, d_R, d_C, d_V, d_K,  d_sumBG, m, BETA_Adjust, d_eigv);
	//lambda = Lead_Vector(N, R, C, V, K,  sumBG, m, BETA_Adjust, eigv);
	lambda += BETA_Adjust;
	// If lambda < 0, calucate the leading eigenvalue for  B - lambda * I
	if (lambda < 0)
		lambda += Lead_Vector_GPU(N, nnz, d_R, d_C, d_V, d_K,  d_sumBG, m, lambda, d_eigv);
		//lambda += Lead_Vector(N, R, C, V, K,  sumBG, m, lambda, eigv);
		
	cout<<"Eigen Value: "<<lambda<<'\t';
	fout<<"Eigen Value: "<<lambda<<'\t';

	// Decide whether this round of partition is successful 
	int subN = 0, subP = 0;
		
	sign_eigv<<<blocknum,threadnum>>>(N,d_eigv);
	double S_sum = 0;
	thrust::device_ptr<double> dev_ptr(d_eigv);
	S_sum = thrust::reduce(dev_ptr, dev_ptr+N, (double) 0, thrust::plus<double>());
	
	subP = (N + S_sum)/2;
	subN = (N - S_sum)/2;

	cout<<"subP: "<<subP<<endl;
	cout<<"subN: "<<subN<<endl;

	
	if (subP+subN != N)
	{
		cout<<"S allocation error!";
		system("pause");
	}
	//calculate dQ;
	double dQ = 0;
	//checkCudaErrors( cudaMemcpy( d_eigv, eigv, sizeof(double) * N, cudaMemcpyHostToDevice)) ;
	dQ = calculate_dQ_GPU(d_eigv, N, nnz, d_R, d_C, d_V, d_K, d_sumBG, m);
	cout<< "dQ = "<<dQ<<endl;
	//dQ = calculate_dQ(S, N, R, C, V, K, sumBG, m);
	//cout<< "dQ = "<<dQ<<endl;
	
	/////////////////////////////////////////////////////////////////////
	//                        fine tune results                        //
	if (dQ>DQ_MIN)  
	{
		 dQ = fine_tune_S_GPU(dQ, N, nnz, d_R, d_C, d_V, d_K, m, d_eigv);	
	}
			
	//dQ = calculate_dQ_GPU(d_eigv, N, nnz, d_R, d_C, d_V, d_K, d_sumBG, m);
	cout<< "after fine tune, dQ = "<<dQ<<endl;

	S_sum = thrust::reduce(dev_ptr, dev_ptr+N, (double) 0, thrust::plus<double>());
	subP = (N + S_sum)/2;
	subN = (N - S_sum)/2;
	cout<<"subP: "<<subP<<endl;
	cout<<"subN: "<<subN<<endl;
	////////////////////////////////////////////////////////////////////////


	bool isSplit = (dQ > DQ_MIN && subP > MIN_GROUP && subN > MIN_GROUP);

	cout<<"Divide?: "<<isSplit<<'\t';
	fout<<"Divide?: "<<isSplit<<'\t';
	// If not divided, return; otherwise update Result and Max_Result
	/*double *S = new double [N];
	if (isSplit)
	{
		checkCudaErrors( cudaMemcpy( S, d_eigv, sizeof(double) * N, cudaMemcpyDeviceToHost)) ;
		(*Num_module) += 1;
		if (subP>subN)
			for (i = 0; i < N; i++)
				Result[ind[i]] = ((S[i] >= 0 ) ? Result[ind[i]] : (*Num_module));
		else
			for (i = 0; i < N; i++)
				Result[ind[i]] = ((S[i] < 0 ) ? Result[ind[i]] : (*Num_module));
	}
	
	delete []S;
	cublasDestroy(handle);*/

	if (isSplit)
	{
		checkCudaErrors( cudaMemcpy( eigv, d_eigv, sizeof(double) * N, cudaMemcpyDeviceToHost)) ;
		(*Num_module) += 1;
		if (subP>subN)
			for (i = 0; i < N; i++)
				Result[ind[i]] = ((eigv[i] >= 0 ) ? Result[ind[i]] : (*Num_module));
		else
			for (i = 0; i < N; i++)
				Result[ind[i]] = ((eigv[i] < 0 ) ? Result[ind[i]] : (*Num_module));
	}
	//delete []S;
	delete []eigv;
	delete []Vd;
	
	checkCudaErrors(cudaFree(d_eigv));
	
	checkCudaErrors(cudaFree(d_R));
	checkCudaErrors(cudaFree(d_C));
	checkCudaErrors(cudaFree(d_V));
	checkCudaErrors(cudaFree(d_sumBG));
	checkCudaErrors(cudaFree(d_K));
	cublasDestroy(handle) ;
	cusparseDestroy(s_handle) ;
	return isSplit;
}


/* 
This function does the partition, no return value.
R and C represent the adjacency matrix in CSR format.
Result stores the partition results.
*/
//void Partition_GPU(long long N, int * R, int * C, float * V, int * Result)
//{
//	long long M = R[N];
//	long long i = 0,j = 0;
//	
//	double * K;
//	K = new double [N];
//	double m = 0;
//	memset(K, 0, sizeof(double)*N);
//	for (i = 0; i < N; i++)
//	{
//		for (j = R[i]; j < R[i+1]; j++)
//		{
//			K[i]+= V[j];
//		}
//		m += K[i];
//	}
//	//double m_inv = 1.0/m;
//	cout<<"check K: "<<K[0]<<'\t'<<K[N/2]<<'\t'<<K[N-1]<<endl; //check
//	cout<<"check m: "<<m<<endl;
//	//memset(Result, 0, sizeof(int) * N);
//	fill(Result, Result+N, 1); //initialize Result
//	//int * Adjust_Result = new int [N];   //can be optimized,old version bfs style
//	int result_idx = 1;
//	int Num_module = 1;
//	int NumG = N;
//	int * index = new int [N];
//	
//	int * ind = new int [N];            //later try to put in the loop, using NumG
//	for (i = 0; i < N; i++)
//		ind[i] = i;
//	//int NumG = 0;
//	int Round = 1;						// The iteration round
//	//int Max_Result = 1;					// Maximum index of modules
//	int * R_new = new int [N+1];
//	memcpy(R_new, R, sizeof(int)*(N+1));
//	int * C_new = new int [M];
//	memcpy(C_new, C, sizeof(int)*(M));
//	float *V_new = new float [M];
//	memcpy(V_new,V,sizeof(float)*M);
//	double *K_new = new double [N];
//	memcpy(K_new, K, sizeof(double)*N);
//	double *sumBG = new double [N];
//	memset(sumBG, 0, sizeof(double)*N);
//	
//	bool isSplit;
//	int ITER = 0;
//	while (Round <= Num_module)
//	{		
//		/**********************************************************************************/
//		/************************************ Partition ***********************************/
//		if (NumG>1)
//		{
//			cout<<"\nRound:\t"<<Round<<'\t'<<"Iter:\t"<<ITER<<endl;
//			cout<<"number of nodes:\t"<<NumG<<'\t';
//			cout<<"number of non-zero elements:\t"<<R_new[NumG]<<'\t';
//			cout<<"density of this submodule:\t"<<R_new[NumG]*1.0/( (double) NumG * NumG)<<'\t';
//			fout<<"\nRound:\t"<<Round<<'\t';
//			fout<<"number of nodes:\t"<<NumG<<'\t';
//			fout.flush();
//			
//			Setup(1);
//			Start(1);
//			//isSplit = Sub_Partition(NumG, ind, R_new, C_new, V_new, K_new, sumBG, m, Result, &Num_module); //if split, Num_module+1;
//			isSplit = Sub_Partition_GPU_test(NumG, ind, R_new, C_new, V_new, K_new, sumBG, m, Result, &Num_module); //if split, Num_module+1;
//			Stop(1);
//			
//			cout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;
//			fout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;
//			//if (!isSplit)					// If divided, record the adjusted result
//			//	Adjust_Result[Round] = result_idx++;
//			//Num_module += isSplit;		// Update the total number modules , old version
//		}
//		if (!isSplit)
//			Round++;
//		/**********************************************************************************/
//		/************************** Find the next sub_module ***************************/
//		NumG = 0;
//		fill(index, index+N, -1);
//		for (i = 0; i < N; i++)
//		{
//			if (Result[i] == Round)
//			{	
//				ind[NumG] = i;
//				index[i] = NumG++;   // index[i] >= 0 if node i is involved in this round
//				//NumG ++;
//			}				
//		}
//		if (!NumG)
//		{
//			cout<<"No voxel in the submodule next round";
//			cout<<"iter: "<<ITER<<", and Num_module: "<<Num_module<<" should be equal.\n";
//			continue;
//		}
//		int ii = 0;
//		int jj = 0;
//		R_new[0] = 0;
//		for(i = 0;i < N; i++)
//		{
//			if(index[i] < 0)
//				continue;
//			K_new[ii] = K[i];                 
//			for (j = R[i];j < R[i+1];j++)
//			{
//				if(index[C[j]] < 0)
//					continue;
//				C_new[jj] = index[C[j]];
//				V_new[jj] = V[j];
//				if(C_new[jj] > NumG)                      //check flag
//					cout<<C_new[jj]<<'\t'<<C[j]<<"C_new exceed NumG!\n";				
//				jj++;
//			}
//			R_new[++ii] = jj;			
//		}
//		if (ii!=NumG)    //check flag
//			cout<<"sub module voxel# not match!";
//
//		//update R_new, C_new, V_new, and k (i.e., bsub);
//		/**********************************************************************************/
//		/******************************** diag(sum(bsub)) *********************************/
//		double temp1 = 0, temp2 = 0;
//		for (i = 0; i < NumG; i++)
//				temp2 += (K_new[i]);
//
//		for (i = 0; i < NumG; i++)
//		{
//			//if (!G[i])
//				//continue;
//			//sumBG[i] = 0;
//			temp1 = 0;
//			for (j = R_new[i]; j < R_new[i+1]; j++)
//				temp1 += V_new[j];			
//			sumBG[i] = temp1 - K_new[i] * temp2 / m;
//		}
//		
//		ITER++;
//
//		//ofstream debugfile;
//		//ostringstream s1;
//		//
//		//s1<<"Round"<<Round<<"_Iter"<<ITER;
//		////string debugfilename = "round";
//		//debugfile.open(s1.str().append("_ind"),ios::binary|ios::out);
//		//debugfile.write((char *)ind, NumG*sizeof(int));
//		//debugfile.close();
//		//int Rlength = NumG+1;
//		//int Clength = R_new[NumG];
//		//debugfile.open(s1.str().append("_csr"),ios::binary|ios::out);
//		//debugfile.write((char *)&Rlength, sizeof(int));
//		//debugfile.write((char *)R_new, Rlength*sizeof(int));
//		//debugfile.write((char *)&Clength, sizeof(int));
//		//debugfile.write((char *)C_new, R_new[NumG]*sizeof(int));
//		//debugfile.write((char *)&Clength, sizeof(int));
//		//debugfile.write((char *)V_new, R_new[NumG]*sizeof(float));
//		//debugfile.close();
//		//debugfile.open(s1.str().append("_k"),ios::binary|ios::out);
//		//debugfile.write((char *)K_new, NumG*sizeof(double));
//		//debugfile.close();
//		//debugfile.open(s1.str().append("_diag"),ios::binary|ios::out);
//		//debugfile.write((char *)sumBG, NumG*sizeof(double));
//		//debugfile.close();
//
//		
//	}
//	
//	double Q = 0;
//	for (i = 0; i < N; i++)
//		for (j = R[i]; j < R[i+1]; j++)
//			Q += V[j] * (Result[i] == Result[C[j]]);
//	for (i = 0; i < N; i++)
//		for (j = 0; j < N; j++)
//			Q -= (Result[i] == Result[j]) * (K[i]) * (K[j]) /m;
//	Q = Q / m;
//
//	cout<<"\nNumber of Modules: "<<Num_module<<",\tQ="<<Q<<endl;
//	fout<<"\nNumber of Modules: "<<Num_module<<",\tQ="<<Q<<endl;
//
//	//for (i = 0; i < N; i++)	
//	//	Result[i] = Adjust_Result[Result[i]];
//	
//	//delete []Adjust_Result;
//	delete []R_new;
//	delete []C_new;
//	delete []V_new;
//	delete []K_new;
//	delete []K;
//	delete []sumBG;
//	delete []ind;
//	delete []index;
//}

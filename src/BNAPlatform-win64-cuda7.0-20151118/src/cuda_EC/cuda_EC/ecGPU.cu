#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <memory.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include "dirent.h" 
#include "device_launch_parameters.h"  //同步函数的波浪线不管他了
#include "device_functions.h"
#include "modularity_GPU.cuh"
#include <cmath>
#include <time.h> 
#include "cublas_v2.h"
#include "cusparse.h"
#define CLEANUP(s)   printf ("%s\n", s)   //cusparse的！每步完成之后都要free东西，等所有的变量都定义完后再看吧                           
#define CUBLAS_ERROR_CHECK(sdata) if(CUBLAS_STATUS_SUCCESS!=sdata){printf("ERROR at:%s:%d\n",__FILE__,__LINE__);}//exit(-1);}  
#pragma comment(lib,"cublas.lib")
#pragma comment(lib,"cusparse.lib")

const int MAX_ITER=10000 ;			// The maximum iteration times in the power method
const int ITERNUMBER=500;
const double BETA_Adjust = 0;		// An optional parameter for quicker convergence. Its effect is uncertain
const double Epsilon = 0.000001;	// If |x - x0| < Epsilon, quit iteraion 
const double LAMBDA = 0.01;		// if labmda > LAMBDA, initiate the division
const int MIN_GROUP = 1;

const int threadnumx = 16;
const int threadnumy = 16;
const int  threadnum = 256;
const int blocknum    = 48;
extern ofstream fout;

__global__ void init_AD(long N, double *AD)
{
	int tid = blockDim.x*blockIdx.x+threadIdx.x;  //都是一维的！
	for (int i = tid; i<N; i+=blockDim.x*gridDim.x) 
		AD[i] = i;
}
__global__ void init_unweightednet(long nnz, double *V)
{
	int tid = blockDim.x*blockIdx.x+threadIdx.x;  
	for (int i = tid; i<nnz; i+=blockDim.x*gridDim.x) 
		V[i] = 1;
}
 //向量相加
 __global__ void vvplus (long N, double *result, double *v0, double alpha, double *v1, double beta)    //计算向量加法
 {
	 const int blockid   = blockIdx.x;
	 const int threadid  = threadIdx.x;
	 int offset;  	 
	 for(offset=threadid+blockid*threadnum; offset<N; offset+=threadnum*gridDim.x)
		 result[offset]=(alpha*v0[offset]+beta*v1[offset]);
 }
/**********函数功能：求对称、稀疏矩阵的特征向量中心度，参数依次分别为行偏移、列号、元素值、Ntemp*Ntemp*********************/
      
double Lead_Vector_GPU(int *R, int *C, int *V, long Ntemp,double *v)  //beta到底用不用考虑 ，beta，d_u0,d_uu关系乱了！！
{
	
	//变量定义
    cudaError_t cudaStat1,cudaStat2,cudaStat3,cudaStat4,cudaStat5,cudaStat6,cudaStat7;
	//int* d_AD_init;
	long M=R[Ntemp]/2;//	Ntemp是矩阵的行列，M是稀疏点的个数的一半！！
	int* d_r; 
    int* d_c;
	double* d_v;
	double* d_u;
	double* d_u0;
	double * d_vector;
	double* y;  
	
	long long i = 0, j = 0;

	cusparseStatus_t status;
	cusparseHandle_t cushandle;
	cusparseMatDescr_t descrA=0;
	const double alpha_mv= 1.0;
	const double beta_mv= 0.0;

	//cublas的一些参数初始化
	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;
	
	stat = cublasCreate(&handle) ;
	//cusparse的一些参数初始化，下面的英文注释是例程中固有的
	//initialize cusparse library 
	status= cusparseCreate(&cushandle); 
	if (status != CUSPARSE_STATUS_SUCCESS) 
	{  printf("CUSPARSE Library initialization failed");
	   cusparseDestroy(cushandle); 
	} 
	//create and setup matrix descriptor 
	status= cusparseCreateMatDescr(&descrA); 
	if (status != CUSPARSE_STATUS_SUCCESS) 
	{  printf("Matrix descriptor initialization failed"); 
	   cusparseDestroyMatDescr(descrA);
	   return 1;                   
	} 
	cusparseSetMatType(descrA,CUSPARSE_MATRIX_TYPE_GENERAL); 
	cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ZERO); 
	// another parameters 
	cusparseOperation_t transA= CUSPARSE_OPERATION_NON_TRANSPOSE;


	double err1 = 1, err2 = 1;
	int ITER = 0;
	double vNorm = 0;
	double temp2= -1;
	double temp1=0;

    double v_k;
	//第一步  分配空间；初始化d_u（相当于迭代公式中的x[k]）；将矩阵的R,C传到GPU上
	cudaStat1= cudaMalloc( (void**) &d_v, sizeof(double) * (2*M));
	cudaStat2= cudaMalloc( (void**) &d_r, sizeof(int) * (Ntemp + 1));
	cudaStat3= cudaMalloc( (void**) &d_c, sizeof(int) *(2*M));
	cudaStat4= cudaMalloc( (void**) &d_u, sizeof(double) * Ntemp);
	cudaStat5= cudaMalloc( (void**) &d_u0, sizeof(double) * Ntemp);
	cudaStat6= cudaMalloc( (void**) &y, sizeof(double) * Ntemp);	
	cudaStat7= cudaMalloc( (void**) &d_vector, sizeof(double) * Ntemp);
	if( (cudaStat1 != cudaSuccess)||
			(cudaStat2 != cudaSuccess)||
			(cudaStat3 != cudaSuccess)||
			(cudaStat4 != cudaSuccess)||
			(cudaStat5 != cudaSuccess)||
			(cudaStat6 != cudaSuccess)||
			(cudaStat7 != cudaSuccess))
	{
	 CLEANUP(" Device malloc failed");
	}	
    init_AD<<<blocknum,threadnum>>>(Ntemp, d_u);
	cudaStat1 = cudaMemcpy(d_r, R, 
                           (size_t)((Ntemp+1)*sizeof(d_r[0])), 
                           cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(d_c, C, 
                           (size_t)(2*M*sizeof(d_c[0])), 
						    cudaMemcpyHostToDevice);
	if ((cudaStat1 != cudaSuccess) ||
        (cudaStat2 != cudaSuccess) 
        ) {
        CLEANUP("Memcpy from Host to Device failed");
        return 1;
    }
	//第二步：判断是否为加权网络，若是，将元素值传到gpu，若不是，在gpu直接生成元素值
	if(*V==NULL)  //不加权 初始化为1
	{
	init_unweightednet<<<blocknum,threadnum>>>(2*M, d_v);
	}
	else   //加权 传进去
	{
		cudaStat1 = cudaMemcpy(d_v, V, 
                           (size_t)(2*M*sizeof(d_v[0])), 
                           cudaMemcpyHostToDevice);
	if (cudaStat1 != cudaSuccess) 
      {
        CLEANUP("Memcpy from Host to Device failed");
    }
	}
	                                         //wocao求最大值都免了！！！！        //soga!难道是前后一致！！！！！
	//第三步：循环;<1>Y(k) = X(k)/║ X(k)║∞;<2>X(k+1) = AY(k) k=0,1,2,…;<3>判断：当k充分大时，或当║ X(k)- X(k+1)║ <ε时，跳出循环;<4>结果：Y(k)≈V1,max |Xj(k)| ≈ λ1 ,1≤j≤n为x(k)的第j个分量
	while (err1 > Epsilon &&  err2 > Epsilon && ITER < MAX_ITER)
	{	  		
          //3.1先把d_u赋给d_u0
		  cublasDcopy(handle, (int) Ntemp, d_u, 1 ,d_u0, 1 );
         //3.2循环第一步 归一化	   
		  cublasDnrm2(handle, Ntemp, d_u, 1, &vNorm);            //思想：不管除以哪种范数，最后的向量收敛，计算结果一定是相同的！
		  temp1=1/vNorm;                                          //du除以du的范数，目的是为了防止精度损失，用哪种范数都行,验证过了是2范数
		  cublasDscal (handle, (int) Ntemp, &temp1, d_u, 1);   //Normalize v, v[i] = v[i]/vNorm 哦！这个是v【i】除上它自己的范数，归一化。      
	     checkCudaErrors( cudaMemcpy(v,d_u, sizeof(double) * Ntemp,  cudaMemcpyDeviceToHost) );
		  //注意上一个费时间，优化时把他优化掉
		 //3.3循环第二步 矩阵向量相乘
		  status= cusparseDcsrmv(cushandle,  transA,  Ntemp,  Ntemp,  2*M, 

		&alpha_mv,  descrA,  d_v,  d_r,      //因为这儿d_v要求必须double 所以前面都得改成double？
		d_c,  d_u,  &beta_mv,  y);
	if (status != CUSPARSE_STATUS_SUCCESS) 
	{ 
		CLEANUP("Matrix-vector multiplication failed");
		//	return 1; 
	} 
	cudaStat1 = cudaMemcpy(d_u, y, 
		(size_t)(Ntemp*sizeof(d_v[0])), 
                           cudaMemcpyDeviceToDevice);
	if (cudaStat1 != cudaSuccess) 
      {
        CLEANUP("Memcpy from Device to Device failed");
    }
       //3.4 判断
		vvplus<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_u, 1.0, d_u0, -1.0);
	    cublasDnrm2(handle, Ntemp, d_vector, 1, &err1);   //xk-x（k-1）的范数
		vvplus<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_u, 1.0, d_u0, 1.0);
		cublasDnrm2(handle, Ntemp, d_vector, 1, &err2);   //xk+x（k-1）的范数                    //这些是收敛条件，重要的参考资料是百度百科
				 
		ITER++;
	}	 
	cout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';
	fout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';
	
	//第四步:释放内存，并求最大值
	cublasDestroy(handle);
	cudaFree(y);
	cusparseDestroyMatDescr(descrA);
    cusparseDestroy(cushandle);
//	double *v = new double [Ntemp];
	double *v0 = new double [Ntemp];
//	checkCudaErrors( cudaMemcpy( v, d_u, sizeof(double) * Ntemp , cudaMemcpyDeviceToHost) ); 
	checkCudaErrors( cudaMemcpy(v0,d_u0, sizeof(double) * Ntemp,  cudaMemcpyDeviceToHost) );
/*		long long max_index = 0;
	for (i = 0; i < Ntemp; i++)
		if (fabs(v0[i]) > fabs(v0[max_index]))
			max_index = i;  */
	return vNorm ;   //what ??????????
	//return v0[max_index];
//	for (i = 0; i < Ntemp; i++)
//		v[i]/=v[max_index];
}
void main(){
	  /* create the following sparse test matrix in CSR format */
    /* |0.0     1.0 1.0|
       |    0.0 1.0    |
       |1.0 1.0 0.0 1.0|
       |1.0     1.0 0.0| */
	int C[10]={2,3,2,0,1,3,0,2};
	int R[5]={0,2,3,6,8};
	float x;
	double *v = new double [4];
	int V=NULL;
	x=Lead_Vector_GPU(R, C, &V, 4,v);
	 for(int i=0;i<4;i++)   {   
		 printf("\n%f\n",v[i]);
	 }
	 printf("eigenvalue is %f\n",x);
}

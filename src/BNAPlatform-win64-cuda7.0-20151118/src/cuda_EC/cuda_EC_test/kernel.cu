#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <memory.h>
#include <fstream>
#include <iostream>
#include <cstring>

#include "device_launch_parameters.h"  //同步函数的波浪线不管他了
#include "device_functions.h"
#include<cuda_runtime.h>
#include <cmath>
#include <time.h> 
#include "cublas_v2.h"
#include "cusparse.h"
#define CLEANUP(s)   printf ("%s\n", s)   //cusparse的！每步完成之后都要free东西，等所有的变量都定义完后再看吧                           
#define CUBLAS_ERROR_CHECK(sdata) if(CUBLAS_STATUS_SUCCESS!=sdata){printf("ERROR at:%s:%d\n",__FILE__,__LINE__);}//exit(-1);}  
#pragma comment(lib,"cublas.lib")
#pragma comment(lib,"cusparse.lib")
void main()
{
double u[4]={2.0,5.0,2.0,4.0};
double* d_u;
double vNorm = 0;
cublasHandle_t handle;
 cudaError_t cudaStat4;
 cublasStatus_t stat;
cudaStat4= cudaMalloc( (void**) &d_u, sizeof(double) * 4);
if( cudaStat4 != cudaSuccess){
CLEANUP(" Device malloc failed");
}
cudaStat4 = cudaMemcpy(d_u, u, 
                           (size_t)(4*sizeof(d_u[0])), 
                           cudaMemcpyHostToDevice);
if( cudaStat4 != cudaSuccess){
 CLEANUP("Memcpy from Host to Device failed");
}
stat = cublasCreate(&handle) ;  
cublasDnrm2(handle, 4, d_u, 1, &vNorm);
  printf("%f",vNorm);
}
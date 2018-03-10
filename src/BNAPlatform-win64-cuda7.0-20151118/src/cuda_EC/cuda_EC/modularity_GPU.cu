#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <memory.h>
#include <fstream>
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


using namespace std;

void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);

#define RANDOM_V0

extern long long N, Ntemp;
extern double * v, *vv; 
extern double * v0, * verr; 
extern	double * sumBG;
extern long long seed;

const int MAX_ITER=10000 ;			// The maximum iteration times in the power method
const int ITERNUMBER=500;
const double BETA_Adjust = 0;		// An optional parameter for quicker convergence. Its effect is uncertain
const double Epsilon = 0.000001;	// If |x - x0| < Epsilon, quit iteraion 
const double LAMBDA = 0.01;		// if labmda > LAMBDA, initiate the division
const int MIN_GROUP = 1;			// The minimum nodes of an allowed module 
extern ofstream fout;	

const int threadnumx = 16;
const int threadnumy = 16;
const int  threadnum = 256;
const int blocknum    = 48;
int* d_AD_init;
int* d_AD;
int* d_r; 
int* d_c;
int* d_orir;
double* d_u;
double* d_u0;
double* d_uu;
double * d_sumBG;
double * d_norm;
//bool* d_G;
double * temp_result;
double * d_vector;
double * d_vector1;
double * d_vector2;
double *d_k;
double *d_orik;

void Partition(int * R, int * C, int * Result);
bool Sub_Partition(int * OriR, int * OriC, int * R, int * C, int M, long long innerM, int * Result, int * Max_Result,int * AD);
double Lead_Vector(int * OriR, int * R, int * C, int M, double * sumBG1, double beta, int *AD, double *v, double *vv);

__global__ void init_AD(long N, int *AD)
{
	int tid = blockDim.x*blockIdx.x+threadIdx.x;  //都是一维的！
	for (int i = tid; i<N; i+=blockDim.x*gridDim.x) // 不着急 慢慢看
		AD[i] = i;
}

__global__ void cal_k_kernel(long N, int *AD,int *d_OriR, double *d_K )
{		
	int offset;
	const int blockid   = blockIdx.x;
	const int threadid = threadIdx.x;
	
	for(offset=(blockid/2)*threadnum*2+threadid*2+blockid%2; offset<N; offset+=blockDim.x*gridDim.x )
		d_K[offset]= (double)(d_OriR[AD[offset]+1]-d_OriR[AD[offset]]);
}

__global__ void sum_kj_kernel(long N, double *d_K ,double *sum_kj)
{	__shared__ int sum[threadnum];
	int temp=0;
	int offset;
   const int threadid =blockIdx.x*blockDim.x + threadIdx.x;
	
	//for(offset=(blockid/2)*threadnum*2+threadid*2+blockid%2; offset<N; offset+=blockDim.x*gridDim.x )
		// temp+= dG[offset]*(d_R[offset+1]-d_R[offset]);
	//for (offset=blockid*threadnum; offset+threadid<N; offset+=blockDim.x*gridDim.x ) 
	 	//  temp +=dG[offset+threadid]*(d_R[offset + threadid+1]-d_R[offset + threadid]);	 		
	
	for(offset=threadid; offset<N; offset+=blockDim.x*gridDim.x)
	{
		temp+=(int)d_K[offset] ;
	}
	sum[threadIdx.x]=temp;
	syncthreads();
	
	for(offset=1;offset+threadIdx.x<threadnum;offset*=2){
			if (threadIdx.x%(2*offset)==0)  sum[threadIdx.x]+=sum[threadIdx.x+offset] ;
			syncthreads();
	}
	if(threadIdx.x==0)
	  sum_kj[blockIdx.x]=(double) sum[threadIdx.x];
}

__global__ void sum_kernel(int size, double *data, double scale)		   //将一些前面操作中用树形结构相加每个block的结果求和。
{	
	__shared__ double sum[blocknum];
	int offset;
	sum[threadIdx.x]=data[threadIdx.x];
	syncthreads();
	
	for(offset=1;offset+threadIdx.x<blocknum; offset*=2){
		if (threadIdx.x%(2*offset)==0)  sum[threadIdx.x]+=sum[threadIdx.x+offset] ;
		syncthreads();
	}
	if(threadIdx.x==0)
	    data[0]=sum[0]*scale; 
}
               //矩阵和向量相乘  接下来搞懂两个对应函数的参数的意思，依次替换
__global__ void spmv_one_thread(long N, long M, double *result, int *R, int *C, double *vv, double *dk, double vk, double *d_sum, double beta, double *v0)	   //计算Ai*v0,   每16个threads计算一行
{
	__shared__ int R_shared[threadnum+1];

	double temp1=0;
	int offset;  
	
	for (offset=blockIdx.x*threadnum+threadIdx.x;offset<N; offset+=gridDim.x*blockDim.x)
	{
		R_shared[threadIdx.x+1] = R[offset+1];
		if(threadIdx.x==0) R_shared[threadIdx.x]=R[offset] ;
		syncthreads();	 
		temp1 =0 ;
		for(int i=R_shared[threadIdx.x]; i<R_shared[threadIdx.x+1];i++)
		{		
			temp1+=vv[C[i]];			
		}		
		temp1-=vk/(2*M)*dk[offset]+(d_sum[offset]-beta)*v0[offset];
		result[offset]=temp1 ;
			//syncthreads();
	}		
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
 
 __global__ void sumBG_kernel (double *sumBG, double *d_orik, double *d_k, int *AD,long N, double innerMd2M)    //计算向量加法
 {
	 //sumBG[i] = R[i+1] - R[i] - (OriR[AD[i]+1] - OriR[AD[i]]) * (double)innerM / 2 / M;
	 const int blockid   = blockIdx.x;
	 const int threadid  = threadIdx.x;
	 int offset;  	 
	 for(offset=threadid+blockid*threadnum; offset<N; offset+=threadnum*gridDim.x)
	 {
		sumBG[offset] = d_k[offset] - (d_orik[AD[offset]]) * innerMd2M;
	 }
 }
 
 __global__ void cal_VV_kernel ( long Ntemp, long N, int *AD, double *vv, double *v)
 {
	 const int blockid   = blockIdx.x;
	 const int threadid  = threadIdx.x;
	 int offset;
	 for(offset=threadid+blockid*threadnum; offset<Ntemp; offset+=threadnum*gridDim.x)
		 vv[AD[offset]] = v[offset];
 }

 __global__ void calc_vector (long N, double *result, double *dk, double vk2m, double *d_sum, double beta, double *v0)	   //计算 ( v_k/2/M * (R[i+1] - R[i])+(sumBG1[i] - beta) * v0[i])，
 {		
	 const int blockid   = blockIdx.x;
	 const int threadid = threadIdx.x;
	 int offset;
	 double temp=0;
	 
	 for(offset=threadnum*blockid+threadid; offset<N; offset+=threadnum*gridDim.x)
	 {	
		  temp=vk2m*dk[offset]+(d_sum[offset]-beta)*v0[offset];
		  result[offset]=temp;
	 }
  }

 /*__global__ void  Norm2_ph1(long N, double *norm,  double *v, bool *dG)	 //求向量的二范数，配合sum_kernel 得到最终结果
 {
	 __shared__  double temp[threadnum];
	 const int blockid   = blockIdx.x;
	 const int threadid = threadIdx.x;
	 double temp1=0;
	 int offset;

	 for(offset=blockid*threadnum+threadid; offset<N; offset+=threadnum*gridDim.x)
		 temp1+=dG[offset]? v[offset]*v[offset] : 0;
	 temp[threadid]=temp1;
	 syncthreads();

	 for(offset=1;offset+threadid<threadnum;offset*=2){
		 if (threadid%(2*offset)==0)  temp[threadid]+=temp[threadid+offset] ;
		 syncthreads();
	 }
	 if (threadid==0)
		 norm[blockid] = temp[threadid];  	 
 }
 __global__ void Norm2_ph2(long N,  double norm, double *v)						  //将向量归一化，若为零向量，则每个元素除以1，保持不变。
 {	 	 
	 for (int offset = threadIdx.x+blockIdx.x*threadnum; offset<N ; offset+= threadnum*gridDim.x  )
		 v[offset]/=(norm? (norm) : 1);  
 } */







/* 
This function returns the norm of the input vector x[G].
G is the logic subscriber and N is the matrix dimension.
*/


double Lead_Vector_GPU(int * OriR, int *R, int *C,  int M, double beta, int *AD)
{
	long long i = 0, j = 0;
	/*double *k = new double [Ntemp];
	for (int p=0; p<Ntemp; p++)
	{
		k[p] = p;
	}*/
	// Initialize v. Two methods are optional. Define RANDOM_V0 if you want to use random starting vector
#ifdef RANDOM_V0
	//srand(time(0));
	srand(seed);
	for (i = 0; i < Ntemp; i++){
		v[i] = AD[i];
		//if(i<20) cout<<v[i]<<endl;
	}	 
	for (i = 0; i < N; i++){
		vv[i] = 0;
		//if(i<20) cout<<v[i]<<endl;
	}
#else
	for (i = 0; i < N && !G[i]; i++)
#endif
	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;
	stat = cublasCreate(&handle) ;
	checkCudaErrors( cudaMemset (d_u, 0, sizeof(double) * (Ntemp)));
	checkCudaErrors( cudaMemcpy( d_u, v, sizeof(double) * Ntemp , cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy( d_uu, vv, sizeof(double) * N , cudaMemcpyHostToDevice) );
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;

	double err1 = 1, err2 = 1;
	int ITER = 0;
	double vNorm = 0;
	double temp2= -1;
	double temp1=0;
	//double *norm_v= new double ;
	//double *check_du= new double [N];
    double v_k;
	
	//int blocknum_spmv = N*HALF_WARP/threadnum+1;

	//spmv_kernel<<<blocknum_spmv, threadnum>>>((long) N, d_r, d_c, d_G, d_u, d_u0);
	dim3 blocknum_spmv ( Ntemp/threadnumy+(Ntemp%threadnumy?1:0) );
	dim3 threadn(threadnumx,threadnumy);
	//cal_k_kernel<<<6,threadnum>>>((long) Ntemp, d_AD, d_orir, d_k);

	while (err1 > Epsilon &&  err2 > Epsilon && ITER < MAX_ITER)
	{	  		
		
	   cublasDcopy(handle, (int) Ntemp, d_u, 1 ,d_u0, 1 );
	   //这里需要计算vv！
	   cal_VV_kernel<<<blocknum,threadnum>>>((long) Ntemp, (long) N, d_AD, d_uu, d_u);
	   
	   cublasDdot (handle, Ntemp, d_u0, 1, d_k, 1, &v_k);
	   //calc_vector<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_k, v_k/(2*M), d_sumBG, beta, d_u0);
	  /* checkCudaErrors( cudaMemcpy( vvector, d_vector, sizeof(double) * (Ntemp), cudaMemcpyDeviceToHost) );
	   for (int ii = 0;ii < Ntemp;ii++)
	   {
		   cout<<"vector["<<ii<<"] = "<<vvector[ii]<<endl;
	   }*/
  		spmv_one_thread<<<blocknum , threadnum>>>((long)Ntemp, (long) M, d_u, d_r, d_c, d_uu, d_k, v_k, d_sumBG, beta, d_u0) ;
		/*checkCudaErrors( cudaMemcpy( v, d_u, sizeof(double) * Ntemp , cudaMemcpyDeviceToHost) ); 
	  for (int ii = 0;ii < Ntemp;ii++)
	   {
		   cout<<"v["<<ii<<"] = "<<v[ii]<<endl;
	   }*/
		//cublasDaxpy(handle, (int) Ntemp, &temp2, d_vector, 1, d_u, 1);
		//checkCudaErrors( cudaMemcpy( v, d_u, sizeof(double) * (Ntemp), cudaMemcpyDeviceToHost) );
	   /*for(int ii=0;ii<Ntemp;ii++)
	   {
			if (v[ii]!=0)
				cout<<"v["<<ii<<"] = "<<v[ii]<<endl;
	   }*/
	   
	    cublasDnrm2(handle, Ntemp, d_u, 1, &vNorm);
		temp1=1/vNorm;                                          //du除以du的范数？
		cublasDscal (handle, (int) Ntemp, &temp1, d_u, 1);   //Normalize v, v[i] = v[i]/vNorm 哦！这个是v【i】除上du的范数？
	
		vvplus<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_u, 1.0, d_u0, -1.0);
	    cublasDnrm2(handle, Ntemp, d_vector, 1, &err1);
		vvplus<<<blocknum, threadnum>>>((long) Ntemp, d_vector, d_u, 1.0, d_u0, 1.0);
		cublasDnrm2(handle, Ntemp, d_vector, 1, &err2);                       //这些是收敛条件
				 
		ITER++;
	}	 
	//system("pause");
	cout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';
	fout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';
	
	checkCudaErrors( cudaMemcpy( v, d_u, sizeof(double) * Ntemp , cudaMemcpyDeviceToHost) ); 
	checkCudaErrors( cudaMemcpy(v0,d_u0, sizeof(double) * Ntemp,  cudaMemcpyDeviceToHost) );
	cublasDestroy(handle);
	long long max_index = 0;
	for (i = 0; i < Ntemp; i++)
		if (fabs(v[i]) > fabs(v[max_index]))
			max_index = i;
	return (v[max_index] * v0[max_index] > 0) ? vNorm: -vNorm; //计算最大的lamda
}



bool Sub_Partition_GPU(int * OriR, int * OriC, int * R, int * C, int M, int innerM, int * Result, int * Max_Result,int * AD)
{
	//double *v = new double [Ntemp];
	//double *vv = new double [Ntemp];
	long long i = 0, j = 0;
	long long temp1,temp2;
    double *sumBG1= new double [Ntemp];
	cublasStatus_t stat;
	cublasHandle_t handle;
	stat = cublasCreate(&handle) ;
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
	//dim3 blocknum_sumBG(  Ntemp/threadnumy+(Ntemp%threadnumy?1:0) ) ;
	//dim3 threadn(threadnumx,threadnumy);
	checkCudaErrors( cudaMemset (d_sumBG, 0, sizeof(double) * (Ntemp)));
	//sum_kj_kernel <<<blocknum, threadnum>>> ((long) N,  d_k, temp_result);
	//sum_kernel<<<1,blocknum>>>(blocknum, temp_result, 1); 
	/*for (i = 0; i < Ntemp; i++)
	{
		sumBG[i] = 0;
		sumBG[i] = R[i+1] - R[i] - (OriR[AD[i]+1] - OriR[AD[i]]) * (double)innerM / 2 / M;
	}*/
	cal_k_kernel<<<blocknum,threadnum>>>((long) Ntemp, d_AD_init, d_r, d_sumBG);
	cal_k_kernel<<<blocknum,threadnum>>>((long) Ntemp, d_AD, d_orir, d_k);
	double innerMd2M = 0.0 - (double) innerM/2/M;
	cublasDaxpy(handle, (int) Ntemp, &innerMd2M, d_k, 1, d_sumBG, 1);
	//sumBG_kernel<<< blocknum, threadnum >>>(d_sumBG, d_orik, d_k, d_AD, (long) Ntemp, (double)innerM/2/M);
	//checkCudaErrors( cudaMemcpy( sumBG1, d_sumBG, sizeof(double)*(Ntemp), cudaMemcpyDeviceToHost)) ;
	//checkCudaErrors(cudaMemcpy( d_sumBG, sumBG, sizeof(double)*(Ntemp), cudaMemcpyHostToDevice)) ;
	//checkCudaErrors( cudaMemcpy( sumBG, d_sumBG, sizeof(double)*(Ntemp), cudaMemcpyDeviceToHost)) ;
	//cout<<"innerM = "<<innerM;
	/*for(int ii=0;ii<Ntemp;ii++)
	{
		if ((sumBG1[ii]-sumBG[ii])>Epsilon || sumBG1[ii]-sumBG[ii]<-Epsilon)
			cout<<"sumBG["<<ii<<"] = "<<sumBG[ii]<<" sumBG1["<<ii<<"] = "<<sumBG1[ii]<<endl;
    }*/
	double lambda = 0;
	lambda = Lead_Vector_GPU(OriR, R, C, M, 0, AD);
	lambda -= BETA_Adjust;
	// If lambda < 0, calucate the leading eigenvalue for  B - lambda * I
	if (lambda < (-1.0)*LAMBDA)
		//lambda += Lead_Vector_GPU(  M,  G, -lambda);
		lambda += Lead_Vector_GPU(OriR, R, C, M, -lambda, AD);
		//lambda = Lead_Vector(R, C, M, sumBG, G, -lambda);
	cout<<"Eigen Value: "<<lambda<<'\t';
	fout<<"Eigen Value: "<<lambda<<'\t';

	// Decide whether this round of partition is successful 
	long long subN = 0, subP = 0;
	for (i = 0; i < Ntemp; i++)
	{
		subP += (v[i] > 0);
		subN += (v[i] <= 0);
	}
	bool Issub = (lambda > LAMBDA && subP > MIN_GROUP && subN > MIN_GROUP);

	cout<<"Divide?: "<<Issub<<"\t\n";
	fout<<"Divide?: "<<Issub<<"\t\n";
	// If not divided, return; otherwise update Result and Max_Result
	if (!Issub)
		return 0;
	for (i = 0; i < Ntemp; i++)
		Result[AD[i]] = *Max_Result + 1 + (v[i] * (subP - subN + 0.5) <= 0);
	// notice: this is wrong  Result[i] = *Max_Result + 1 + (v[i] * (subP - subN) >= 0);
	(*Max_Result) += 2;
	//delete []v;
	//delete []vv;
	return Issub;
}


/* 
This function does the partition, no return value.
R and C represent the adjacency matrix in CSR format.
Result stores the partition results.
*/
void Partition_GPU(int * R, int * C, int * Result)
{   
	int devID;
	cudaDeviceProp deviceProps;
	devID = findCudaDevice();
	// get number of SMs on this GPU
	checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
	//printf("CUDA device [%s] has %d Multi-Processors\n", deviceProps.name, deviceProps.multiProcessorCount);
   

	
	Setup(0);
	Start(0);
	int M = (R[N] - R[0]) / 2;			// The total Number connection in the network
	int Round = 0;						// The iteration round
	int Max_Result = 0;					// Maximum index of modules
	int Module_Num = 1;					// Used for adjust module index.
	bool Issub;							// Return by function Sub_Partition()
	bool * G = new bool [N];			// G[i] = 1 if node i is involved in this round of partition
	memset(Result, 0, sizeof(int) * N);
	int * Adjust_Result = new int [N];	// Map the results to consecutive intergers starting from 1 
	int Index = 1;						// Used in the adjusted results, starting from 1, increase by 1 at each successful division
	int * NewRow = new int [N+1];
	int * NewCol = new int [R[N]];
	int * Index_Result = new int [N];  // Used for matching the order number for every round partition
	long long i = 0, j = 0;
	int NumG = 0;						//
	long long innerM = 2*M;
	int Newtemp1 = 0;
	int	Newtemp2 = 0;
	while (Round <= Max_Result)	   //???????
	{
		innerM =2*M;
		NumG = 0;
		for (int i = 0; i < N; i++)
		{
			G[i] = (Result[i] == Round); 
			if (G[i])
			{
				Index_Result[NumG] = i;
				NumG ++;				// G[i] = 1 if node i is involved in this round
				continue;
			}
			innerM -= R[i+1] - R[i];
		}
		//Select all the involved node to form new row and col
		Newtemp1 = 0;
		Newtemp2 = 0;
		NewRow[Newtemp2] = Newtemp1;
		Newtemp2++;
		for(int i = 0;i < N; i++)
		{
			if(!G[i])
				continue;
			for (int j = R[i];j < R[i+1];j++)
			{
				if(!G[C[j]])
					continue;
				NewCol[Newtemp1] = C[j];
				Newtemp1++;
			}
			NewRow[Newtemp2] = Newtemp1;
			Newtemp2++;
		}
		Ntemp = Newtemp2 - 1;
		//main part of the partition
		checkCudaErrors( cudaMemcpy( d_r, NewRow, sizeof(int) * (N + 1), cudaMemcpyHostToDevice)) ;
		checkCudaErrors( cudaMemcpy( d_c, NewCol, sizeof(int) * (2*M), cudaMemcpyHostToDevice)) ;
		checkCudaErrors( cudaMemcpy( d_AD, Index_Result, sizeof(int) * (N), cudaMemcpyHostToDevice)) ;
		if (NumG)						 
		{
			cout<<"\nRound:\t"<<Round<<'\t';
			cout<<"number of nodes:\t"<<NumG<<'\t';
			fout<<"\nRound:\t"<<Round<<'\t';
			fout.flush();
			fout<<"number of nodes:\t"<<NumG<<'\t';	
			Setup(1);
			Start(1);
			if (NumG>25000)
				Issub = Sub_Partition_GPU(R,C,NewRow, NewCol, M, innerM, Result, &Max_Result, Index_Result);
			else
				Issub = Sub_Partition(R,C,NewRow, NewCol, M, innerM, Result, &Max_Result, Index_Result);
										// call Sub_Partition() for this round of division
			Stop(1);
			cout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;
			fout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;

			if (!Issub)					// If divided, record the adjusted result
				Adjust_Result[Round] = Index++;
			Module_Num += Issub;		// Update the total number modules 
		}
		Round++;
	}
	Stop(0);

	// calculate Q
	double Q = 0;
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
			Q += 1.0 * (Result[i] == Result[C[j]]);
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			Q -= 1.0 * (Result[i] == Result[j]) * (R[i+1]-R[i]) * (R[j+1]-R[j]) / 2 / M;

	Q /= 2 * M;

	cout<<"\nNumber of Modules: "<<Module_Num<<",\tQ="<<Q<<endl;
	fout<<"\nNumber of Modules: "<<Module_Num<<",\tQ="<<Q<<endl;
	cout<<"Elapsed time:   "<<GetElapsedTime(0)<<"s"<<endl;
	fout<<"Elapsed time:   "<<GetElapsedTime(0)<<"s"<<endl;

	// Adjust the results
	for (i = 0; i < N; i++)	
		Result[i] = Adjust_Result[Result[i]];
	delete []Adjust_Result;
	delete []G;
	delete []Index_Result;
	delete []NewRow;
	delete []NewCol;
	return;
}

int main(int argc, char * argv[])
{
	ofstream flog("BNA_time_log", ios::app);//ofstream是stream的子类，从内存到硬盘；ios：：app以追加方式打开文件
	clock_t total_time = clock();
	if (argc != 3) 
	{
		cerr<<"Input format: .\\Modularity.exe dir_for_csr num_of_random_networks \nFor example: .\\Modularity_CPU.exe d:\\data 10"<<endl;
		exit(1);	
	   //cerr与cout的主要区分就是,cout输出的信息可以重定向,而cerr只能输出到标准输出(显示器)上。
	}

	DIR *dp;
	struct dirent *dirp;
	if (NULL == (dp = opendir(argv[1])))
	{
		printf("can't open %s", argv[1]);
		exit (1);
	}                 //都是些文件操作而已，现在先不管他
	int FileNumber = 0;
	string filenametmp;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);

		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			FileNumber++;
		}
	}
	cout<<FileNumber<<" files to be processed."<<endl;

	closedir(dp);
	string *filename = new string[FileNumber];
	dp = opendir(argv[1]);
	long long i = 0;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			filename[i++] = filenametmp;
		}
	}

	int max_iso_n = 0;
	for (long long i = 0; i < FileNumber; i++)
	{
		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\nModular analysis for "<<a.c_str()<<" ..."<<endl;
		ifstream fin(a.c_str(), ios_base::binary);
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}
		// Read x.csr
		int Rlength = 0, Clength = 0;
		fin.read((char*)&Rlength, sizeof(int));
		int * R = new int [Rlength];
		fin.read((char*)R, sizeof(int) * Rlength);
		fin.read((char*)&Clength, sizeof(int));
		int * C = new int [Clength];
		fin.read((char*)C, sizeof(int) * Clength);
		fin.close();
		N = Rlength - 1;

		int isolated_n = 0;
		for (int j = 0; j < N; j++)
			if (R[j]==R[j+1])
				isolated_n++;
		if (isolated_n > max_iso_n)
			max_iso_n = isolated_n;

		// allocate buffers used in the iteration
		v = new double [N];
		v0 = new double [N];
		vv = new double [N];
		verr = new double [N];
		sumBG = new double [N];
		// allocate the result buffer and call Partition()
		int * Modu_Result = new int [N];	

		// Parse file name
		string X_modu = a.substr(0, a.find_last_of('.') + 1).append("modu");
		string X_cp_mas = a.substr(0, a.find_last_of('.')).append("_modu.txt");
		fout.open(X_cp_mas.c_str(), ios::out);	// Open the log file

		checkCudaErrors( cudaMalloc( (void**) &d_AD, sizeof(int) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_AD_init, sizeof(int) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_orir, sizeof(int) * (N + 1)));
		checkCudaErrors( cudaMalloc( (void**) &d_r, sizeof(int) * (N + 1)));
		checkCudaErrors( cudaMalloc( (void**) &d_c, sizeof(int) * (R[N]-R[0])));
		checkCudaErrors( cudaMalloc( (void**) &d_u, sizeof(double) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_uu, sizeof(double) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_u0, sizeof(double) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_sumBG, sizeof(double) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_k, sizeof(double) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_orik, sizeof(double) * (N)));
		checkCudaErrors( cudaMalloc( (void**) &d_vector, sizeof(double) * N));
		
	// copy host memory to device
		checkCudaErrors( cudaMemcpy( d_orir, R, sizeof(int) * (N + 1), cudaMemcpyHostToDevice)) ;
		init_AD<<<blocknum,threadnum>>>((long) N, d_AD_init);
		cal_k_kernel<<<blocknum,threadnum>>>((long) N, d_AD_init, d_orir, d_orik);
		//double *k= new double [N];
		//checkCudaErrors( cudaMemcpy( k, d_orik, sizeof(double) * (N), cudaMemcpyDeviceToHost)) ;
		//for (int x=0; x < N; x++)
		//	if ((R[x+1]-R[x]-k[x])!=0) cout<<"k["<<x<<"] = "<<R[x+1]-R[x]<<"  ;  " <<k[x] <<endl;  


		seed=time(NULL);
        Setup(0);
		Start(0);
		//Partition(R, C, Modu_Result);
		Partition_GPU(R, C, Modu_Result);
		Stop(0);
		flog<<"Modularity\t"<<a.c_str()<<"CPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;

		cout<<"Save partition results as "<<X_modu.c_str()<<endl;
		ofstream fresult;
		fresult.open(X_modu.c_str(), ios::binary|ios::out);
		fresult.write((char*)&N, sizeof(int));
		fresult.write((char*)Modu_Result, sizeof(int) * N);
		fresult.close();

		// Analysis for random networks
		int Maslov_num = atoi(argv[2]);
		cout<<"Modular analysis for random networks..."<<endl;

		int * R_dst = new int [Rlength];
		int * C_dst = new int [Clength];
		Setup(0);
		Start(0);
		for (long long l = 0; l < Maslov_num; l++)
		{
			Maslov(R_dst, C_dst, R, C, Rlength, Clength);
			checkCudaErrors( cudaMemcpy( d_r, R_dst, sizeof(int) * (N + 1), cudaMemcpyHostToDevice)) ;
			checkCudaErrors( cudaMemcpy( d_c, C_dst, sizeof(int) * ((R_dst[N]-R_dst[0])), cudaMemcpyHostToDevice)) ;
			cal_k_kernel<<<6,threadnum>>>((long) N, d_AD, d_r, d_k);
			
			//Partition(R_dst, C_dst, Modu_Result);
			Partition_GPU(R_dst, C_dst, Modu_Result);
		}
		Stop(0);
		flog<<"Modularity\tRandom"<<"CPU\t(Maslov+kernel) time\t"<<GetElapsedTime(0)<<"s"<<endl;
		// Clean up
		fout.close();
		delete []Modu_Result;
		delete []sumBG;
		delete []verr;
		delete []v0;
		delete []vv;
		delete []v;
		delete []R;
		delete []C;
		delete []R_dst;
		delete []C_dst;
		checkCudaErrors(cudaFree(d_r));
		checkCudaErrors(cudaFree(d_k));
		checkCudaErrors(cudaFree(d_orik));
		checkCudaErrors(cudaFree(d_c));
		checkCudaErrors(cudaFree (d_u));
		checkCudaErrors(cudaFree (d_uu));
		checkCudaErrors(cudaFree (d_u0));
		checkCudaErrors(cudaFree (d_sumBG));
		//checkCudaErrors(cudaFree (d_G));
		checkCudaErrors(cudaFree (temp_result));
		checkCudaErrors(cudaFree (d_vector));
		checkCudaErrors(cudaFree (d_vector1));
		checkCudaErrors(cudaFree (d_vector2));	
		checkCudaErrors(cudaFree (d_norm));
		checkCudaErrors(cudaFree (d_AD));
		checkCudaErrors(cudaFree (d_AD_init));
	}
	cout<<"==========================================================="<<endl;
	total_time = clock() - total_time;
	flog<<"Modularity\tCPU\ttotal time\t"<<1.0*total_time/1000<<"s"<<endl;
	flog<<"seed="<<seed<<endl;
	cout<<"max isolated voxel number: "<<max_iso_n;
	flog<<endl;
	flog.close();
	//system("pause");
	delete[]filename;
	return 0;
}	   

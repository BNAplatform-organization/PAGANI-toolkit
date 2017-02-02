#include "cublas_v2.h"
#include "cusparse.h"
#include "cuda_runtime.h"
#include "memory.h"
#include <iostream>
#include <ctime>
#include <fstream>
#include <vector>
#include <Windows.h>
#include<iomanip>
#pragma comment(lib,"cusparse.lib")
//#include "Utilities.cuh"

using namespace std;

#define ep  1e-6  //third question

#pragma comment(lib,"cublas.lib")
typedef float real__t;
typedef unsigned int uint__t;

#define TOM(byteValue) (byteValue/1024/1024)
extern void cusparseSafeCall(cusparseStatus_t err);

//#define CPUCormat 0

typedef struct cv
		{  
		 int column;
		 real__t value;
		} ColumnValueInfo;   



const int thread_num = 256;
const int block_num = 48;
const int blocksize = 1024*1024*48;

void select(real__t *A,long long n,long long k);
void MatrixMultiplication(real__t * BOLD_t1, real__t * BOLD_t2,real__t * out,int Batch_size,int L);

void Thrust(vector <vector<ColumnValueInfo>>::iterator begin, real__t *out, int ii, int Batch_size, real__t r_thresh, real__t er);
void ThrustAsymmetrical(vector <vector<ColumnValueInfo>>::iterator begin, real__t *out, int ii, int jj, int Batch_size, real__t r_thresh, real__t er);

__global__ void standardAndThresholdingKernel(real__t* devCormat, int Batch_size, bool diagnoal, real__t thres)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = blockDim.x * gridDim.x;
	while(i<Batch_size*Batch_size) 
	{
		if (!(devCormat[i]>(thres-ep) && devCormat[i]<=(1+ep)) ) //flag!
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

__global__ void initialone(real__t *vec,int N)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = blockDim.x * gridDim.x;
	while(i<N) 
	{
		vec[i] = 1.0;
		i += offset;
	}
}
/*************************************************************************************************/
/* GPU-based Cormat function,
output the CSR format connectivity matrix,
and functional connectivity strength.*/
/*************************************************************************************************/
int CorMat_gpu(string OutCor, real__t * BOLD_t, const int N, const int N0, const int Num_Blocks, const int L, const int Batch_size,real__t *r_thresh)
{			
	
	/*************************************************************************************************/
	/*							  Setup CUBLAS and CUSPARSE parameters                               */
	/*************************************************************************************************/
	cudaError_t cudaStat;
	cublasStatus_t stat;
	cusparseStatus_t  sparseStat;
	cublasHandle_t handle;
	cusparseHandle_t sparseHandle;
	cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;
	cusparseMatDescr_t descrA = 0;
	sparseStat= cusparseCreate(&sparseHandle);
    if (sparseStat != CUSPARSE_STATUS_SUCCESS)
		return sparseStat;
	
	/* create and setup matrix descriptor */ 
	sparseStat= cusparseCreateMatDescr(&descrA); 
	if (sparseStat != CUSPARSE_STATUS_SUCCESS) 
		return sparseStat; 
	cusparseSetMatType(descrA,CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ZERO);  
	//sparseStat = cusparse_create_mat_descr(descrA); 
	
	/***********************     Correlation matrix variables    ************************/
	/*          GPU variables        */
	real__t * devBOLD, * devCormat,*csrValA;
	int *nnzPerRowColumn,*csrRowPtrA,*csrColIndA;
	
	cudaStat = cudaMalloc ((void**)&nnzPerRowColumn, sizeof(int) * Batch_size) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
		return cudaStat;
	cudaStat = cudaMalloc ((void**)&csrRowPtrA, sizeof(int) * (Batch_size+1) ) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
		return cudaStat;
//	stat = cublasAlloc(L*N0, sizeof(real__t), (void**)&devBOLD);
	cudaStat = cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * N0) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
		return cudaStat;
//	stat = cublasAlloc(Batch_size * Batch_size, sizeof(real__t), (void**)&devCormat);		
	cudaStat = cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
		return cudaStat;
	stat = cublasSetMatrix(N0, L, sizeof(real__t), BOLD_t, N0, devBOLD, N0);
//	cudaStat = cudaMemcpy(devBOLD, BOLD_t, sizeof(real__t) * L * N0, cudaMemcpyHostToDevice);
	stat = cublasCreate(&handle) ;
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
	
	/*          CPU variables        */
	uint__t Overall_Num_Blocks = Num_Blocks * Num_Blocks;
	uint__t **Column = new uint__t* [Overall_Num_Blocks];
	real__t **Value = new real__t* [Overall_Num_Blocks];
	uint__t **Rown = new uint__t* [Overall_Num_Blocks];
	int *nnzOfEachBlock = new int [Overall_Num_Blocks];
	long long totalNonzero = 0;

	long long *Row = new long long [Num_Blocks * Batch_size+1];
	memset(Row,0,sizeof(long long)*(Num_Blocks * Batch_size+1));
		
	/****************  Functional connectivity strength variables  ****************/
	real__t *fcsGPU,*vec;
	
	cudaStat = cudaMalloc ((void**)&fcsGPU, sizeof(real__t) * Batch_size * Num_Blocks) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
		return cudaStat;	
	cudaMemset(fcsGPU,0,sizeof(real__t) * Batch_size * Num_Blocks);
	
	cudaStat = cudaMalloc ((void**)&vec, sizeof(real__t) * Batch_size) ;
	if (cudaStat != CUBLAS_STATUS_SUCCESS) 
		return cudaStat;
	initialone<<<block_num,thread_num>>>(vec,Batch_size);			                                       

	cout<<"matrix block number: "<<Num_Blocks<<endl;
	const float alpha = 1.0;
	const float beta = 0;
	//const float gamma = 1.0;
	
	clock_t correlationTime = clock();
		
	/*************************************************************************************************/
	/*						         Start the correlation computation                               */
	/*************************************************************************************************/
	for (int ii = 0; ii < Num_Blocks; ii++)
	{
		for (int jj = ii; jj < Num_Blocks; jj++)
		{
			
//#ifdef CPUCormat
//          real__t * BOLD_t1, * BOLD_t2;
//			BOLD_t1 = BOLD_t + ii * Batch_size * L;
//			BOLD_t2 = BOLD_t + jj * Batch_size * L;
//			//  real__t *v425 = new real__t[L];
//			MatrixMultiplication(BOLD_t1, BOLD_t2, out, Batch_size, L);
//#endif																																			//has been amending					
			stat = cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L,  &alpha, devBOLD + ii * Batch_size * L, L, devBOLD + jj * Batch_size * L, L, &beta, devCormat, Batch_size);//virtually kernel
			if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;
						
			//Calculating FCS
			standardAndThresholdingKernel<<<block_num,thread_num>>>(devCormat, Batch_size, ii==jj,0);  
				
			stat = cublasSgemv(handle, CUBLAS_OP_N, Batch_size, Batch_size, &alpha, devCormat, Batch_size, vec, 1, &alpha, fcsGPU + ii * Batch_size, 1);
			if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;					
			if (ii!=jj)
			{
				stat = cublasSgemv(handle, CUBLAS_OP_T, Batch_size, Batch_size, &alpha, devCormat, Batch_size, vec, 1, &alpha, fcsGPU + jj * Batch_size, 1);
				if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;
			}
			
			//so 1.thresholding 2.dense2csr 3.transfer 4.free
			//start 1.thresholding	
			standardAndThresholdingKernel<<<block_num,thread_num>>>(devCormat, Batch_size, ii==jj, *r_thresh);
			
			//2. dense2csr
			/*sparseStat = cusparseSnnz(sparseHandle, CUSPARSE_DIRECTION_ROW, Batch_size, Batch_size, descrA, devCormat,  Batch_size, nnzPerRowColumn, &nnzTotalDevHostPtr);
			if (sparseStat !=CUSPARSE_STATUS_SUCCESS)
				return sparseStat;
			cudaStat = cudaMalloc ((void**)&csrValA, sizeof(real__t) * nnzTotalDevHostPtr) ;
			if (cudaStat != CUBLAS_STATUS_SUCCESS) 
			return cudaStat;
			cudaStat = cudaMalloc ((void**)&csrColIndA, sizeof(int) * nnzTotalDevHostPtr) ;
			if (cudaStat != CUBLAS_STATUS_SUCCESS) 
				return cudaStat;
			sparseStat = cusparseSdense2csr(sparseHandle, Batch_size, Batch_size, descrA, devCormat, Batch_size, nnzPerRowColumn, csrValA, csrRowPtrA, csrColIndA);
			if (sparseStat != CUSPARSE_STATUS_SUCCESS)
				return sparseStat;*/
			int nnzTotalDevHostPtr = 0;
			sparseStat = cusparseSnnz(sparseHandle, CUSPARSE_DIRECTION_ROW, Batch_size, Batch_size, descrA, devCormat,  Batch_size, nnzPerRowColumn, &nnzTotalDevHostPtr);
			if (sparseStat !=CUSPARSE_STATUS_SUCCESS)
				return sparseStat;
			nnzOfEachBlock[ii * Num_Blocks + jj] = nnzTotalDevHostPtr;
			if (ii!=jj)
				nnzOfEachBlock[jj * Num_Blocks + ii] = nnzTotalDevHostPtr;
			
			Column[ii * Num_Blocks + jj] = new uint__t [nnzTotalDevHostPtr];
			Value[ii * Num_Blocks + jj] = new real__t [nnzTotalDevHostPtr];
			Rown[ii * Num_Blocks + jj] = new uint__t [Batch_size+1];
			if (ii!=jj)
			{
				Column[jj * Num_Blocks + ii] = new uint__t [nnzTotalDevHostPtr];
				Value[jj * Num_Blocks + ii] = new real__t [nnzTotalDevHostPtr];
				Rown[jj * Num_Blocks + ii] = new uint__t [Batch_size+1];
			}
			if (nnzTotalDevHostPtr==0)
			{
				//Rown[ii * Num_Blocks + jj] = new uint__t [Batch_size+1];
				for (int i = 0; i < (Batch_size + 1); i++)
				{
					Rown[ii * Num_Blocks + jj][i] = 0;
				}
				if (ii!=jj)
				{
					//Rown[jj * Num_Blocks + ii] = new uint__t [Batch_size+1];
					for (int i = 0; i < (Batch_size + 1); i++)
					{
						Rown[jj * Num_Blocks + ii][i] = 0;
					}
				}
				continue;
			}
			
			//malloc GPU csr column index and value
			cudaStat = cudaMalloc ((void**)&csrValA, sizeof(real__t) * nnzTotalDevHostPtr) ;
			if (cudaStat != CUBLAS_STATUS_SUCCESS)
				return cudaStat;
			cudaStat = cudaMalloc ((void**)&csrColIndA, sizeof(int) * nnzTotalDevHostPtr) ;
			if (cudaStat != CUBLAS_STATUS_SUCCESS)
				return cudaStat;
			sparseStat = cusparseSdense2csr(sparseHandle, Batch_size, Batch_size, descrA, devCormat, Batch_size, nnzPerRowColumn, csrValA, csrRowPtrA, csrColIndA);
			if (sparseStat != CUSPARSE_STATUS_SUCCESS)
				return sparseStat;

			//3. transfer
						
			cudaStat = cudaMemcpy(Column[ii * Num_Blocks + jj], csrColIndA, sizeof(int) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess)
				return cudaStat;
				 
			cudaStat = cudaMemcpy(Value[ii * Num_Blocks + jj], csrValA, sizeof(real__t) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess)
				return cudaStat;
			
			cudaStat = cudaMemcpy(Rown[ii * Num_Blocks + jj], csrRowPtrA, sizeof(int) * (Batch_size+1), cudaMemcpyDeviceToHost);
			if (cudaStat != cudaSuccess)
				return cudaStat;
			
			if(Rown[ii * Num_Blocks + jj][Batch_size]!=nnzTotalDevHostPtr)
			{
				cout<<"checking error diagnoal:"<<Rown[ii * Num_Blocks + jj][Batch_size]<<endl;
			}

			//3.5 need transposition if ii!=jj
			if (ii!=jj)
			{
				real__t *cscVal;
				int *cscRowInd, *cscColPtr;
				cudaStat = cudaMalloc ((void**)&cscVal, sizeof(real__t) * nnzTotalDevHostPtr) ;
				if (cudaStat != CUBLAS_STATUS_SUCCESS)
					return cudaStat;
				cudaStat = cudaMalloc ((void**)&cscColPtr, sizeof(int) * (Batch_size + 1)) ; //transposed R
				if (cudaStat != CUBLAS_STATUS_SUCCESS)
					return cudaStat;
				cudaStat = cudaMalloc ((void**)&cscRowInd, sizeof(int) *  nnzTotalDevHostPtr) ; //transposed C
				if (cudaStat != CUBLAS_STATUS_SUCCESS)
					return cudaStat;
				sparseStat = cusparseScsr2csc(sparseHandle, Batch_size, Batch_size, nnzTotalDevHostPtr, csrValA, csrRowPtrA, csrColIndA, cscVal, cscRowInd, cscColPtr, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
				if (sparseStat !=CUSPARSE_STATUS_SUCCESS)
					return sparseStat;
									
				cudaStat = cudaMemcpy(Column[jj * Num_Blocks + ii], cscRowInd, sizeof(int) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost);
				if (cudaStat != cudaSuccess)
					return cudaStat;
				cudaStat = cudaMemcpy(Value[jj * Num_Blocks + ii], cscVal, sizeof(real__t) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost);
				if (cudaStat != cudaSuccess)
					return cudaStat;
				cudaStat = cudaMemcpy(Rown[jj * Num_Blocks + ii], cscColPtr, sizeof(int) * (Batch_size+1), cudaMemcpyDeviceToHost);
				if (cudaStat != cudaSuccess)
				   return cudaStat;
				if (Rown[jj * Num_Blocks + ii][Batch_size]!=nnzTotalDevHostPtr)
				{
					cout<<"checking error non-diagnoal:"<<Rown[jj * Num_Blocks + ii][Batch_size]<<endl;
				}
				/*	for (int i = 0; i < nnzTotalDevHostPtr; i++)
					{
						if (Value[jj * Num_Blocks + ii][i] == 0)
						{
							cout<<"How come?"<<endl;
						}
					}*/
				cudaFree(cscVal);
				cudaFree(cscRowInd);
				cudaFree(cscColPtr);
			}

			if(ii==jj)
				totalNonzero += nnzTotalDevHostPtr;
			else
			{
				totalNonzero += nnzTotalDevHostPtr * 2 ;
			}

			//4.1 free GPU CSR column index and value.
			cudaFree(csrValA);
			cudaFree(csrColIndA);
		}
		cout<<"Fulfill the "<<ii+1<<"th disposition."<<endl;
	}

	//4.2 free nnzPerRowColumn and csrRowPtrA.
	cudaFree(nnzPerRowColumn);
	cudaFree(csrRowPtrA);
	sparseStat = cusparseDestroyMatDescr(descrA);
	if (sparseStat != CUSPARSE_STATUS_SUCCESS)
		return sparseStat;
	sparseStat = cusparseDestroy(sparseHandle);
	if (sparseStat != CUSPARSE_STATUS_SUCCESS)
		return sparseStat;


	long long Rcheck = 0;
	for (int ii = 0; ii < Num_Blocks; ii++)
	{
		for (uint__t jj = 0; jj < Num_Blocks; jj++)
		{
			Rcheck += Rown[ii*Num_Blocks+jj][Batch_size];
			for (int x = 0; x < Batch_size+1; x++)
			{
				Row[x+ii*Batch_size] += Rown[ii*Num_Blocks+jj][x];
			}
		}
		for (uint__t y = ii*Batch_size+Batch_size+1; y < (Num_Blocks * Batch_size+1); y++)
		{
			Row[y] = Row[Batch_size+ii*Batch_size];
		}
	}
	cout<<"Rcheck:"<<Rcheck<<endl;
	cout<<"Row[N]:"<<Row[N]<<endl;
	if (Row[N] != totalNonzero )
	{
		cout<<"error:"<<"R values abnormal!"<<endl;
		return false;
	}

	/****************    Write CSR networks   ****************/
	correlationTime = clock() - correlationTime;
	cout<<"correlation time: "<<correlationTime<<"ms"<<endl;
	//cout<<"overall time for histogram plus correlation: "<<*aggregrate<<"ms"<<endl;
	
	unsigned int FreeMem = 0;
	MEMORYSTATUS MemStat;
	MemStat.dwLength = sizeof(MEMORYSTATUS);
	GlobalMemoryStatus(&MemStat);
	FreeMem = TOM(MemStat.dwAvailPhys);
	cout << "bytes of physical memory: " << TOM(MemStat.dwTotalPhys) <<"M" <<endl;
	cout << "percent of memory in use: " << MemStat.dwMemoryLoad <<"%" <<endl;
	cout << "free physical memory bytes: " << TOM(MemStat.dwAvailPhys) <<"M" <<endl;
	cout<<"number of non-zero elements: "<<Row[N]<<endl;
	long long M1 = (N-1);
	M1 *= N;	
	real__t spa = 100.0 * Row[N] / M1;
	char sparsity[100];
	sprintf(sparsity, "_spa%.3f%%_cor%.3f", spa,*r_thresh);
	string Outfilename = OutCor;
	Outfilename.append(string(sparsity)).append("_weighted.csr");
	ofstream fout;
	cout<<"generating "<<Outfilename.c_str()<< "..."<<endl;
	fout.open(Outfilename.c_str(), ios::binary | ios::out);
	if (!fout)
	{
		cout<<"create outfile unsuccessfully. error code:  "<<GetLastError()<<endl;
		exit(false);
	}
	uint__t Rlength = N+1;
	fout.write((char*)&Rlength, sizeof(uint__t));
	fout.write((char*)Row, sizeof(long long) * Rlength);
	long long Clength = Row[N];
	fout.write((char*)&Clength, sizeof(long long));
	for (int ii = 0; ii < Num_Blocks; ii++)
	{
		for (int i = 0; i < Batch_size; i++)
		{
			for (int jj = 0; jj < Num_Blocks; jj++)
			{
				if (Rown[ii*Num_Blocks+jj][i]==Rown[ii*Num_Blocks+jj][i+1])
					continue;
				else
				{
					for (uint__t j = Rown[ii*Num_Blocks+jj][i]; j < Rown[ii*Num_Blocks+jj][i+1]; j++)
					{
						//cout<<Column[ii*Num_Blocks+jj][j]<<endl;
						Column[ii*Num_Blocks+jj][j] += jj * Batch_size;
						fout.write((char*)&Column[ii*Num_Blocks+jj][j], sizeof(uint__t));
					}
				}
			}
		}
	}
	fout.write((char*)&Clength, sizeof(long long));
	for (int ii = 0; ii < Num_Blocks; ii++)
	{
		for (int i = 0; i < Batch_size; i++)
		{
			for (int jj = 0; jj < Num_Blocks; jj++)
			{
				if (Rown[ii*Num_Blocks+jj][i]==Rown[ii*Num_Blocks+jj][i+1])
					continue;
				else
				{
					for (uint__t j = Rown[ii*Num_Blocks+jj][i]; j < Rown[ii*Num_Blocks+jj][i+1]; j++)
					{
						fout.write((char*)&Value[ii*Num_Blocks+jj][j], sizeof(real__t));
					}
				}
			}
		}
	}
	fout.close();
	cout<<"Transmition finished."<<endl;
	
	/****************    Write FCS information   ****************/
	real__t *fcs = new real__t[N];
	memset(fcs,0,sizeof(real__t)*N);	
	cudaStat = cudaMemcpy(fcs, fcsGPU, sizeof(real__t) * N, cudaMemcpyDeviceToHost);
	if (cudaStat != cudaSuccess) 
		return cudaStat;
	string nmOutfilename = OutCor;
	nmOutfilename.append("_weighted_fcs.nm");
	ofstream fot;
	fot.open(nmOutfilename.c_str(), ios::binary | ios::out);
	if (!fot)
	{
		cout<<"create unsuccessfully. error code:  "<<GetLastError()<<endl;
		exit(false);
	}
	int length = N;
	fot.write((char*)&length, sizeof(int));
	for (int i = 0; i < N; i++)
	{
		fot.write((char*)&fcs[i], sizeof(real__t));
	}
	fot.close();
	delete[] fcs;
	cudaFree(fcsGPU);

	return true;
}

void MatrixMultiplication(real__t * BOLD_t1, real__t * BOLD_t2,real__t * out,int Batch_size,int L)
{
	long kk = 0;
	for (int k = 0; k < Batch_size; k++)
	{
		for (int i = 0; i < Batch_size; i++)
		{   
			double sum3 = 0.0;
			for (int j = 0; j < L; j++)
			{
				sum3 += 1.0*BOLD_t1[k*L+j] * BOLD_t2[i*L+j];
			}
			out[kk++] = sum3;
		}
	}
	
}








	
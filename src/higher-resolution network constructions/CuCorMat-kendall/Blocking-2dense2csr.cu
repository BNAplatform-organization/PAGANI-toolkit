#include "cublas_v2.h"
#include "cusparse.h"
#include "cuda_runtime.h"
#include "memory.h"
#include "data_type.h"
#include <iostream>
#include <ctime>
#include <fstream>
#include <vector>
#include <Windows.h>
#include <iomanip>
#include "help_func.cuh"
#include "pre_process.cuh"
#include <thrust/device_vector.h>
#include <thrust/copy.h>

#pragma comment(lib,"cusparse.lib")
#pragma comment(lib,"cublas.lib")
using namespace std;

extern __global__ void standardAndThresholdingKernel(real__t* devCormat, int Batch_size, bool diagnoal, double thres);

extern __global__ void initialone(real__t *vec,int N);

extern void updataPointer(real__t** pointer, real__t* devBOLD, int L, int Batch_size, const int Transferblocknum);

extern void updateDevBOLD(int jj, real__t* devBOLD, real__t* BOLD_t, int L, int Batch_size, int* hosPoint, int* count, real__t* pointer, const int Num_Blocks);

/*************************************************************************************************/
/* GPU-based Cormat function,
output the CSR format connectivity matrix,
and functional connectivity strength.*/
/*************************************************************************************************/
int CorMat_gpu_blocking(string OutCor, real__t * BOLD_t, const int &N, const int &N0, const int &Num_Blocks, const int &L, const int &Batch_size, real__t* r_thresh, const int &NumS, const int Transferblocknum)
{			
	size_t L2 = L * ( L - 1 ) / 2;
	/*************************************************************************************************/
	/*							  Setup CUBLAS and CUSPARSE parameters                               */
	/*************************************************************************************************/
	//cout<<"*r_thresh:"<<*r_thresh<<endl;
	cudaError_t cudaStat;
	cublasStatus_t stat;
	cusparseStatus_t  sparseStat;
	cublasHandle_t handle;
	cusparseHandle_t sparseHandle;
	//cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;
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
	
	/*****************     Correlation matrix variables for each block   ******************/
	/*          GPU variables        */
	real__t * devBOLD, * devCormat, * devBOLD_x, * devBOLD_y, * csrValA;
	size_t * tiecount;
	real__t * pointer;
	int count = 0; //used to control pointer
	int hosPoint = Transferblocknum - 1;
	int *nnzPerRowColumn,*csrRowPtrA,*csrColIndA;
	checkCudaErrors (cudaMalloc ((void**)&nnzPerRowColumn, sizeof(int) * Batch_size)) ;
	checkCudaErrors (cudaMalloc ((void**)&csrRowPtrA, sizeof(int) * (Batch_size+1)) ) ;
	checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * Batch_size * Transferblocknum)) ;
	checkCudaErrors (cudaMalloc ( (void**)&devCormat, sizeof(real__t) * Batch_size * Batch_size)) ;
	checkCudaErrors ( cudaMalloc ((void**)&tiecount, sizeof(size_t) * N0) );
	checkCudaErrors ( cudaMalloc ((void**)&devBOLD_x, sizeof(real__t) * L2 * Batch_size) );
	checkCudaErrors ( cudaMalloc ((void**)&devBOLD_y, sizeof(real__t) * L2 * Batch_size) );
	checkCudaErrors( cudaMemcpy(devBOLD, BOLD_t, sizeof(real__t) * L * Batch_size * Transferblocknum, cudaMemcpyHostToDevice) );
	pointer = devBOLD;
	stat = cublasCreate(&handle) ;
	if (stat != CUBLAS_STATUS_SUCCESS)
		return stat;
	
	/*          CPU variables        */
	uint__t Overall_Num_Blocks = Num_Blocks * Num_Blocks;
	uint__t **Column = new uint__t* [Overall_Num_Blocks];
	real__t **Value = new real__t* [Overall_Num_Blocks];
	uint__t **Rown = new uint__t* [Overall_Num_Blocks];
	//int *nnzOfEachBlock = new int [Overall_Num_Blocks];
	R_type totalNonzero = 0;
			
	/****************  Functional connectivity strength variables  ****************/
	real__t *fcsGPU,*vec;
	
	checkCudaErrors (cudaMalloc ((void**)&fcsGPU, sizeof(real__t) * N0)) ;
	checkCudaErrors (cudaMemset(fcsGPU,0,sizeof(real__t) * N0));
	checkCudaErrors (cudaMalloc ((void**)&vec, sizeof(real__t) * Batch_size)) ;
	
	initialone<<<block_num,thread_num>>>(vec,Batch_size);			                                       

	cout<<"matrix block number: "<<Num_Blocks<<endl;
	//const float gamma = 1.0;
	
	clock_t correlationTime = clock();
		
	/*************************************************************************************************/
	/*						         Start the correlation computation                               */
	/*************************************************************************************************/
    bool flag = true;
	size_t shareSize = sizeof(real__t) * L + sizeof(size_t) * L;
	dim3 Grid(Batch_size/thread_num2D, Batch_size/thread_num2D), Block(thread_num2D,thread_num2D);
	real__t* biPointer =  devBOLD_y;
	for (int ii = 0; ii < Num_Blocks; ii++)
	{
#ifdef myDebug
		cout<<"loop flag:"<<ii<<" "<<endl;
#endif
		if ( ii > 0 )
			flag = false;
		for (int jj = ii; jj < Num_Blocks; jj++)
		{
			//1. Matrix multiplication																	//very different from step1!											
			biPointer = ii == jj ? devBOLD_y : devBOLD_x;
#ifdef figure
			int thread_num = L;
			int block_num = 180;
#endif
			pre_process <<<block_num, thread_num, shareSize>>>(biPointer, pointer, L, L2, Batch_size, tiecount + jj * Batch_size, flag);
			checkCublasErrors( cublasSgemm(handle, CUBLAS_OP_T,  CUBLAS_OP_N, Batch_size, Batch_size, L2,  &alpha, devBOLD_y, L2, biPointer, L2, &beta, devCormat, Batch_size) ); //Is this really right?
			dividedByDenominator<<<Grid, Block>>>(devCormat, Batch_size, L, tiecount + ii * Batch_size, tiecount + jj * Batch_size);
			updateDevBOLD(jj, devBOLD, BOLD_t, L, Batch_size, &hosPoint, &count, pointer, Num_Blocks); //precedence order need to be gone by
			updataPointer(&pointer, devBOLD, L, Batch_size, Transferblocknum);
			//gpuOutput(devCormat, sizeof(real__t) * Batch_size * Batch_size, OutCor, true);
#ifdef myDebug
			/*if (jj == 5)
			{
				cout<<"ii jj:"<<ii<<" "<<jj<<endl;
				thrust::copy(thrust::device_pointer_cast(devCormat + Batch_size * 837 ), thrust::device_pointer_cast(devCormat + Batch_size * 837 + 10), std::ostream_iterator<real__t>(std::cout, " "));
				cout<<"breakpoint"<<endl;	
			}*/
#endif

			//Calculating FCS
			standardAndThresholdingKernel<<<block_num, thread_num>>>(devCormat, Batch_size, ii==jj,0);  
#ifdef myDebug
			gpuOutput(devCormat, sizeof(real__t) * Batch_size * Batch_size, OutCor, false);	
			cout<<"loop flag:"<<ii<<" "<<jj<<endl;
#endif
			stat = cublasSgemv(handle, CUBLAS_OP_N, Batch_size, Batch_size, &alpha, devCormat, Batch_size, vec, 1, &alpha, fcsGPU + ii * Batch_size, 1);
			if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;					
			if (ii!=jj)
			{
				stat = cublasSgemv(handle, CUBLAS_OP_T, Batch_size, Batch_size, &alpha, devCormat, Batch_size, vec, 1, &alpha, fcsGPU + jj * Batch_size, 1);
				if (stat != CUBLAS_STATUS_SUCCESS)
				return stat;
			}
//#ifdef myDebug
//			gpuOutput(devCormat, sizeof(real__t) * Batch_size * Batch_size, OutCor, false);
//			cout<<"ii jj:"<<ii<<" "<<jj<<endl;
//#endif
			//2.thresholding	
#ifdef myLiteDebug
			//r_thresh[0] = 0.381512;
#endif
			standardAndThresholdingKernel<<<block_num,thread_num>>>(devCormat, Batch_size, ii==jj, r_thresh[0]);
			//standardAndThresholdingKernel<<<block_num,thread_num>>>(devCormat, Batch_size, ii==jj, 0);
			//3. dense2csr
			int nnzTotalDevHostPtr = 0;
			sparseStat = cusparseSnnz(sparseHandle, CUSPARSE_DIRECTION_ROW, Batch_size, Batch_size, descrA, devCormat,  Batch_size, nnzPerRowColumn, &nnzTotalDevHostPtr);
			if (sparseStat !=CUSPARSE_STATUS_SUCCESS)
				return sparseStat;
			//nnzOfEachBlock[ii * Num_Blocks + jj] = nnzTotalDevHostPtr;
			//if (ii!=jj)
			//	nnzOfEachBlock[jj * Num_Blocks + ii] = nnzTotalDevHostPtr;
			
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
			//malloc GPU csr column index and value,cusparseSdens2csr
			checkCudaErrors (cudaMalloc ((void**)&csrValA, sizeof(real__t) * nnzTotalDevHostPtr)) ;
			checkCudaErrors (cudaMalloc ((void**)&csrColIndA, sizeof(int) * nnzTotalDevHostPtr)) ;
			
			sparseStat = cusparseSdense2csr(sparseHandle, Batch_size, Batch_size, descrA, devCormat, Batch_size, nnzPerRowColumn, csrValA, csrRowPtrA, csrColIndA);
			if (sparseStat != CUSPARSE_STATUS_SUCCESS)
				return sparseStat;

			//3. transfer			
			checkCudaErrors (cudaMemcpy(Column[ii * Num_Blocks + jj], csrColIndA, sizeof(int) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost));
			checkCudaErrors (cudaMemcpy(Value[ii * Num_Blocks + jj], csrValA, sizeof(real__t) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost));
			checkCudaErrors (cudaMemcpy(Rown[ii * Num_Blocks + jj], csrRowPtrA, sizeof(int) * (Batch_size+1), cudaMemcpyDeviceToHost));
						
			if(Rown[ii * Num_Blocks + jj][Batch_size]!=nnzTotalDevHostPtr)
			{
				cout<<"checking error diagnoal:"<<Rown[ii * Num_Blocks + jj][Batch_size]<<endl;
			}

			//need transposition if ii!=jj
			if (ii!=jj)
			{
				real__t *cscVal;
				int *cscRowInd, *cscColPtr;
				
				checkCudaErrors (cudaMalloc ((void**)&cscVal, sizeof(real__t) * nnzTotalDevHostPtr)) ;
				checkCudaErrors (cudaMalloc ((void**)&cscColPtr, sizeof(int) * (Batch_size + 1))) ; //transposed R
				checkCudaErrors (cudaMalloc ((void**)&cscRowInd, sizeof(int) *  nnzTotalDevHostPtr)) ; //transposed C
				
				sparseStat = cusparseScsr2csc(sparseHandle, Batch_size, Batch_size, nnzTotalDevHostPtr, csrValA, csrRowPtrA, csrColIndA, cscVal, cscRowInd, cscColPtr, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
				if (sparseStat !=CUSPARSE_STATUS_SUCCESS)
					return sparseStat;
									
				checkCudaErrors (cudaMemcpy(Column[jj * Num_Blocks + ii], cscRowInd, sizeof(int) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost));
				checkCudaErrors (cudaMemcpy(Value[jj * Num_Blocks + ii], cscVal, sizeof(real__t) * nnzTotalDevHostPtr, cudaMemcpyDeviceToHost));
				checkCudaErrors (cudaMemcpy(Rown[jj * Num_Blocks + ii], cscColPtr, sizeof(int) * (Batch_size+1), cudaMemcpyDeviceToHost));
				
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
				checkCudaErrors (cudaFree(cscVal));
				checkCudaErrors (cudaFree(cscRowInd));
				checkCudaErrors (cudaFree(cscColPtr));
			}

			if(ii==jj)
				totalNonzero += nnzTotalDevHostPtr;
			else
			{
				totalNonzero += nnzTotalDevHostPtr * 2 ;
			}

			//4.1 free GPU CSR column index and value.
			checkCudaErrors (cudaFree(csrValA));
			checkCudaErrors (cudaFree(csrColIndA));
		}
		//cout<<"Fulfill the "<<ii+1<<"th disposition."<<endl;
	}

	//4.2 free nnzPerRowColumn and csrRowPtrA.
	checkCudaErrors (cudaFree(nnzPerRowColumn));
	checkCudaErrors (cudaFree(csrRowPtrA));

	sparseStat = cusparseDestroyMatDescr(descrA);
	if (sparseStat != CUSPARSE_STATUS_SUCCESS)
		return sparseStat;
	sparseStat = cusparseDestroy(sparseHandle);
	if (sparseStat != CUSPARSE_STATUS_SUCCESS)
		return sparseStat;

	correlationTime = clock() - correlationTime;
	cout<<"correlation time: "<<correlationTime<<"ms"<<endl;
	//cout<<"overall time for histogram plus correlation: "<<*aggregrate<<"ms"<<endl;

	/****************    Write FCS information   ****************/
	real__t *fcs = new real__t[N];
	memset(fcs,0,sizeof(real__t)*N);	
	checkCudaErrors (cudaMemcpy(fcs, fcsGPU, sizeof(real__t) * N, cudaMemcpyDeviceToHost));
	
	ofstream fcs_fout;
	string fcs_out_str = OutCor;
	fcs_out_str.append("_fcs.nm");
	fcs_fout.open(fcs_out_str.c_str(), ios::binary | ios::out);
	if (!fcs_fout)
	{
		cout<<"create unsuccessfully. error code:  "<<GetLastError()<<endl;
		exit(false);
	}
	int length = N;
	fcs_fout.write((char*)&length, sizeof(int));
	for (int i = 0; i < N; i++)
	{
		fcs_fout.write((char*)&fcs[i], sizeof(real__t));
	}
	fcs_fout.close();
	delete[] fcs;
	checkCudaErrors (cudaFree(fcsGPU));
	checkCudaErrors (cudaFree(vec));

		
	/************************** multiple thresholds ***************************/
	
	//5. vector Column Index and Value
	R_type *Row = new R_type[N+1];
	memset(Row,0,sizeof(R_type)*(N+1));
	vector<C_type> C;
	vector<V_type> V;

	//first r_threshold
	
	//R_type Rcheck = 0;
	/*for (uint__t ii = 0; ii < Num_Blocks; ii++)
	{
		for (uint__t jj = 0; jj < Num_Blocks; jj++)
		{
			Rcheck += Rown[ii*Num_Blocks+jj][Batch_size];
			for (uint__t x = 0; x < Batch_size+1; x++)
			{
				Row[x+ii*Batch_size] += Rown[ii*Num_Blocks+jj][x];
			}
		}
		for (uint__t y = ii*Batch_size+Batch_size+1; y < (Num_Blocks * Batch_size+1); y++)
		{
			Row[y] = Row[Batch_size+ii*Batch_size];
		}
	}*/

			
	
	if ( totalNonzero > C.max_size() )
	{
		cout<<"error:"<<"Vector max_size exceeds!"<<endl;
		return false;
	}	

	C.reserve(totalNonzero);
	V.reserve(totalNonzero);
	
	for (int ii = 0; ii < Num_Blocks; ii++)
	{
		for (int i = 0; i < Batch_size && ii*Batch_size+i < N; i++)
		{			
			for (int jj = 0; jj < Num_Blocks; jj++)
			{
				if (Rown[ii*Num_Blocks+jj][i]==Rown[ii*Num_Blocks+jj][i+1])
					continue;
				else
					for (uint__t j = Rown[ii*Num_Blocks+jj][i]; j < Rown[ii*Num_Blocks+jj][i+1]; j++)
					{
						C.push_back(Column[ii*Num_Blocks+jj][j] + (C_type) jj*Batch_size);
						V.push_back(Value[ii*Num_Blocks+jj][j]);						
					}
			}		
			Row[ii*Batch_size+i+1] = C.size();			
		}
		for (int jj = 0; jj < Num_Blocks; jj++)
		{
			delete[] Rown[ii*Num_Blocks+jj];
			delete[] Column[ii*Num_Blocks+jj];
			delete[] Value[ii*Num_Blocks+jj];
		}
	}

	delete[] Rown;
	delete[] Column;
	delete[] Value;
	
	
	cout<<"Row[N]:"<<Row[N]<<endl;
	if (Row[N] != totalNonzero )
	{
		cout<<"error:"<<"R values abnormal!"<<endl;
		return false;
	}
	//checking point!
	if (Row[N] != C.size() || Row[N] != V.size() )
	{
		cout<<"error:"<<"R values abnormal!"<<endl;
		return false;
	}

	MEMORYSTATUS MemStat;
	MemStat.dwLength = sizeof(MEMORYSTATUS);
	GlobalMemoryStatus(&MemStat);	
	cout << "bytes of physical memory: " << TOM(MemStat.dwTotalPhys) <<"M" <<endl;
	cout << "percent of memory in use: " << MemStat.dwMemoryLoad <<"%" <<endl;
	cout << "free physical memory bytes: " << TOM(MemStat.dwAvailPhys) <<"M" <<endl;
	cout<<"number of non-zero elements: "<<Row[N]<<endl;	
	cout<<"Transmition finished."<<endl;
	
	long long M1 = (N-1);
	M1 *= N;	
	real__t spa = 100.0 * Row[N] / M1;
	cout<<"sparsity: "<<spa<<endl;
	char sparsity[30];
	sprintf(sparsity, "_spa%.3f%%_cor%.3f", spa, r_thresh[0]);
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
#ifdef figure
	uint__t Rlength = N+1;
	fout.write((const char*)&Rlength, sizeof(uint__t));
	fout.write((const char*)Row, sizeof(R_type)*Rlength);
	R_type nnzlength = C.size();
	fout.write((const char*)&nnzlength, sizeof(R_type));
	fout.write((const char*)&C[0],sizeof(C_type)*nnzlength);
	fout.write((const char*)&nnzlength, sizeof(R_type));
	fout.write((const char*)&V[0],sizeof(V_type)*nnzlength);
#endif
	//Other thresholds
	

	for (int s = 1; s != NumS; s++)
	{
#ifdef figure
		int idx = 0;
		R_type j = Row[0];
		for (int i = 0; i != N ; i++)
		{
			for ( ; j != Row[i+1]; j++)
				if (V[j] > r_thresh[s]-ep)
				{
					C[idx] = C[j];
					V[idx] = V[j];
					++idx;
				}
			Row[i+1] = idx;
		}	

		spa = 100.0 * Row[N] / M1;
		cout<<"sparsity: "<<spa<<endl;
		char sparsity[30];
		sprintf(sparsity, "_spa%.3f%%_cor%.3f", spa, r_thresh[s]);
		Outfilename = OutCor;
		Outfilename.append(string(sparsity)).append("_weighted.csr");
		ofstream fout;
		cout<<"generating "<<Outfilename.c_str()<< "..."<<endl;
		fout.open(Outfilename.c_str(), ios::binary | ios::out);
		if (!fout)
		{
			cout<<"create outfile unsuccessfully. error code:  "<<GetLastError()<<endl;
			exit(false);
		}	
		fout.write((const char*)&Rlength, sizeof(uint__t));
		fout.write((const char*)Row, sizeof(R_type)*Rlength);
		nnzlength = Row[N];
		fout.write((const char*)&nnzlength, sizeof(R_type));
		fout.write((const char*)&C[0],sizeof(C_type)*nnzlength);
		fout.write((const char*)&nnzlength, sizeof(R_type));
		fout.write((const char*)&V[0],sizeof(V_type)*nnzlength);
#endif
	}

	C.clear();
	V.clear();
	checkCudaErrors ( cudaFree (devBOLD));
	checkCudaErrors ( cudaFree (devCormat));
	checkCudaErrors ( cudaFree (devBOLD_x));
	checkCudaErrors ( cudaFree (devBOLD_y));
	checkCudaErrors ( cudaFree (tiecount));
	fout.close();
	
	return 1;

}









	
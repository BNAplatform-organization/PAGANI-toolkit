#ifndef PRE_PROCESS_
#define PRE_PROCESS_

#include "cuda_runtime.h"

#include<string>
#include<ostream>
#include "help_func.cuh"
#include <fstream>
using namespace std;

typedef float real__t;
typedef unsigned int uint__t;

#define WARP 32
#define ep  1e-6  
#define TOM(byteValue) (byteValue/1024/1024)

const int thread_num = 1024;//for L<120,so thread_num ought to be less than 120! 
const int block_num = 30;
const int thread_num2D = 32;
const float alpha = 1.0;
const float beta = 0;

__global__ void pre_process (real__t * devBOLD, real__t * BOLD_ori, int L, int L2, int Batch_size, size_t * tiecount, bool sumtieflag );
//hahahaha! another demo with fine or coar grained comparison.
__global__ void dividedByDenominatorAndStandardedKernel(real__t* devCormat, int Batch_size, int L, real__t * tieAddr1, real__t * tieAddr2, bool diagnoal);
// 32 * 32 threads and 32 block. 
//You need try without shared memory!
__global__ void dividedByDenominatorAndStandardedKernelWith2DBlock(real__t* devCormat, int Batch_size, int L, size_t * tieAddr1, size_t * tieAddr2, bool diagnoal);
//without any set 0 operation; will also copare with nonshared memory kernel.
__global__ void dividedByDenominator(real__t* devCormat, int Batch_size, int L, size_t * tieAddr1, size_t * tieAddr2);

void gpuOutput( real__t* gpuAddr, unsigned int byteNo, string OutCor, bool nameFlag);

#endif
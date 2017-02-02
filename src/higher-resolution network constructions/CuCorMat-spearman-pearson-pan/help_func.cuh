#ifndef HELP_FUNC_CUH
#define HELP_FUNC_CUH

#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include "device_functions.h"
#include "cublas_v2.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////////
// These are CUDA Helper functions

// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
	if(cudaSuccess != err)
	{
		fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString( err ) );
		//exit(-1);
		system("pause");
	}
}

// This will output the proper error string when calling cudaGetLastError
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


// Initialization code to find the best CUDA Device
inline int findCudaDevice()
{
	cudaDeviceProp deviceProp;
	int devID = 0;
	checkCudaErrors( cudaSetDevice( devID ) );
	checkCudaErrors( cudaGetDeviceProperties(&deviceProp, devID) );
	printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);

	return devID;
}

//check cublas error
inline const char* cublasGetErrorString(cublasStatus_t status)
{
    switch(status)
    {
        case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE"; 
        case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH"; 
        case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED"; 
        case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR"; 
    }
    return "unknown error";
}

#define checkCublasErrors(err)  __checkCublasErrors (err, __FILE__, __LINE__)

inline void __checkCublasErrors(cublasStatus_t err, const char *file, const int line )
{
	if(CUBLAS_STATUS_SUCCESS != err)
	{
		fprintf(stderr, "%s(%i) : Cublas error %d: %s.\n",file, line, (int)err, cublasGetErrorString( err ) );
		//exit(-1);
		system("pause");
	}
}

#endif
// end of CUDA Helper Functions
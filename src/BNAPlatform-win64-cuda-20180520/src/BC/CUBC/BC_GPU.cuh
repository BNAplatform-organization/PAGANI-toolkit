#include <iostream>
#include <ctime>
#include <cuda_runtime.h>
#include "device_functions.h"

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
		exit(-1);        
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
	printf("GPU Device %d: \"%s\" with compute capability %d.%d\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);

	return devID;
}
// end of CUDA Helper Functions


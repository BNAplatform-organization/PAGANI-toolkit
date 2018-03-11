#include <iostream>
#include <ctime>
#include <cuda_runtime.h>
#include <cuda.h>
#include "device_functions.h"
#include "Timer.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////////
// These are CUDA Helper functions

// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
/*　函数: ctime 　　功 能: 把日期和时间转换为字符串 　　
用 法: char *ctime(const time_t *time); 　　
创建CTime对象， 使他的时间为当前时间。 　　
类函数： 　　GetMinute() 得到分钟. 　　GetSecond() 得到秒; 　　GetHour() 得到小时; 　　
GetDay() 得到 CTime持有的"天" ; 　　GetMonth() 得到月; 　　
GetDayOfWeek() 得到 CTime持有的"天"是一星期中的那一天 ; 　　GetYear() 得到年; 　　
GetTime() 返回用 __time32_t 表示的时间;
*/
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
	printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);

	return devID;
}
// end of CUDA Helper Functions


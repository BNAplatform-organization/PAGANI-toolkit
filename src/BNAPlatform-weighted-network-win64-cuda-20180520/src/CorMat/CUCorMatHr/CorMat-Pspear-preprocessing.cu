#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/functional.h>
#include <cmath>
#include <ctime>
#include "help_func.cuh"
#include "histogram.h"
using namespace std;
typedef float real__t;
const int thread_num = 1024; //maybe redefinition
const int block_num = 30; 

// square<T> computes the square of a number f(x) -> x*x 
template <typename T> 
struct square { 
	T m;
	square(T _m){ m = _m; }
	__host__ __device__ T operator()(const T& x) 
		const { return (x-m) * (x-m); } 
}; 
template<class T>
struct normalize_functor{
	T m,n;
    
	normalize_functor(T _m,T _n){
        m = _m;
		n = _n;
    }

    __host__ __device__ T operator()(T &x) const{
		return (x - m)/n;
    }
};

__global__ void assignrank(real__t* L_begin, int L, int*  addr, real__t *devTimeSeriesRank)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = blockDim.x * gridDim.x;
	//record ties
	while(i<L)
	{
		
		if (i<(L-1)&&L_begin[i]<L_begin[i+1]&&L_begin[i]==L_begin[i-1] || i ==(L-1)&&L_begin[i]==L_begin[i-1]) //be care that all elements are identical.
		{
			int rightBound = i ;
			do{i--;}
			while(L_begin[i]== L_begin[rightBound]&&i!=0);
            int leftBound = i ;
			i = threadIdx.x + blockIdx.x * blockDim.x;
			real__t averageRank = (leftBound + rightBound + 1 + 1) / 2.0; //rank is 1-based index
			for (int j = leftBound; j <= rightBound; j++)
			{
				devTimeSeriesRank[addr[j]] = averageRank; 
			}
		}
		i += offset;
	}
	i = threadIdx.x + blockIdx.x * blockDim.x;
	while(i<L)
	{
		if (devTimeSeriesRank[addr[i]] == 0 )
		{
			devTimeSeriesRank[addr[i]] = i + 1; //rank is 1-based index
		}		    
		i += offset;
	}
}

//1 block for 1 sequence. L must greater than 1024!!!!!! sacrifice generality.
__global__ void assignrankBroadcast(real__t* begin, int L, int N)
{
	int I = threadIdx.x + blockIdx.x * blockDim.x;
	for(int k = blockIdx.x; k < N; k += gridDim.x)
	{
		int i = threadIdx.x ;
		int offset = blockDim.x ;
		//real__t* addr = begin + blockIdx.x * L;
		real__t* addr = begin + k * L;
		int lessOne = 0, lessTwo = 0;
		int equalOne = 0, equalTwo = 0;
		float benchmarkOne = 0, benchmarkTwo = 0;//interval is blockDim.x

		//1.update two benchmark	
		benchmarkOne = *(addr + i);
		
		benchmarkTwo = i+offset < L ? *(addr + i + offset):0;
		
		//2.update both equal and less
# ifdef myDebug
# if __CUDA_ARCH__>=200 //requires computing capability greater than 2.0
		if (I == 0 && k == 0)
		{
			printf("benchmarkOne is %d \n", benchmarkOne);
		}
# endif
# endif
		
		for (int j = 0; j < L; j++)
		{
			real__t temp = *(addr + j);
# ifdef myDebug
# if __CUDA_ARCH__>=200
			if (I == 0 && k == 0 && j == 0)
			{
				printf("temp is %d \n", temp); //note: error usage here since temp is real__t type.
			}
# endif
# endif
			lessOne += ( temp < benchmarkOne ? 1:0 );
			equalOne += ( temp == benchmarkOne ? 1:0 );
			lessTwo += ( temp < benchmarkTwo ? 1:0 );
			equalTwo += ( temp == benchmarkTwo ? 1:0 );
		}  
		//3. compute rank and assign
		__syncthreads();
# ifdef myDebug
# if __CUDA_ARCH__>=200
		if (I == 0 && k == 0)
		{
			printf("lessOne is %d \n", lessOne);
			printf("equalOne is %d \n", equalOne);
		}
# endif
# endif
		*(addr + i) = lessOne + 0.5 * ( 1.0 + (float)equalOne );
		if (i+offset<L)
		{
			*(addr + i + offset) = lessTwo + 0.5 * ( 1.0 + (float)equalTwo );
		}
		
	}		
}
//error scope for comparison to matlab results:less than 1e-5
//coud make comparison to algorithm of single variable handle but multiple all-broadcast.
__global__ void assignrankBroadcastSharedMemory(real__t* begin, int L, int N)
{
	extern __shared__ real__t s[];
	for(int k = blockIdx.x; k < N; k += gridDim.x)
	{
		int i = threadIdx.x ;
		int offset = blockDim.x ;
		//real__t* addr = begin + blockIdx.x * L;
		real__t* addr = begin + k * L;
		int lessOne = 0, lessTwo = 0;
		int equalOne = 0, equalTwo = 0;
		real__t benchmarkOne = 0, benchmarkTwo = 0;//interval is blockDim.x

		//1.transfer data from global memory to shared memory
		s[i] = addr[i];
		if (i+offset<L)
		{
			s[i + offset] = addr[i + offset];
		}
		 __syncthreads();
		//2.update two benchmark	
		benchmarkOne = s[i];
		
		benchmarkTwo = i+offset < L ? s[i + offset]:0;
		
		//3.update both equal and less
		
		for (int j = 0; j < L; j++)
		{
			real__t temp = s[j];
			
			lessOne += temp < benchmarkOne ? 1:0;
			equalOne += temp == benchmarkOne ? 1:0;
			lessTwo += temp < benchmarkTwo ? 1:0;
			equalTwo += temp == benchmarkTwo ? 1:0;
		}  
		//4. compute rank and assign
		__syncthreads();
		*(addr + i) = lessOne + 0.5 * ( 1.0 + (float)equalOne );
		if (i+offset<L)
		{
			*(addr + i + offset) = lessTwo + 0.5 * ( 1.0 + (float)equalTwo );
		}
		
	}

}

void normalization(real__t* devBOLD, int i, int L) 
{
	thrust::plus<float> binary_op;
	float init = 0;
	real__t mean = ( thrust::reduce(thrust::device_pointer_cast(devBOLD + i * L), thrust::device_pointer_cast(devBOLD + i * L + L)) )/ L; //caution:may be integers!
	real__t norm = sqrt( thrust::transform_reduce(thrust::device_pointer_cast(devBOLD + i * L), thrust::device_pointer_cast(devBOLD + i * L + L), square<real__t>(mean), init, binary_op) );
	thrust::transform(thrust::device_pointer_cast(devBOLD + i * L), thrust::device_pointer_cast(devBOLD + i * L + L), thrust::device_pointer_cast(devBOLD + i * L), normalize_functor<real__t>(mean,norm));

}

//note: this func may exists error cause matlab check accuracy is only 1e-4 at the most.
void SpearmanAssignmentAndNormalization(real__t* BOLD_t, int N, int L)
{
	
	/**********************   allocate and transfer data   *************************/	
	real__t * devTimeSeriesRank,*devBOLD;
	checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * N)) ; //N0 to N, may cause trouble.
	checkCudaErrors (cudaMalloc ( (void**)&devTimeSeriesRank, sizeof(real__t) * L));
	checkCudaErrors (cudaMemcpy (devBOLD, BOLD_t, sizeof(real__t)* L * N, cudaMemcpyHostToDevice) );
	thrust::device_vector<int> address(L);
   
	/**********************   sort\assignRank\normalization   **********************/	
	clock_t spearTime = clock();
	for(int i = 0; i<N; i++)
	{
		//1.sorting
		thrust::sequence(address.begin(), address.end());
		thrust::sort_by_key(thrust::device_pointer_cast(devBOLD + i * L), thrust::device_pointer_cast(devBOLD + (i+1) * L), address.begin());
		//2.assign rank
		cudaMemset(devTimeSeriesRank, 0, sizeof(real__t) * L);
		assignrank<<<block_num,thread_num>>>(devBOLD + i * L, L, thrust::raw_pointer_cast(address.data()), devTimeSeriesRank);//revisit block_num later.
		checkCudaErrors(cudaMemcpy(devBOLD + i * L,devTimeSeriesRank, sizeof(real__t) * L,cudaMemcpyDeviceToDevice));
		// normalization
		normalization(devBOLD,i,L);
		
	}
	spearTime = clock() - spearTime;
	cout<<"extra-elapsed time for calculating Spearman coefficient of correlation (rank assignment plus normalization): "<<spearTime/1000.0<<"s."<<endl;
	//system("pause");

	thrust::device_vector<int>().swap(address);
	checkCudaErrors (cudaFree(devTimeSeriesRank));
	/************************  pass back data   ******************************/	
	checkCudaErrors (cudaMemcpy (BOLD_t, devBOLD, sizeof(real__t)* L * N, cudaMemcpyDeviceToHost) );
	cudaFree(devBOLD);
} 

real__t* SpearmanAssignmentAndNormalizationPointer(real__t* BOLD_t, int N, int L)
{
	
	real__t * devTimeSeriesRank,*devBOLD;
	checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * N)) ; //N0 to N, may cause trouble.
	checkCudaErrors (cudaMalloc ( (void**)&devTimeSeriesRank, sizeof(real__t) * L));
	checkCudaErrors (cudaMemcpy (devBOLD, BOLD_t, sizeof(real__t)* L * N, cudaMemcpyHostToDevice) );
	thrust::device_vector<int> address(L);
   
	/**********************   sort\assignRank\normalization   **********************/	
	clock_t spearTime = clock();
	for(int i = 0; i<N; i++)
	{
		//1.sorting
		thrust::sequence(address.begin(), address.end());
		thrust::sort_by_key(thrust::device_pointer_cast(devBOLD + i * L), thrust::device_pointer_cast(devBOLD + (i+1) * L), address.begin());
		//2.assign rank
		cudaMemset(devTimeSeriesRank, 0, sizeof(real__t) * L);
		assignrank<<<block_num,thread_num>>>(devBOLD + i * L, L, thrust::raw_pointer_cast(address.data()), devTimeSeriesRank);//revisit block_num later.
		checkCudaErrors(cudaMemcpy(devBOLD + i * L,devTimeSeriesRank, sizeof(real__t) * L,cudaMemcpyDeviceToDevice));
		// normalization
		normalization(devBOLD,i,L);
		
	}
	spearTime = clock() - spearTime;
	cout<<"extra-elapsed time for calculating Spearman coefficient of correlation (rank assignment plus normalization): "<<spearTime/1000.0<<"s."<<endl;
	
	thrust::device_vector<int>().swap(address);
	checkCudaErrors (cudaFree(devTimeSeriesRank));
	/************************  pass back pointer   ******************************/	
	return devBOLD;
} 

void SpearmanAssignmentBroadcast(real__t* BOLD_t, int N, int L, bool blocking)
{
	int blocknum = 0;
    if(blocking == true)
		blocknum = 2; //maybe future you need to define sth applied in entire project.

	if (!blocking)
	{
		real__t *devBOLD;
		checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * N)) ; //N0 to N, may cause trouble.
		checkCudaErrors (cudaMemcpy (devBOLD, BOLD_t, sizeof(real__t)* L * N, cudaMemcpyHostToDevice) );

		/**********************   sort\assignRank\normalization   **********************/	
		cout<<"assigning rank..."<<endl;
		cudaEvent_t start, stop;
		float elapsedTime;

		cudaEventCreate(&start);
		cudaEventRecord(start,0);

		//assignrankBroadcast<<<block_num,thread_num>>>(devBOLD , L, N );
#ifdef figure
		int thread_num = L;
		int block_num = 180;
#endif
		assignrankBroadcastSharedMemory<<<block_num,thread_num, L * sizeof(real__t)>>>(devBOLD, L, N);
	
		cudaEventCreate(&stop);
		cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);

		cudaEventElapsedTime(&elapsedTime, start,stop);
		printf("assignment time: %f ms\n" ,elapsedTime);

		/************************  pass back pointer   ******************************/	
		checkCudaErrors ( cudaMemcpy (BOLD_t, devBOLD, sizeof(real__t)* L * N, cudaMemcpyDeviceToHost) );
		checkCudaErrors ( cudaFree(devBOLD) );
	}else
	{
		//grnerally speaking, the smaller blocknum is, the better. 
		int generalDataWidth = N / blocknum , endDataWidth = generalDataWidth;
		real__t *devBOLD;
		if (N % blocknum != 0)
		{
			endDataWidth = N - ( blocknum - 1 ) *  generalDataWidth; //joint short tail
			checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * endDataWidth)) ; 
		}else
		{
			checkCudaErrors (cudaMalloc ((void**)&devBOLD, sizeof(real__t) * L * generalDataWidth)) ; 
		}
		
		/**********************   blocked rank assignment   **************************/
		cout<<"assigning rank (blocked transmission)..."<<endl;
		cudaEvent_t start, stop;
		float elapsedTime;

		cudaEventCreate(&start);
		cudaEventRecord(start,0);
		for (int blockid = 0; blockid < blocknum; blockid++)
		{
			int dataWidth = ( blockid == ( blocknum - 1 ) ) ?  endDataWidth : generalDataWidth;
			real__t* hostAddress = BOLD_t + blockid * dataWidth * L;
			checkCudaErrors (cudaMemcpy (devBOLD, hostAddress, sizeof(real__t)* L * dataWidth, cudaMemcpyHostToDevice) );
			
#ifdef figure
			int thread_num = L;
			int block_num = 180;
#endif
			assignrankBroadcastSharedMemory<<<block_num,thread_num, L * sizeof(real__t)>>>(devBOLD, L, dataWidth );
			//assignrankBroadcast<<<block_num,thread_num>>>(devBOLD , L, dataWidth );

			checkCudaErrors ( cudaMemcpy ( hostAddress, devBOLD, sizeof(real__t)* L * dataWidth, cudaMemcpyDeviceToHost) );

		}
		cudaEventCreate(&stop);
		cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);

		cudaEventElapsedTime(&elapsedTime, start,stop);
		printf("assignment time: %f ms\n" ,elapsedTime);

		checkCudaErrors ( cudaFree(devBOLD) );

	}
	

	
} 


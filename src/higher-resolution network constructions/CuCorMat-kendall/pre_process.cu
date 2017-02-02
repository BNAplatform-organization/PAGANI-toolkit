#include "pre_process.cuh" 

__global__ void pre_process (real__t * devBOLD, real__t * BOLD_ori, int L, int L2, int Batch_size, size_t * tiecount, bool sumtieflag )
{
	extern __shared__ real__t share[];
	real__t* BOLD_v = share;
	size_t* tcount = (size_t*)( &BOLD_v[L] );
	int tid = threadIdx.x; //each thead in distinct block
	int current_offs;
              //each block have a v   
	for(int v = blockIdx.x; v < Batch_size; v += gridDim.x)
	{
		long long offset_v_obj = v*L2; 
	
		size_t tie_count = 0;
		real__t tmp = 0;
		      
		BOLD_v[tid] = BOLD_ori[v*L+tid];//distinct block handle distinct quantity.
		if (tid + thread_num < L)
		{
			BOLD_v[tid + thread_num] = BOLD_ori[v*L + tid + thread_num];;
		}
		syncthreads();
	
		for (int ii = tid/WARP; ii <= (L-2)/2; ii+=thread_num/WARP) //question: two loop? //elements NO. 1 warp for 1 elements
		{
			current_offs = (2*L-ii-1)*ii/2;//offset address fornula: (L- 1 + L - ii) * ii / 2.0 
			for (int j = tid%WARP; j < L-1-ii; j+=WARP) //number of pair of each elements
			{
				tmp = (real__t)(BOLD_v[j+ii+1]>BOLD_v[ii]) - (BOLD_v[j+ii+1]<BOLD_v[ii]); // greater than benchmark 1; less than benchmark -1; equal to benchmark:0;
				tie_count += (tmp==0);
				devBOLD[offset_v_obj + current_offs + j] =  tmp;			
			}		
		}
	                 
		for (int ii = L-2-tid/WARP; ii >(L-2)/2; ii-=thread_num/WARP) 
		{
			current_offs = (2*L-1-ii)*ii/2;
			for (int j = tid%WARP; j < L-1-ii; j+=WARP)
			{
				tmp = (real__t) (BOLD_v[ii+j+1]>BOLD_v[ii]) - (BOLD_v[ii+j+1]<BOLD_v[ii]);
				tie_count +=  (tmp==0);
				devBOLD[offset_v_obj + current_offs + j] = tmp;			
			}
		}
		syncthreads();

		if(sumtieflag) // two definition; in order to be immune to repetitive calculation.
		{		
			tcount[tid] = tie_count;
			syncthreads();
			for (int i = thread_num/2; i > 0; i /= 2) //add together like tree
			{
				if (tid<i) tcount[tid] += tcount[tid + i];
				syncthreads();
			}
			if (tid==0)
				tiecount[v] = tcount[tid];
		}		
		syncthreads();
	}
}
//hahahaha! another demo with fine or coar grained comparison.
__global__ void dividedByDenominatorAndStandardedKernel(real__t* devCormat, int Batch_size, int L, real__t * tieAddr1, real__t * tieAddr2, bool diagnoal)
{
	int tid = threadIdx.x;
	//may cause error.
	size_t n1 = 0; // obtained in broadcasted way
	size_t n2 = 0;
	size_t n0 = L * ( L - 1 ) / 2.0; //batch_size equals to L; This is best!
	size_t temp1 = 0;
	real__t temp2 = 0;

	if ( tid >= Batch_size ) {return;}  //may cause conflict or broadcast??
	n2 = tieAddr2[tid];
	
	for(int v = blockIdx.x; v < Batch_size; v += gridDim.x)//distinct block handle distinct quantity.
	{
		real__t* nominatorAddr = devCormat + v * Batch_size;
		n1 = tieAddr1[v];
		temp1 = ( n0 - n1 ) * ( n0 - n2 );
		temp2 = nominatorAddr[tid] / temp1;
		if ( temp1 == 0 || ( diagnoal && tid == v ) || temp2<0.0f || temp2 >= (1+ep) ) //three condition need to artificially set 0: denominator is 0; diagnoal elements; thrsholding
			temp2 = 0;
		nominatorAddr[tid]  = temp2; 
	}
}
// 32 * 32 threads and 32 block. 
//You need try without shared memory!                                                                        //other kernels also need to be modified    
__global__ void dividedByDenominatorAndStandardedKernelWith2DBlock(real__t* devCormat, int Batch_size, int L, size_t * tieAddr1, size_t * tieAddr2, bool diagnoal)
{
	
	__shared__ size_t rowTile[thread_num2D];
	__shared__ size_t colTile[thread_num2D];
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	/*rowTile[threadIdx.y] = tieAddr1[row*thread_num2D+threadIdx.y];
	colTile[threadIdx.x] = tieAddr2[col*thread_num2D+threadIdx.x];*/
	
	rowTile[threadIdx.y] = tieAddr2[row];
	colTile[threadIdx.x] = tieAddr1[col];

	__syncthreads();

	size_t n0 = L * ( L - 1 ) / 2.0; 
	size_t temp1 = 0;
	real__t temp2 = 0;
	real__t* nominator = devCormat + row * Batch_size + col;

#ifdef myDebug
	if ( row == 0 && col == 4 )
	{
		printf("nominator is %f \n", *nominator);
	}
#endif
	
	temp1 = ( n0 - rowTile[threadIdx.y] ) * ( n0 - colTile[threadIdx.x] );
	temp2 = ( *nominator ) / sqrt( (real__t) temp1 );

	if ( temp1 == 0 || ( diagnoal && row == col ) || temp2<0.0f || temp2 >= (1+ep) ) //three condition need to manually set 0: denominator is 0; diagnoal elements; out of range.
		temp2 = 0;
	
	*nominator = temp2;

#ifdef myDebug
	if ( row == 0 && col == 4 )
	{
		printf("nominator after operation is %f \n", *nominator);
		printf("temp1 is %d \n", temp1);
		printf("n0 is %d \n", n0);
		printf("rowTile[threadIdx.y] is %d \n", rowTile[threadIdx.y]);
		printf("colTile[threadIdx.x] is %d \n", colTile[threadIdx.x]);

	}
#endif

}
//without any set 0 operation; will also copare with nonshared memory kernel.
__global__ void dividedByDenominator(real__t* devCormat, int Batch_size, int L, size_t * tieAddr1, size_t * tieAddr2) 
{
	__shared__ size_t rowTile[thread_num2D];
	__shared__ size_t colTile[thread_num2D];
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	/*rowTile[threadIdx.y] = tieAddr1[row*thread_num2D+threadIdx.y];
	colTile[threadIdx.x] = tieAddr2[col*thread_num2D+threadIdx.x];
	*/
	
	rowTile[threadIdx.y] = tieAddr2[row];
	colTile[threadIdx.x] = tieAddr1[col];
	
	__syncthreads();
	
	size_t n0 = L * ( L - 1 ) / 2.0; 
	size_t temp1 = 0;
	real__t* nominator = devCormat + row * Batch_size + col;

	temp1 = ( n0 - rowTile[threadIdx.y] ) * ( n0 - colTile[threadIdx.x] );

	(*nominator) /= sqrt ( (real__t) temp1);
}

void gpuOutput( real__t* gpuAddr, unsigned int byteNo, string OutCor, bool nameFlag)
{
	ofstream fout;
	real__t* cpuAddr = new real__t[byteNo/sizeof(real__t)];
	checkCudaErrors ( cudaMemcpy(cpuAddr, gpuAddr, byteNo, cudaMemcpyDeviceToHost) );
	string filename;
	if(!nameFlag)
		filename = OutCor.append("_wrong.matrix");
	else
		filename = OutCor.append("_right.matrix");
	fout.open(filename.c_str(), ios::binary | ios::out);
	if (!fout)
	{
		cout<<"create outfile(gpu) unsuccessfully. error code:  "<<"fighting!"<<endl;		
		system("pause");
	}	
	fout.write((const char*)cpuAddr,byteNo);
	fout.close();
	delete[] cpuAddr;
}
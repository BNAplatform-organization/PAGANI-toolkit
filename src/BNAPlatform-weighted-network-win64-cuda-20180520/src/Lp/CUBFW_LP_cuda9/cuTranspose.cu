__global__ void Transpose_ker(float * dst, float * src, int size)
{
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size && j < size)
		dst[i * size + j] = src[j * size + i];
}

void cuTranspose(float * dst, float * src, int size)
{
	size = (size + 16 - 1) / 16 * 16;
	dim3 dimBlock(16, 16); 
	dim3 dimGrid(size / 16, size / 16);
	Transpose_ker<<<dimGrid, dimBlock>>>(dst, src, size);
}
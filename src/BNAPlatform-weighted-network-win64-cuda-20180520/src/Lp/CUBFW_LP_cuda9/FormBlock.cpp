# include <memory.h>
# include <iostream>
using namespace std;

bool FormBlock(float * B, float * A, int N, int Block_size)
{
	int Num_block = N / Block_size;
	int i, j, l;
	float * src, * dst;
	dst = B;
	for (i = 0; i < Num_block; i++)
		for (j = 0; j < Num_block; j++)
		{
			src = A + ((long long)i*N+j) * Block_size;
			for (l = 0; l < Block_size; l++)
			{
				memcpy(dst, src, sizeof(float) * Block_size);
				dst += Block_size;
				src += N;
			}
		}
	return true;
}

bool DeFormBlock(float * A, float * B, int N, int Block_size)
{
	int Num_block = N / Block_size;
	int i, j, l;
	float * src, * dst;
	src = B;
	for (i = 0; i < Num_block; i++)
		for (j = 0; j < Num_block; j++)
		{
			dst = A + ((long long)i*N+j) * Block_size;
			for (l = 0; l < Block_size; l++)
			{
				memcpy(dst, src, sizeof(float) * Block_size);
				src += Block_size;
				dst += N;
			}
		}
	return true;
}

void init_block(float *costmat, const int *row, const int *col, const float *power, const int numVertices, const int block_size)
{
	int k = 0;
	int	block_cnt = numVertices / block_size;
	for (long long i = 0; i < (long long)numVertices * numVertices; i++)
		costmat[i] = numeric_limits<float>::infinity();

	for ( int i = 0; i < numVertices; i ++ ) 
	{
		for ( int j = row[i]; j < row[i+1]; j ++ ) 
		{
			int ab_row 		= i;
			int ab_col 		= col[j];
			int block_row 	= ab_row / block_size;
			int block_col	= ab_col / block_size;
			int block_i		= ab_row % block_size;
			int block_j		= ab_col % block_size;
			long long offset 		= block_row * block_cnt + block_col;	
			offset *= block_size*block_size;
			costmat[block_i * block_size + block_j + offset] = 1/power[j]; //definition of distance
		}
	}

	for (int k = 0; k < block_cnt; k++)
	{
		for (int i = 0; i < block_size; i++)
		{
			costmat[((long long)k * block_cnt + k) * block_size * block_size + i * block_size + i] = 0;
		}
	}
}

void post_block (float *output, float *costmat, int dim, int block_size ) 
{
	int block_cnt	= (dim + block_size - 1)/ block_size;
	long long  index = 0;
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
		{
			int block_row = i / block_size;
			int block_col = j / block_size;
			int block_i   = i % block_size;
			int block_j   = j % block_size;
			long long  offset    = block_row * block_cnt + block_col;
			offset *= block_size*block_size;
			output[index] = costmat[block_i*block_size + block_j + offset];
			index ++;
		}
}

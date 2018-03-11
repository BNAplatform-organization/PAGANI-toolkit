void init_block(float *costmat, const int *row, const int *col,const float *power, const int numVertices, const int sizeBlock);
void post_block (float *output, float *costmat, int dim, int block_size);

bool FormBlock(float * B, float * A, int N, int Block_size);
bool DeFormBlock(float * A, float * B, int N, int Block_size);
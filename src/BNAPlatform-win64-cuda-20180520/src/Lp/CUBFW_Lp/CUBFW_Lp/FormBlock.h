#ifndef __FORM_BLOCK__
#define __FORM_BLOCK__

#include "data_type.h"

void init_block(float *costmat, const R_type *row, const C_type *col, const int numVertices, const int sizeBlock);
void post_block (float *output, float *costmat, int dim, int block_size);

bool FormBlock(float * B, float * A, int N, int Block_size);
bool DeFormBlock(float * A, float * B, int N, int Block_size);

#endif
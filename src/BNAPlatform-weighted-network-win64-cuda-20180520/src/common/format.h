#ifndef FORMAT_H
#define FORMAT_H

typedef struct {
	int c;
	int r;
} Sparse;

bool sparse_less ( const Sparse &, const Sparse & );
void CSR_sort ( float*, float*, int, int, int*, int*, int* );


#endif
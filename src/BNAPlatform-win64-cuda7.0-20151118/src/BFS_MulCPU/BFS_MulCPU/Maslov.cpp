
#include <stdlib.h>
#include <memory.h>
#include <ctime>
//#include <iostream>
using namespace std;

void Rewire(bool ** G, int * R, int * C, int N, int M)
/*
input:
	G :		The adjacency matrix
	R & C :	Rows and columns of each edge
	N :		Total number of nodes
	M :		Total number of edges	
*/

{
	int V1, V2, V3, V4;			// Four nodes
	int Try = 2 * M;			// Total shuffling times
	srand(time(0));
	int i = 0, j = 0, k = 0;
	int Shuffle_times = 0;
	for (k = 0; k < Try; k++)
	{
		i = (rand() * (0x7fff+1) + rand()) % M;	// rand() returns just from 0 to 0x7fff on windows. Too narrow for M
		V1 = R[i];	
		V2 = C[i];
		j = (rand() * (0x7fff+1) + rand()) % M;
		V3 = R[j];	
		V4 = C[j];				// Select two edges 1-2 3-4
		if (V1 == V3 || V1 == V4 || V2 == V3 || V2 == V4)
			continue;		
		// If two edges have a node in common, this shuffling effort is aborted
		if (rand() % 2)
		{	
		// 50% Shuffles to 1-3 2-4	
			if (G[V1][V3] || G[V2][V4])
				continue;		// If 1-3 or 2-4 already connects, this shuffling effort is aborted
			G[V1][V2] = G[V2][V1] = 0;
			G[V3][V4] = G[V4][V3] = 0;
			G[V1][V3] = G[V3][V1] = 1;
			G[V2][V4] = G[V4][V2] = 1;
			C[i] = V3;
			R[j] = V2;
		}
		else
		{
		// 50% Shuffles to 1-4 3-2	
			if (G[V1][V4] || G[V2][V3])
				continue;		// If 1-4 or 2-3 already connects, this shuffling effort aborted
			G[V1][V2] = G[V2][V1] = 0;
			G[V3][V4] = G[V4][V3] = 0;
			G[V1][V4] = G[V4][V1] = 1;
			G[V2][V3] = G[V3][V2] = 1;
			C[i] = V4;
			C[j] = V2;
		}
		Shuffle_times++;
	}
	
	return;
}

void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength)
{
	int i = 0, j = 0, k = 0;
	int N = Rlength - 1;
	bool ** G = new bool * [N];
	for (i = 0; i < N; i++)
	{
		G[i] = new bool [N];
		memset(G[i], 0, sizeof(bool) * (N));
	}
	for (i = 0; i < N; i++)
		for (j = R_src[i]; j < R_src[i+1]; j++)
			G[i][C_src[j]] = 1;

	// Convert to RC format (R1, C1)
	int M = Clength / 2;
	int *R1 = new int [M];
	int *C1 = new int [M];
	for (i = 0; i < N; i++)
		for (j = i+1; j < N; j++)
			if (G[i][j])
			{
				R1[k] = i;
				C1[k] = j;
				k++;
			}
	
	Rewire(G, R1, C1, N, M);
	delete []R1;
	delete []C1;

	// Generate new CSR
	k = 0;
	for (i = 0; i < N; i++)	
	{
		R_dst[i] = R_src[i];
		for (j = 0; j < N; j++)
			if (G[i][j])
				C_dst[k++] = j;			// In fact, R stays unchanged. Update C
	}
	R_dst[N] = R_src[N];
	for (i = N-1; i >= 0; i--)
		delete []G[i];
	delete []G;
	return;
	
}
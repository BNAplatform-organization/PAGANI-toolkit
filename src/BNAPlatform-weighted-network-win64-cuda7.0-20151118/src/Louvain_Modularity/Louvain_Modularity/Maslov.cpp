#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <ctime>
//#include <iostream>
using namespace std;

struct edge
{
     int R;
     int C;
	 float V;
};
typedef struct edge EDGE;

void randperm(EDGE * E, int M)
{	
	srand(time(NULL));
	for (int i=M-1; i>0; --i)
	{  
		int idx = (rand() * (0x7fff+1) + rand()) % (i+1);  
		if (idx == i) continue;  
		float tmp=E[idx].V;  
		E[idx].V = E[i].V;  
		E[i].V = tmp;
	}
}


int cmp1(const void * a, const void * b)
{
     EDGE *aa=(EDGE *)a;
     EDGE *bb=(EDGE *)b;
     if(aa->R!=bb->R)
         return(((aa->R)>(bb->R))?1:-1);
     else
         return(((aa->C)>(bb->C))?1:-1);
}

int cmp2(const void * a, const void * b)
{
     EDGE *aa=(EDGE *)a;
     EDGE *bb=(EDGE *)b;
     if(aa->C!=bb->C)
         return(((aa->C)>(bb->C))?1:-1);
     else
         return(((aa->R)>(bb->R))?1:-1);
}

void Rewire(bool ** G, EDGE * E, int N, int M)
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
		V1 = E[i].R;	
		V2 = E[i].C;
		j = (rand() * (0x7fff+1) + rand()) % M;
		V3 = E[j].R;	
		V4 = E[j].C;				// Select two edges 1-2 3-4
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
			E[i].C = V3;
			E[j].R = V2;
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
			E[i].C = V4;
			E[j].C = V2;
		}
		Shuffle_times++;
	}
	
	return;
}


void Maslov_weighted_1(int * R_dst, int * C_dst, float * V_dst, int * R_src, int * C_src, float * V_src, int Rlength, int Clength)
{
	int i = 0, j = 0, k = 0;
	int N = Rlength - 1;
	bool ** G = new bool * [N];
	for (i = 0; i < N; i++)
	{
		G[i] = new bool [N];
		memset(G[i], 0, sizeof(bool) * (N));
	}
	
	int M = Clength / 2;
	float * vu, * vl;
	vu = new float [M];
	vl = new float [M];
	EDGE * E;
	E = (EDGE *) malloc(M*sizeof(EDGE));
	
	for (i = 0; i < N; i++)
		for (j = R_src[i]; j < R_src[i+1]; j++)
		{
			G[i][C_src[j]] = 1;
			if (C_src[j] > i)
			{
				E[k].R = i;
				E[k].C = C_src[j];
				E[k].V = V_src[j];
				k++;
			}
		}
	// Convert to RC format (R1, C1)
	
	
	//int *R1 = new int [M];
	//int *C1 = new int [M];	
	
	
	Rewire(G, E, N, M);
	//printf("size: %d\n",sizeof(EDGE));
	for (i = 0; i < M; i++)
		if(E[i].R>E[i].C)
		{	
			int tmp = E[i].R;
			E[i].R = E[i].C;
			E[i].C = tmp;
		}
	qsort (E,M,sizeof(EDGE),cmp1);
	for (i = 0; i < M; i++)
	{		
		vu[i] = E[i].V;
	}
	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
	qsort (E,M,sizeof(EDGE),cmp2);
	for (i = 0; i < M; i++)
	{		
		vl[i] = E[i].V;
	}
	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
	free(E);
	//delete []R1;
	//delete []C1;
	
	
	// Generate new CSR
	k = 0;
	int u = 0;
	int l = 0;
	for (i = 0; i < N; i++)	
	{
		R_dst[i] = R_src[i];
		for (j = 0; j < N; j++)
			if (G[i][j])
			{
				C_dst[k] = j;			// In fact, R stays unchanged. Update C
				if (i < j)
					V_dst[k] = vu[u++];
				else
					V_dst[k] = vl[l++];
				k++;
			}
	}
	R_dst[N] = R_src[N];
	for (i = N-1; i >= 0; i--)
		delete []G[i];
	delete []G;

	delete []vu;
	delete []vl;
	return;
	
}


void Maslov_weighted_2(int * R_dst, int * C_dst, float * V_dst, int * R_src, int * C_src, float * V_src, int Rlength, int Clength)
{
	int i = 0, j = 0, k = 0;
	int N = Rlength - 1;
	bool ** G = new bool * [N];
	for (i = 0; i < N; i++)
	{
		G[i] = new bool [N];
		memset(G[i], 0, sizeof(bool) * (N));
	}
	
	int M = Clength / 2;
	float * vu, * vl;
	vu = new float [M];
	vl = new float [M];
	EDGE * E;
	E = (EDGE *) malloc(M*sizeof(EDGE));
	
	for (i = 0; i < N; i++)
		for (j = R_src[i]; j < R_src[i+1]; j++)
		{
			G[i][C_src[j]] = 1;
			if (C_src[j] > i)
			{
				E[k].R = i;
				E[k].C = C_src[j];
				E[k].V = V_src[j];
				k++;
			}
		}
	// Convert to RC format (R1, C1)
	
	
	//int *R1 = new int [M];
	//int *C1 = new int [M];	
	
	
	Rewire(G, E, N, M);
	randperm(E, M);
	printf("size: %d\n",sizeof(EDGE));
	for (i = 0; i < M; i++)
		if(E[i].R>E[i].C)
		{	
			int tmp = E[i].R;
			E[i].R = E[i].C;
			E[i].C = tmp;
		}
	qsort (E,M,sizeof(EDGE),cmp1);
	for (i = 0; i < M; i++)
	{		
		vu[i] = E[i].V;
	}
	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
	qsort (E,M,sizeof(EDGE),cmp2);
	for (i = 0; i < M; i++)
	{		
		vl[i] = E[i].V;
	}
	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
	free(E);
	//delete []R1;
	//delete []C1;
	
	
	// Generate new CSR
	k = 0;
	int u = 0;
	int l = 0;
	for (i = 0; i < N; i++)	
	{
		R_dst[i] = R_src[i];
		for (j = 0; j < N; j++)
			if (G[i][j])
			{
				C_dst[k] = j;			// In fact, R stays unchanged. Update C
				if (i < j)
					V_dst[k] = vu[u++];
				else
					V_dst[k] = vl[l++];
				k++;
			}
	}
	R_dst[N] = R_src[N];
	for (i = N-1; i >= 0; i--)
		delete []G[i];
	delete []G;

	delete []vu;
	delete []vl;
	return;
	
}
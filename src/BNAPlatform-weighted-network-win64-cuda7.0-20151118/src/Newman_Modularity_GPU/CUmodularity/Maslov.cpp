#include "Maslov.h"
#include "bitmap.h"
#include "hashset.h"
#include "hybridmap.h"
//#include "data_type.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

//typedef HashSet<_ULonglong> HashGraph;

const double loadfactor_inv = 2;
//typedef struct KeyValue<_ULonglong, _ULonglong> KeyVal;

//typedef unsigned int u_int;

//struct edge
//{
//     unsigned int R;
//     unsigned int C;
//	 float V;
//};
//typedef struct edge EDGE;


struct Edge {
	u_int x;
	u_int y;
	V_type v;
	Edge(u_int i = 0, u_int j = 0) : x(i), y(j) {};
	bool operator < (const Edge& objstruct) const
	{
		return (x < objstruct.x || (x == objstruct.x && y < objstruct.y));
	}
	bool operator == (const Edge& objstruct) const
	{
		return (x == objstruct.x && y == objstruct.y);
	}
};

void randperm(vector<Edge> &E, int M)
{	
	srand(time(NULL));
	for (int i=M-1; i>0; --i)
	{  
		int idx = (rand() * (0x7fff+1) + rand()) % (i+1);  
		if (idx == i) continue;  
		float tmp=E[idx].v;  
		E[idx].v = E[i].v;  
		E[i].v = tmp;
	}
}


//int cmp1(const void * a, const void * b)
//{
//     EDGE *aa=(EDGE *)a;
//     EDGE *bb=(EDGE *)b;
//     if(aa->R!=bb->R)
//         return(((aa->R)>(bb->R))?1:-1);
//     else
//         return(((aa->C)>(bb->C))?1:-1);
//}
//
//int cmp2(const void * a, const void * b)
//{
//     EDGE *aa=(EDGE *)a;
//     EDGE *bb=(EDGE *)b;
//     if(aa->C!=bb->C)
//         return(((aa->C)>(bb->C))?1:-1);
//     else
//         return(((aa->R)>(bb->R))?1:-1);
//}

//void Rewire(bool ** G, Edge * E, int N, int M)
///*
//input:
//	G :		The adjacency matrix
//	R & C :	Rows and columns of each edge
//	N :		Total number of nodes
//	M :		Total number of edges	
//*/
//
//{
//	int V1, V2, V3, V4;			// Four nodes
//	int Try = 2 * M;			// Total shuffling times
//	srand(time(0));
//	int i = 0, j = 0, k = 0;
//	int Shuffle_times = 0;
//	for (k = 0; k < Try; k++)
//	{
//		i = (rand() * (0x7fff+1) + rand()) % M;	// rand() returns just from 0 to 0x7fff on windows. Too narrow for M
//		V1 = E[i].R;	
//		V2 = E[i].C;
//		j = (rand() * (0x7fff+1) + rand()) % M;
//		V3 = E[j].R;	
//		V4 = E[j].C;				// Select two edges 1-2 3-4
//		if (V1 == V3 || V1 == V4 || V2 == V3 || V2 == V4)
//			continue;		
//		// If two edges have a node in common, this shuffling effort is aborted
//		if (rand() % 2)
//		{	
//		// 50% Shuffles to 1-3 2-4	
//			if (G[V1][V3] || G[V2][V4])
//				continue;		// If 1-3 or 2-4 already connects, this shuffling effort is aborted
//			G[V1][V2] = G[V2][V1] = 0;
//			G[V3][V4] = G[V4][V3] = 0;
//			G[V1][V3] = G[V3][V1] = 1;
//			G[V2][V4] = G[V4][V2] = 1;
//			E[i].C = V3;
//			E[j].R = V2;
//		}
//		else
//		{
//		// 50% Shuffles to 1-4 3-2	
//			if (G[V1][V4] || G[V2][V3])
//				continue;		// If 1-4 or 2-3 already connects, this shuffling effort aborted
//			G[V1][V2] = G[V2][V1] = 0;
//			G[V3][V4] = G[V4][V3] = 0;
//			G[V1][V4] = G[V4][V1] = 1;
//			G[V2][V3] = G[V3][V2] = 1;
//			E[i].C = V4;
//			E[j].C = V2;
//		}
//		Shuffle_times++;
//	}
//	
//	return;
//}

void Rewire(HybridMap &G, vector<Edge> &G_coo, _ULonglong M)
/*
input:
G :		The adjacency matrix
R & C :	Rows and columns of each edge
N :		Total number of nodes
M :		Total number of edges
*/
{
	size_t V1, V2, V3, V4;			// Four nodes
	_ULonglong Try = 2 * M;			// Total shuffling times
	srand(time(0));
	_ULonglong i = 0, j = 0;
	_ULonglong k = 0;
	int Shuffle_times = 0;
	for (k = 0; k < Try; k++)
	{
		i = (rand() * (0x7fff + 1) + rand()) % M;	// rand() returns just from 0 to 0x7fff on windows. Too narrow for M
		V1 = G_coo[i].x;
		V2 = G_coo[i].y;
		j = (rand() * (0x7fff + 1) + rand()) % M;
		V3 = G_coo[j].x;
		V4 = G_coo[j].y;				// Select two edges 1-2 3-4
		if (V1 == V3 || V1 == V4 || V2 == V3 || V2 == V4)
			continue;
		// If two edges have a node in common, this shuffling effort is aborted
		if (rand() % 2)
		{
			// 50% Shuffles to 1-3 2-4	
			//if (G[V1][V3] || G[V2][V4])
			if (G.Get(V1, V3) || G.Get(V2, V4))
				continue;		// If 1-3 or 2-4 already connects, this shuffling effort is aborted
			G.Del(V1, V2);
			//G.Del(V2, V1);
			G.Del(V3, V4);
			//G.Del(V4, V3);

			G.Set(V1, V3);
			//G.Set(V3, V1);
			G.Set(V2, V4);
			//G.Set(V4, V2);

			/*G[V1][V2] = G[V2][V1] = 0;
			G[V3][V4] = G[V4][V3] = 0;
			G[V1][V3] = G[V3][V1] = 1;
			G[V2][V4] = G[V4][V2] = 1;*/

			G_coo[i].y = V3;
			G_coo[j].x = V2;
		}
		else
		{
			// 50% Shuffles to 1-4 3-2	
			//if (G[V1][V4] || G[V2][V3])
			if (G.Get(V1, V4) || G.Get(V2, V3))
				continue;		// If 1-4 or 2-3 already connects, this shuffling effort aborted
			G.Del(V1, V2);
			//G.Del(V2, V1);
			G.Del(V3, V4);
			//G.Del(V4, V3);

			G.Set(V1, V4);
			//G.Set(V4, V1);
			G.Set(V2, V3);
			//G.Set(V3, V2);

			/*G[V1][V2] = G[V2][V1] = 0;
			G[V3][V4] = G[V4][V3] = 0;
			G[V1][V4] = G[V4][V1] = 1;
			G[V2][V3] = G[V3][V2] = 1;*/

			G_coo[i].y = V4;
			G_coo[j].y = V2;
		}
		Shuffle_times++;
	}
	//cout<<Shuffle_times<<endl;
	return;
}

//void Maslov_weighted_1(R_type * R_dst, C_type * C_dst, V_type * V_dst, R_type * R_src, C_type * C_src, V_type * V_src, unsigned int Rlength, R_type Clength)
//{
//	int i = 0, j = 0, k = 0;
//	int N = Rlength - 1;
//	bool ** G = new bool * [N];
//	for (i = 0; i < N; i++)
//	{
//		G[i] = new bool [N];
//		memset(G[i], 0, sizeof(bool) * (N));
//	}
//	
//	int M = Clength / 2;
//	float * vu, * vl;
//	vu = new float [M];
//	vl = new float [M];
//	EDGE * E;
//	E = (EDGE *) malloc(M*sizeof(EDGE));
//	
//	for (i = 0; i < N; i++)
//		for (j = R_src[i]; j < R_src[i+1]; j++)
//		{
//			G[i][C_src[j]] = 1;
//			if (C_src[j] > i)
//			{
//				E[k].R = i;
//				E[k].C = C_src[j];
//				E[k].V = V_src[j];
//				k++;
//			}
//		}
//	// Convert to RC format (R1, C1)
//	
//	
//	//int *R1 = new int [M];
//	//int *C1 = new int [M];	
//	
//	
//	Rewire(G, E, N, M);
//	//printf("size: %d\n",sizeof(EDGE));
//	for (i = 0; i < M; i++)
//		if(E[i].R>E[i].C)
//		{	
//			int tmp = E[i].R;
//			E[i].R = E[i].C;
//			E[i].C = tmp;
//		}
//	qsort (E,M,sizeof(EDGE),cmp1);
//	for (i = 0; i < M; i++)
//	{		
//		vu[i] = E[i].V;
//	}
//	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
//	qsort (E,M,sizeof(EDGE),cmp2);
//	for (i = 0; i < M; i++)
//	{		
//		vl[i] = E[i].V;
//	}
//	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
//	free(E);
//	//delete []R1;
//	//delete []C1;
//	
//	
//	// Generate new CSR
//	k = 0;
//	int u = 0;
//	int l = 0;
//	for (i = 0; i < N; i++)	
//	{
//		R_dst[i] = R_src[i];
//		for (j = 0; j < N; j++)
//			if (G[i][j])
//			{
//				C_dst[k] = j;			// In fact, R stays unchanged. Update C
//				if (i < j)
//					V_dst[k] = vu[u++];
//				else
//					V_dst[k] = vl[l++];
//				k++;
//			}
//	}
//	R_dst[N] = R_src[N];
//	for (i = N-1; i >= 0; i--)
//		delete []G[i];
//	delete []G;
//
//	delete []vu;
//	delete []vl;
//	return;
//	
//}

void Maslov_weighted_1(int * R_dst, int * C_dst, float * V_dst, int * R_src, int * C_src, float * V_src, unsigned int Rlength, int Clength)
{
	size_t i = 0, j = 0;
	size_t k = 0;
	int N = Rlength - 1;

	size_t Gsize = ((size_t)N * N - 1) / 8 + 1;
	//cout<<Gsize<<endl;
	clock_t time_t = clock();
	int M = Clength / 2;
	cout << M << endl;


	vector<Edge> G_coo(M);

	//HashGraph G(Clength,loadfactor_inv);
	HybridMap G(R_src, C_src, N, Clength, loadfactor_inv);
	k = 0;
	for (i = 0; i < N; i++) {
		for (j = R_src[i]; j < R_src[i + 1]; j++)
		{
			if (C_src[j]>i) {
				//_ULonglong key = GetKey(i,C_src[j]);
				//if(!G.Insert(key)) {
				//	cout<<"hashtable insert failure!";
				//}			
				G_coo[k].x = i;
				G_coo[k].y = C_src[j];
				G_coo[k].v = V_src[j];
				++k;
			}
		}
	}
	/*cout<<G.size()<<' '<<G_coo.size()<<endl;
	for(size_t i = 0; i < G_coo.size(); i++)
	{
	_ULonglong key = GetKey(G_coo[i].x, G_coo[i].y);
	if(!G.exist(key))
	cout<<"hash and COO not match: "<<G_coo[i].x<<' '<<G_coo[i].y<<endl;
	}*/

	/*BitMap BG(Gsize);
	for(i = 0; i < N; i++)
	for(j = R_src[i]; j < R_src[i+1]; j++)
	BG.bitmapSet(i*N+C_src[j]);*/



	// Convert to RC format (R1, C1)	

	//cout<<k;
	Rewire(G, G_coo, M);   //hybridMap 

						   //BG.clear();

						   //Rewire(G,G_coo);    //hashgraph
	assert(G_coo.size() * 2 == Clength);
	//G_coo.reserve(Clength);
	for (auto & elem : G_coo) {
		if (elem.x < elem.y)
			swap(elem.x, elem.y);
	}
	sort(G_coo.begin(), G_coo.end());


	//delete []R1;
	//delete []C1;

	// Generate new CSR
	//k = 0;

	time_t = clock() - time_t;
	cout << "rewire time: " << time_t / 1000.0 << endl;
	/*vector<_ULonglong> vec;
	vec.reserve(2*M);
	G.GetValVector(vec);
	cout<<vec.size()<<endl;
	for (size_t i = 0; i < M; i++){
	_ULonglong temp = vec[i]<<32;
	temp = temp | (vec[i]>>32);
	vec.push_back(temp);
	}*/

	//sort(vec.begin(),vec.end());

	//for (i = 0; i < vec.size(); i++)	
	//{
	//	R_dst[i] = R_src[i];
	//	for (j = 0; j < N; j++)
	//	{	//if (G[i][j])
	//		//	C_dst[k++] = j;
	//		_ULonglong key = GetKey(i,j);
	//		if (G.Find(key,key) != -1)	
	//			C_dst[k++] = j;			// In fact, R stays unchanged. Update C				
	//	}
	//}

	memcpy(R_dst, R_src, (N + 1) * sizeof(R_type));
	//cout<<G.size()<<' '<<G_coo.size()<<endl;
	//cout<<G.total_count<<' '<<G.find_cnt<<endl;
	//vector<u_int> deg(N,0);

	for (i = 0; i < G_coo.size(); i++)
	{
		if (i>0) assert(!(G_coo[i] == G_coo[i - 1]));
		u_int x = G_coo[i].x;
		u_int y = G_coo[i].y;
		float v = G_coo[i].v;
		assert(x != y);
		//++deg[x];
		//++deg[y];
		C_dst[R_src[x]] = y;
		C_dst[R_src[y]] = x;
		V_dst[R_src[x]++] = v;
		V_dst[R_src[y]++] = v;
	}
#ifdef DEBUG
	for(i = 0; i < N; i++){
		if(deg[i]!=R_dst[i+1]-R_dst[i])
			cout<<"The results not match! "<<i<<"  "<<deg[i]<<"  "<<R_dst[i+1]-R_dst[i]<<endl;;
		if(R_src[i]!=R_dst[i+1])
			cout<<"The results not match! "<<R_src[i]<<"  "<<R_dst[i+1]<<endl;
	}
#endif
	memcpy(R_src, R_dst, (N + 1) * sizeof(R_type));

	G_coo.clear();
	G_coo.swap(vector<Edge>(0));

	//for (i = N-1; i >= 0; i--)
	//	delete []G[i];
	//delete []G;
	return;
}

void Maslov_weighted_2(int * R_dst, int * C_dst, float * V_dst, int * R_src, int * C_src, float * V_src, unsigned int Rlength, int Clength)
{
	size_t i = 0, j = 0;
	size_t k = 0;
	int N = Rlength - 1;

	size_t Gsize = ((size_t)N * N - 1) / 8 + 1;
	//cout<<Gsize<<endl;
	clock_t time_t = clock();
	int M = Clength / 2;
	cout << M << endl;
	
	vector<Edge> G_coo(M);

	//HashGraph G(Clength,loadfactor_inv);
	HybridMap G(R_src, C_src, N, Clength, loadfactor_inv);
	k = 0;
	for (i = 0; i < N; i++) {
		for (j = R_src[i]; j < R_src[i + 1]; j++)
		{
			if (C_src[j]>i) {
				//_ULonglong key = GetKey(i,C_src[j]);
				//if(!G.Insert(key)) {
				//	cout<<"hashtable insert failure!";
				//}			
				G_coo[k].x = i;
				G_coo[k].y = C_src[j];
				G_coo[k].v = V_src[j];
				++k;
			}
		}
	}
	/*cout<<G.size()<<' '<<G_coo.size()<<endl;
	for(size_t i = 0; i < G_coo.size(); i++)
	{
	_ULonglong key = GetKey(G_coo[i].x, G_coo[i].y);
	if(!G.exist(key))
	cout<<"hash and COO not match: "<<G_coo[i].x<<' '<<G_coo[i].y<<endl;
	}*/

	/*BitMap BG(Gsize);
	for(i = 0; i < N; i++)
	for(j = R_src[i]; j < R_src[i+1]; j++)
	BG.bitmapSet(i*N+C_src[j]);*/



	// Convert to RC format (R1, C1)	

	//cout<<k;
	Rewire(G, G_coo, M);   //hybridMap 
	//Rewire(G,G_coo);    //hashgraph
	randperm(G_coo, M);

	assert(G_coo.size() * 2 == Clength);
	//G_coo.reserve(Clength);
	for (auto & elem : G_coo) {
		if (elem.x < elem.y)
			swap(elem.x, elem.y);
	}
	sort(G_coo.begin(), G_coo.end());


	//delete []R1;
	//delete []C1;

	// Generate new CSR
	//k = 0;

	time_t = clock() - time_t;
	cout << "rewire time: " << time_t / 1000.0 << endl;
	/*vector<_ULonglong> vec;
	vec.reserve(2*M);
	G.GetValVector(vec);
	cout<<vec.size()<<endl;
	for (size_t i = 0; i < M; i++){
	_ULonglong temp = vec[i]<<32;
	temp = temp | (vec[i]>>32);
	vec.push_back(temp);
	}*/

	//sort(vec.begin(),vec.end());

	//for (i = 0; i < vec.size(); i++)	
	//{
	//	R_dst[i] = R_src[i];
	//	for (j = 0; j < N; j++)
	//	{	//if (G[i][j])
	//		//	C_dst[k++] = j;
	//		_ULonglong key = GetKey(i,j);
	//		if (G.Find(key,key) != -1)	
	//			C_dst[k++] = j;			// In fact, R stays unchanged. Update C				
	//	}
	//}

	memcpy(R_dst, R_src, (N + 1) * sizeof(int));
	//cout<<G.size()<<' '<<G_coo.size()<<endl;
	//cout<<G.total_count<<' '<<G.find_cnt<<endl;
	//vector<u_int> deg(N,0);

	for (i = 0; i < G_coo.size(); i++)
	{
		if (i>0) assert(!(G_coo[i] == G_coo[i - 1]));
		u_int x = G_coo[i].x;
		u_int y = G_coo[i].y;
		float v = G_coo[i].v;
		assert(x != y);
		//++deg[x];
		//++deg[y];
		C_dst[R_src[x]] = y;
		C_dst[R_src[y]] = x;
		V_dst[R_src[x]++] = v;
		V_dst[R_src[y]++] = v;
	}
#ifdef DEBUG
	for (i = 0; i < N; i++) {
		if (deg[i] != R_dst[i + 1] - R_dst[i])
			cout << "The results not match! " << i << "  " << deg[i] << "  " << R_dst[i + 1] - R_dst[i] << endl;;
		if (R_src[i] != R_dst[i + 1])
			cout << "The results not match! " << R_src[i] << "  " << R_dst[i + 1] << endl;
	}
#endif
	memcpy(R_src, R_dst, (N + 1) * sizeof(int));

	G_coo.clear();
	G_coo.swap(vector<Edge>(0));

	//for (i = N-1; i >= 0; i--)
	//	delete []G[i];
	//delete []G;
	return;
}

//void Maslov_weighted_2(R_type * R_dst, C_type * C_dst, V_type * V_dst, R_type * R_src, C_type * C_src, V_type * V_src, unsigned int Rlength, R_type Clength)
//{
//	int i = 0, j = 0, k = 0;
//	int N = Rlength - 1;
//	bool ** G = new bool * [N];
//	for (i = 0; i < N; i++)
//	{
//		G[i] = new bool [N];
//		memset(G[i], 0, sizeof(bool) * (N));
//	}
//	
//	int M = Clength / 2;
//	float * vu, * vl;
//	vu = new float [M];
//	vl = new float [M];
//	EDGE * E;
//	E = (EDGE *) malloc(M*sizeof(EDGE));
//	
//	for (i = 0; i < N; i++)
//		for (j = R_src[i]; j < R_src[i+1]; j++)
//		{
//			G[i][C_src[j]] = 1;
//			if (C_src[j] > i)
//			{
//				E[k].R = i;
//				E[k].C = C_src[j];
//				E[k].V = V_src[j];
//				k++;
//			}
//		}
//	// Convert to RC format (R1, C1)
//	
//	
//	//int *R1 = new int [M];
//	//int *C1 = new int [M];	
//	
//	
//	Rewire(G, E, N, M);
//	randperm(E, M);
//	printf("size: %d\n",sizeof(EDGE));
//	for (i = 0; i < M; i++)
//		if(E[i].R>E[i].C)
//		{	
//			int tmp = E[i].R;
//			E[i].R = E[i].C;
//			E[i].C = tmp;
//		}
//	qsort (E,M,sizeof(EDGE),cmp1);
//	for (i = 0; i < M; i++)
//	{		
//		vu[i] = E[i].V;
//	}
//	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
//	qsort (E,M,sizeof(EDGE),cmp2);
//	for (i = 0; i < M; i++)
//	{		
//		vl[i] = E[i].V;
//	}
//	//printf("edge 0 : (%d,%d) %f\n",E[0].R,E[0].C,E[0].V);
//	free(E);
//	//delete []R1;
//	//delete []C1;
//	
//	
//	// Generate new CSR
//	k = 0;
//	int u = 0;
//	int l = 0;
//	for (i = 0; i < N; i++)	
//	{
//		R_dst[i] = R_src[i];
//		for (j = 0; j < N; j++)
//			if (G[i][j])
//			{
//				C_dst[k] = j;			// In fact, R stays unchanged. Update C
//				if (i < j)
//					V_dst[k] = vu[u++];
//				else
//					V_dst[k] = vl[l++];
//				k++;
//			}
//	}
//	R_dst[N] = R_src[N];
//	for (i = N-1; i >= 0; i--)
//		delete []G[i];
//	delete []G;
//
//	delete []vu;
//	delete []vl;
//	return;
//	
//}
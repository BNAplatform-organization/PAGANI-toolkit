#include "bitmap.h"
#include "data_type.h"
//#include "hashtable.h"
#include "hashset.h"
#include "hybridmap.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <memory.h>
#include <ctime>
#include <vector>
#include <algorithm>
//#include <iostream>
using namespace std;

typedef HashSet<_ULonglong> HashGraph;
const double loadfactor_inv = 2;
//typedef struct KeyValue<_ULonglong, _ULonglong> KeyVal;

//typedef unsigned int u_int;

struct Edge{
	u_int x;
	u_int y;
	Edge(u_int i = 0, u_int j = 0): x(i), y(j){};
	bool operator < (const Edge& objstruct) const
	{
		return (x < objstruct.x || (x == objstruct.x && y < objstruct.y));		
	}
	bool operator == (const Edge& objstruct) const
	{
		return (x == objstruct.x && y == objstruct.y);		
	}
}; 



//_ULonglong GetVal(_ULonglong v1, _ULonglong v2)
//{
//	return ((v1<<32) | v2);
//}

void Rewire(BitMap &G, vector<Edge> &G_coo, int N, int M)
/*
input:
	G :		The adjacency matrix
	R & C :	Rows and columns of each edge
	N :		Total number of nodes
	M :		Total number of edges	
*/
{
	size_t V1, V2, V3, V4;			// Four nodes
	int Try = 2 * M;			// Total shuffling times
	srand(time(0));
	int i = 0, j = 0, k = 0;
	int Shuffle_times = 0;
	for (k = 0; k < Try; k++)
	{
		i = (rand() * (0x7fff+1) + rand()) % M;	// rand() returns just from 0 to 0x7fff on windows. Too narrow for M
		V1 = G_coo[i].x;	
		V2 = G_coo[i].y;
		j = (rand() * (0x7fff+1) + rand()) % M;
		V3 = G_coo[j].x;	
		V4 = G_coo[j].y;				// Select two edges 1-2 3-4
		if (V1 == V3 || V1 == V4 || V2 == V3 || V2 == V4)
			continue;		
		// If two edges have a node in common, this shuffling effort is aborted
		if (rand() % 2)
		{	
		// 50% Shuffles to 1-3 2-4	
			//if (G[V1][V3] || G[V2][V4])
			if ( G.bitmapGet(V1*N + V3) || G.bitmapGet(V2*N + V4) )
				continue;		// If 1-3 or 2-4 already connects, this shuffling effort is aborted
			G.bitmapDel(V1*N + V2);
			G.bitmapDel(V2*N + V1);
			G.bitmapDel(V3*N + V4);
			G.bitmapDel(V4*N + V3);

			G.bitmapSet(V1*N + V3);
			G.bitmapSet(V3*N + V1);
			G.bitmapSet(V2*N + V4);
			G.bitmapSet(V4*N + V2);

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
			if ( G.bitmapGet(V1*N + V4) || G.bitmapGet(V2*N + V3) )
				continue;		// If 1-4 or 2-3 already connects, this shuffling effort aborted
			G.bitmapDel(V1*N + V2);
			G.bitmapDel(V2*N + V1);
			G.bitmapDel(V3*N + V4);
			G.bitmapDel(V4*N + V3);

			G.bitmapSet(V1*N + V4);
			G.bitmapSet(V4*N + V1);
			G.bitmapSet(V2*N + V3);
			G.bitmapSet(V3*N + V2);
			
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
		i = (rand() * (0x7fff+1) + rand()) % M;	// rand() returns just from 0 to 0x7fff on windows. Too narrow for M
		V1 = G_coo[i].x;	
		V2 = G_coo[i].y;
		j = (rand() * (0x7fff+1) + rand()) % M;
		V3 = G_coo[j].x;	
		V4 = G_coo[j].y;				// Select two edges 1-2 3-4
		if (V1 == V3 || V1 == V4 || V2 == V3 || V2 == V4)
			continue;		
		// If two edges have a node in common, this shuffling effort is aborted
		if (rand() % 2)
		{	
		// 50% Shuffles to 1-3 2-4	
			//if (G[V1][V3] || G[V2][V4])
			if ( G.Get(V1, V3) || G.Get(V2, V4) )
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
			if ( G.Get(V1, V4) || G.Get(V2, V3) )
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

void Rewire(HashGraph &G, vector<Edge> &G_coo)
/*
input:
	G :		The adjacency matrix stored in hashgraph
	R & C :	Rows and columns of each edge
	N :		Total number of nodes
	M :		Total number of edges	
*/
{
	u_int V1, V2, V3, V4;			// Four nodes
	size_t M = G.size();
	size_t Try = 2 * M;			// Total shuffling times
	srand(1);
	size_t i = 0, j = 0, k = 0;
	size_t Shuffle_times = 0;
	clock_t time = clock();
	for (k = 0; k < Try; ++k)
	{
		/*for(size_t debug_i = 0; debug_i < G_coo.size(); debug_i++)
			{
				_ULonglong key = GetKey(G_coo[debug_i].x, G_coo[debug_i].y);
				
				if(!G.exist(key))
					cout<<Shuffle_times<<" hash and COO not match: "<<G_coo[debug_i].x<<' '<<G_coo[debug_i].y<<endl;
				
			}*/
		//if (k%1000000==0)
		//	cout<<k<<" time: "<<clock()-time<<endl;
		i = (rand() * (0x7fff+1) + rand()) % M;
		//i =  M * ((double) rand()/RAND_MAX) ;	// rand() returns just from 0 to 0x7fff on windows. Too narrow for M
		V1 = G_coo[i].x;	
		V2 = G_coo[i].y;
		_ULonglong k1 = GetKey(V2,V1);
		if( !G.exist(k1) )
			cout<<Shuffle_times<<" error: "<<V1<<" "<<V2<<endl;

		j = (rand() * (0x7fff+1) + rand()) % M;
		//j = M * ((double) rand()/RAND_MAX) ;
		V3 = G_coo[j].x;	
		V4 = G_coo[j].y;			// Select two edges 1-2 3-4
		_ULonglong k2 = GetKey(V4,V3);
		if(!G.exist(k2))
			cout<<Shuffle_times<<" error: "<<V3<<" "<<V4<<endl;

		if (V1 == V3 || V1 == V4 || V2 == V3 || V2 == V4)
			continue;		
		//cout<<V1<<' '<<V2<<' '<<V3<<' '<<V4<<" to be check"<<endl;

		// If two edges have a node in common, this shuffling effort is aborted
		if (rand() % 2)
		{	
		// 50% Shuffles to 1-3 2-4	
			/*if (G[V1][V3] || G[V2][V4])
			continue;*/

			k1 = GetKey(V1,V3);
			k2 = GetKey(V2,V4);
			if ( G.exist(k1) || G.exist(k2))
				continue;		// If 1-3 or 2-4 already connects, this shuffling effort is aborted
			
			//k1 = GetKey(V1,V3);
			if(!G.Insert(k1))
				cout<<V1<<' '<<V3<<endl;
			//debug
			//if(!G.exist(k1))
			//	cout<<"did not insert"<< V1<<' '<<V3<<endl;

			//k2 = GetKey(V2,V4);
			if(!G.Insert(k2) )
				cout<<V2<<' '<<V4<<endl;
			
			//debug
			//if(!G.exist(k2))
			//	cout<<"did not insert"<< V2<<' '<<V4<<endl;

			k1 = GetKey(V1,V2);
			if(!G.Remove(k1))
				cout<<V1<<' '<<V2<<endl;
			//if( G.exist(k1) )
			//	cout<<"did not remove"<< V1<<' '<<V2<<endl;
			
			k2 = GetKey(V3,V4);
			if(!G.Remove(k2))
				cout<<V3<<' '<<V4<<endl;
			//if( G.exist(k2) )
			//	cout<<"did not remove"<< V3<<' '<<V4<<endl;
			
			
			G_coo[i].y = V3;
			G_coo[j].x = V2;	
			//k1 = GetKey(G_coo[i].x,G_coo[i].y);
			//k2 = GetKey(G_coo[j].x,G_coo[j].y);
			
			//cout<<V1<<' '<<V3<<' '<<V2<<' '<<V4<<" shuffled"<<endl;
			/*cout<<"G size"<<G.size()<<endl;
			if( !G.exist(k1) )
				cout<<"did not find"<<G_coo[i].x<<' '<<G_coo[i].y<<' '<<V1<<' '<<V3<<endl;
			if( !G.exist(k2) )
				cout<<"did not find"<<G_coo[i].x<<' '<<G_coo[i].y<<' '<<V2<<' '<<V4<<endl;*/
			
		}
		else
		{
		// 50% Shuffles to 1-4 3-2	
			//if (G[V1][V4] || G[V2][V3])
			//continue;
			
			k1 = GetKey(V1,V4);
			k2 = GetKey(V2,V3);
			if ( G.exist(k1) || G.exist(k2))
				continue;		// If 1-3 or 2-4 already connects, this shuffling effort is aborted
			
			//k1 = GetKey(V1,V3);
			if(!G.Insert(k1))
				cout<<V1<<' '<<V4<<endl;
			//if(!G.exist(k1))
			//	cout<<"did not insert"<< V1<<' '<<V3<<endl;

			//k2 = GetKey(V2,V4);
			if(!G.Insert(k2) )
				cout<<V2<<' '<<V3<<endl;
			//if(!G.exist(k2))
			//	cout<<"did not insert"<< V2<<' '<<V3<<endl;

			k1 = GetKey(V1,V2);
			if(!G.Remove(k1))
				cout<<V1<<' '<<V2<<endl;
			//if( G.exist(k1) )
			//	cout<<"did not remove"<< V1<<' '<<V2<<endl;
			
			k2 = GetKey(V3,V4);
			if(!G.Remove(k2))
				cout<<V3<<' '<<V4<<endl;
			//if( G.exist(k2) )
			//	cout<<"did not remove"<< V3<<' '<<V4<<endl;
			
			G_coo[i].y = V4;
			G_coo[j].y = V2;	
			///*G[V1][V2] = G[V2][V1] = 0;
			//G[V3][V4] = G[V4][V3] = 0;
			//G[V1][V4] = G[V4][V1] = 1;
			//G[V2][V3] = G[V3][V2] = 1;*/

			//C[i] = V4;
			//C[j] = V2;
		}
		Shuffle_times++;
	}
	//cout<<Shuffle_times<<endl;
	return;
}

void Maslov(R_type * R_dst, C_type * C_dst, R_type * R_src, C_type * C_src, C_type Rlength, R_type Clength)
{
	size_t i = 0, j = 0;
	size_t k = 0;
	int N = Rlength - 1;
	
	/*bool ** G = new bool * [N];
	for (i = 0; i < N; i++)
	{
		G[i] = new bool [N];
		memset(G[i], 0, sizeof(bool) * (N));
	}
	for (i = 0; i < N; i++)
		for (j = R_src[i]; j < R_src[i+1]; j++)
			G[i][C_src[j]] = 1;	*/
	
	size_t Gsize = ((size_t) N * N - 1)/8 +1; 
	//cout<<Gsize<<endl;
	clock_t time_t = clock();
	int M = Clength / 2;
	cout<<M<<endl;
	
		
	//int *R1 = new int [M];
	//int *C1 = new int [M];
	
	
	vector<Edge> G_coo(M);
	
	//HashGraph G(Clength,loadfactor_inv);
	HybridMap G(R_src, C_src, N, Clength,loadfactor_inv);	
	k = 0;
	for (i = 0; i < N; i++){
		for (j = R_src[i]; j < R_src[i+1]; j++)
		{
			if(C_src[j]>i){
				//_ULonglong key = GetKey(i,C_src[j]);
				//if(!G.Insert(key)) {
				//	cout<<"hashtable insert failure!";
				//}			
				G_coo[k].x = i;
				G_coo[k].y = C_src[j];
				++k;
			}
		}
	}
	//cout<<G.size()<<' '<<G_coo.size()<<endl;
	/*for(size_t i = 0; i < G_coo.size(); i++)
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
	assert(G_coo.size()*2 == Clength);
	//G_coo.reserve(Clength);
	for(auto & elem : G_coo){
		if(elem.x < elem.y)
			swap(elem.x,elem.y);
	}
	sort(G_coo.begin(),G_coo.end());
	
	//G_coo.clear();
	//G_coo.swap(vector<Edge>(0));
	//delete []R1;
	//delete []C1;

	// Generate new CSR
	//k = 0;

	time_t = clock()-time_t;
	cout<<"rewire time: "<<time_t/1000.0<<endl;
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
	memcpy(R_dst, R_src, (N+1)*sizeof(R_type));
	//cout<<G.size()<<' '<<G_coo.size()<<endl;
	//cout<<G.total_count<<' '<<G.find_cnt<<endl;
	vector<u_int> deg(N,0);

	for(i = 0; i < G_coo.size(); i++)
	{		
		if(i>0) assert(!(G_coo[i] == G_coo[i-1])); 
		u_int x = G_coo[i].x;
		u_int y = G_coo[i].y;
		assert(x!=y);
		++deg[x];
		++deg[y];
		//C_dst[R_src[x]++] = y;
		//C_dst[R_src[y]++] = x;		
	}
	//for(i = 0; i < N; i++){
	//	if(deg[i]!=R_dst[i+1]-R_dst[i]) 
	//		cout<<"The results not match! "<<i<<"  "<<deg[i]<<"  "<<R_dst[i+1]-R_dst[i]<<endl;;
		//if(R_src[i]!=R_dst[i])
		//	cout<<"The results not match! "<<R_src[i]<<"  "<<R_dst[i+1]<<endl;
	//}
	memcpy(R_src, R_dst, (N+1)*sizeof(R_type));
	//for (i = N-1; i >= 0; i--)
	//	delete []G[i];
	//delete []G;
	return;	
}

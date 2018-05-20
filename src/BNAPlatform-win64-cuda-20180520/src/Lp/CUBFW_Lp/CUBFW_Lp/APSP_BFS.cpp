#include <stack>
#include <queue>
#include <list>
#include <algorithm>
#include <Windows.h>
#include <process.h>
#include <iostream>
#include "data_type.h"
using namespace std;

extern float *Li_result;

typedef struct arg
{
	int id;
	unsigned int numVertices;
	int numThread;
	R_type *row;
	C_type *col;
	int *deg;
	int *sort_indices;
	//float *dist;
	float *Li;	
	double *Lpsum;
} ARG;

void APSP_thread(void *i);

//double APSP_BFS(float *APSP_output, int * row, int * col, const int numVertices)
double APSP_BFS(R_type * row, C_type * col, const unsigned int numVertices)
{
	SYSTEM_INFO siSysInfo;
	GetSystemInfo(&siSysInfo); 
	int numThread = (int)siSysInfo.dwNumberOfProcessors/2;
	cout<<numThread<<endl;
	//numThread = 2;
	//int numThread = sysconf(_SC_NPROCESSORS_CONF);
	//int numVertices = dim;
	int 	*deg = new int[numVertices];	// Cormat degree
	for(int i = 0; i < numVertices; i++)
		deg[i] = row[i+1] - row[i];
	
	//long long size = numVertices;
	//size *= numVertices;
	//float *APSP_output = new float[size]; 			//output
	//memset(APSP_output, 0, sizeof(float)*size);

	int* sort_degree = new int[numVertices];
	int* sort_indices = new int[numVertices];
	if( sort_degree == NULL || sort_indices == NULL)
	{
		cerr<<"allocate memory failure!"<<endl;
		exit(1);
	}
	//sort
	for(int i = 0; i < numVertices; i++) 
	{
		sort_indices[i] = i;
		sort_degree[i] = deg[i];
	}
	
	for(int i = 0; i < numVertices; i++)
	{
		int max_index = i;
		for(int j = i; j < numVertices; j++)
		{
			if (sort_degree[max_index] < sort_degree[j])
				max_index = j;
		}
		int temp = sort_degree[max_index];
		sort_degree[max_index] = sort_degree[i];
		sort_degree[i] = temp;
		
		temp = sort_indices[max_index];
		sort_indices[max_index] = sort_indices[i];
		sort_indices[i] = temp;
	}
	//pthread_t *thread = new pthread_t[numThread];
	
	memset(Li_result, 0, sizeof(float)*numVertices);
	//float **APSP_dist = new float *[numThread]; 
	double *Lpsum = new double [numThread];
	//memset(APSP_dist, 0, sizeof(float)*numVertices*numThread);

	HANDLE *tHandle = new HANDLE[numThread];
	ARG *arg = new ARG[numThread];

	
	for (int i = 0; i < numThread; i++)
	{
		
		arg[i].id = i;
		arg[i].numVertices = numVertices;
		arg[i].numThread = numThread;
		arg[i].row = row;
		arg[i].col = col;
		arg[i].deg = deg;
		arg[i].sort_indices = sort_indices;
		//APSP_dist[i] = new float [numVertices];
		//arg[i].dist = APSP_dist[i];
		arg[i].Li = Li_result;
		arg[i].Lpsum = Lpsum+i;
	}
	printf("APSP computing started...\n");
	
	for (int i = 0; i < numThread; i++)
	{
		ARG *temp = arg + i;
		//pthread_create(&thread[i], NULL, APSP, (void*)(temp));
		tHandle[i] = (HANDLE) _beginthread(APSP_thread, 0, (char *)temp);
	}	

	for (int i = 0; i < numThread; i++)
		//pthread_join(thread[i], NULL);
		WaitForSingleObject(tHandle[i], INFINITE);

	/*for (int i = 0; i < numThread; i++)
		delete []APSP_dist[i];
	delete []APSP_dist;*/
	double Lp_result = 0.0;
	//double Lp_check  = 0.0;
	
	
	int  isolated_vertices = 0;
	//bool unconnected_mark = false;
	
	//for (int i = 0; i < numVertices; i++)
	//{
	//	for (int j = i + 1; j < numVertices; j++)
	//	{
	//		if (APSP_output[(long long)i * numVertices + j]<1)
	//			unconnected_mark = true;
	//		//cout<<"ii: "<<ii<<"\tjj: "<<jj<<"\ti: "<<i<<"\tj: "<<j<<"\tindex: "<<index<<"\tcostmatPaded[index]:"<<costmatPaded[index]<<endl;
	//		else
	//		{
	//			Lp_result += 1.0f /APSP_output[(long long)i * numVertices + j];
	//			Li_result[i] -= 1.0f/APSP_output[(long long)i * numVertices + j];
	//			Li_result[j] -= 1.0f/APSP_output[(long long)i * numVertices + j];
	//		}
	//	}
	//	Li_result[i] /= (numVertices - 1);
	//	//if (i > 10000 && i < 10020)
	//		//cout<<Li_result[i]<<endl;
	//	if (Li_result[i] == 0) 
	//	{	
	//		isolated_vertices++; 
	//		//cout<<i<<"th node is isolated;\n"; 
	//	}
	//	//else	Li_result[i] = 1.0/Li_result[i];
	//}
	for (int i = 0; i < numVertices; i++)
		if (Li_result[i] == 0)
			isolated_vertices++;
	
	for (int i = 0; i < numThread; i++)
		Lp_result += Lpsum[i];

	Lp_result /= numVertices-1;
	Lp_result /= numVertices;
	Lp_result = 1 / Lp_result;
	if (isolated_vertices > 0 )
		cout<<"\nisolated vertices number : "<<isolated_vertices<<endl;
	//else if (unconnected_mark)
	//	cout<<"\nThe network is unconnected. "<<endl;
	cout<<"Lp_result: "<<Lp_result<<endl;

	delete []sort_degree;
	delete []sort_indices;
	//delete []thread;
	delete []tHandle;
	delete []arg;
	delete []deg;
	return Lp_result;	
}


void APSP_thread(void *i)
{
	ARG *temp = (ARG*)i;
	int 	id = temp->id;
	unsigned int 	numVertices = temp->numVertices;
	int 	numThread = temp->numThread;
	R_type *row = temp->row;
	C_type *col = temp->col;
	int	*deg = temp->deg;
	int	*sort_indices = temp->sort_indices;
	//float	*APSP = temp->APSP;
	float *Li = temp->Li;
	//float *dist = temp->dist; 
	double *Lpsum = temp->Lpsum;
	float *dist = new float[numVertices];
	double Litemp = 0;
	
	queue<int> Q;	
	
	*Lpsum = 0;
	//int timer_thread = (id+1)*2;
	//int timer_iteration = timer_thread + 1;
	for (int s = id; s < numVertices; s += numThread) 
	{		
		for (int i = 0; i < numVertices; i++) {
			dist[i] = -1.0; 
		}

		dist[sort_indices[s]] = 0;
		Q.push(sort_indices[s]);
		int w;
		while (!Q.empty())
		{
			int v = Q.front();
			Q.pop();			
			for(int j = row[v]; j < row[v+1]; j++)
			{
				w = col[j];
				if (w >= numVertices)
					cout<<"id = "<<id<<"; w = "<<w<<"; j = "<<j<<"; col[j] = "<<col[j]<<endl;
				if (dist[w] < 0)
				{
					dist[w] = dist[v] + 1;
					Q.push(w);
					//APSP[(long long)sort_indices[s]*numVertices+w]=APSP[(long long)sort_indices[s]*numVertices+v]+1;				
				}
			}
		}
		Litemp = 0;
		for ( int i = 0; i < numVertices; i++)
		{
			if (dist[i]<1)
				continue;
			else
				Litemp += 1.0/dist[i];
				//Li_result[sort_indices[s]] += 1.0f/dist[i];
		}
		*Lpsum += Litemp;
		Li_result[sort_indices[s]] = Litemp/(numVertices - 1);
	}	
	
	delete []dist;
	return;
}

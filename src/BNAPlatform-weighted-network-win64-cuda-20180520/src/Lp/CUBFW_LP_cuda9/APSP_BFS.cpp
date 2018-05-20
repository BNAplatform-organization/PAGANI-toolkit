#include <stack>
#include <queue>
#include <list>
#include <algorithm>
#include <Windows.h>
#include <process.h>
#include <iostream>
using namespace std;

typedef struct arg
{
	int id;
	int numVertices;
	int numThread;
	int *row;
	int *col;
	int *deg;
	int *sort_indices;
	float *APSP;
	
} ARG;

void APSP_thread(void *i);

double APSP_BFS(float *APSP_output, int * row, int * col, const int numVertices)
{
	SYSTEM_INFO siSysInfo;
	GetSystemInfo(&siSysInfo); 
	int numThread = (int)siSysInfo.dwNumberOfProcessors*2; 
	//int numThread = sysconf(_SC_NPROCESSORS_CONF);
	//int numThread = 1;
	int 	*deg = new int[numVertices];	// Cormat degree
	for(int i = 0; i < numVertices; i++)
		deg[i] = row[i+1] - row[i];
	
	long long size = numVertices;
	size *= numVertices;
	//float *APSP_output = new float[size]; 			//output
	memset((void*)APSP_output, 0, sizeof(float)*size);

	int* sort_degree = new int[numVertices];
	int* sort_indices = new int[numVertices];
	if(APSP_output == NULL || sort_degree == NULL || sort_indices == NULL)
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
		arg[i].APSP = APSP_output;
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

	double Lp_result = 0.0;
	for (int i = 0; i < numVertices; i++)
	{
		for (int j = i + 1; j < numVertices; j++)
		{
			if (APSP_output[(long long)i * numVertices + j]<1)
				;
			//cout<<"ii: "<<ii<<"\tjj: "<<jj<<"\ti: "<<i<<"\tj: "<<j<<"\tindex: "<<index<<"\tcostmatPaded[index]:"<<costmatPaded[index]<<endl;
			else
				Lp_result += 1.0f /APSP_output[(long long)i * numVertices + j];

		}
	}

	Lp_result /= numVertices;
	Lp_result /= (numVertices - 1);
	Lp_result = 1 / Lp_result;
	


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
	int 	numVertices = temp->numVertices;
	int 	numThread = temp->numThread;
	int	*row = temp->row;
	int	*col = temp->col;
	int	*deg = temp->deg;
	int	*sort_indices = temp->sort_indices;
	float	*APSP = temp->APSP;
	int 	*dist 		= new int[numVertices];

	queue<int> Q;	
	int timer_thread = (id+1)*2;
	int timer_iteration = timer_thread + 1;
	for (int s = id; s < numVertices; s += numThread) {
		for (int j = 0; j < numVertices; j++) {
			dist[j] = -1; 
		}
		dist[sort_indices[s]] = 0;
		Q.push(sort_indices[s]);
		while (!Q.empty())
		{
			int v = Q.front();
			Q.pop();	
			int w;
			for(int i = row[v]; i < row[v+1]; i++)
			{
				w = col[i];
				if (dist[w] < 0)
				{
					dist[w] = dist[v] + 1;
					Q.push(w);
					APSP[(long long)sort_indices[s]*numVertices+w]=APSP[(long long)sort_indices[s]*numVertices+v]+1;
				}
			}
		}
	}	
	
	delete []dist;
	return;
}

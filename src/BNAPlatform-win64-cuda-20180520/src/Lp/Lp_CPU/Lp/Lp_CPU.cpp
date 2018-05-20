/*
** calculate APSP for undirected unweighted network using BFS(dijkstra) algorithm on multicore CPU.
** 2011.7.7 by xumo
** 2011.7.16 modified by xumo 
** transplant to linux OS. Different: how to get the number of CPU cores.
*/
#include "Lp_CPU.h"
inline int diagindex(int n, int i, int j)
{
	return (i<j)?(i*(2*n-1-i)/2+j-i):(j*(2*n-1-j)/2+i-j-1);		
};

double Lp_CPU(int *row, int *col, int numVertices, int numEdges)
{
	SYSTEM_INFO siSysInfo;
	GetSystemInfo(&siSysInfo); 
	int numThread = (int)siSysInfo.dwNumberOfProcessors; 
	//int numThread = sysconf(_SC_NPROCESSORS_CONF);
	
	int 	*deg = new int[numVertices];	// Cormat degree
	for(int i = 0; i < numVertices; i++)
		deg[i] = row[i+1] - row[i];
	
	long long size = numVertices;
	size *= numVertices;
	float *APSP_output = new float[size]; 			//output
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
	Setup(0);	
	Start(0);
	
	for (int i = 0; i < numThread; i++)
	{
		ARG *temp = arg + i;
		//pthread_create(&thread[i], NULL, APSP, (void*)(temp));
		tHandle[i] = (HANDLE) _beginthread(APSP, 0, (char *)temp);
	}	

	for (int i = 0; i < numThread; i++)
		//pthread_join(thread[i], NULL);
		WaitForSingleObject(tHandle[i], INFINITE);
	Stop(0);
	cout<<"APSP calculation completed! Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;
	delete []sort_degree;
	delete []sort_indices;
	//delete []thread;
	delete []tHandle;
	delete []arg;
	delete []deg;
	
	double Lp_result = 0;

	float *temp = APSP_output + 1;	
	
	for (int i = 0; i < numVertices; i++)
		for (int j = i + 1; j < numVertices; j++)
		{
			if (!(fabs(APSP_output[i * numVertices + j])<1e-12))
			{
				//APSP_output[i * numVertices + j] = INFINITE;//__FLT_MAX__;
				Lp_result += 1.0f / APSP_output[i * numVertices + j];
			}
		}
	
	Lp_result /= numVertices;
	Lp_result /= (numVertices-1);
	Lp_result *= 2;
	Lp_result = 1 / Lp_result;

	//ofstream fout("test.txt");
	//for (int i = 0; i < 1000000; i ++)
	//{
	//	fout<<APSP_output[i]<<endl;
	//}
	//fout.close();

	delete []APSP_output;
	
	//string inputfile = string(argv[1]);
	//string APSPoutfile =  inputfile.append(".apsp"); 
	//ofstream fAPSP(APSPoutfile.c_str());
	//if (!fAPSP.good())
	//{
	//	cerr<<"cannot write .apsp file!"<<endl;
	//	exit(1);
	//}
	//cout<<"The APSP is saved as "<<APSPoutfile<<"..."<<endl;
	//float *temp = APSP_output + 1;	
	//
	//for (int i = 0; i < numVertices; i++)
	//	for (int j = i; j < numVertices; j++)
	//		if (fabs(APSP_output[i * numVertices + j])<1e-12)
	//			APSP_output[i * numVertices + j] = INFINITE;//__FLT_MAX__;

	//fAPSP.write((char*)&numVertices, 4);
	//for (int i = 0; i < numVertices; i++)
	//{
	//	fAPSP.write((char*)(temp), sizeof(float) * (numVertices - i - 1));
	//	temp += numVertices + 1;
	//}	
	//fAPSP.close();

	return Lp_result;
}

void APSP(void *i)
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
			for(int i = row[v]; i < row[v+1]; i++)
			{
				int w = col[i];
				if (dist[w] < 0)
				{
					dist[w] = dist[v] + 1;
					Q.push(w);
					APSP[sort_indices[s]*numVertices+w]=APSP[sort_indices[s]*numVertices+v]+1;
				}
			}
		}
	}	

	delete []dist;
	return;
}

# include "data_type.h"
# include <windows.h>
# include <iostream>
# include <process.h>
# include <fstream>
# include <math.h>
using namespace std;
// This struct is used for passing parameters to different threads
struct Cp_ARG
{
	int id;	
	C_type * C;			
	R_type * R;			// R and C represents the network in CSR format
	V_type * V;
	int N;				// The size of the network
	double *	Cpsum;		
	float * Cp;
	float * deg;
};

int THREADNUM = 4;

void Cp_single_thread(void * arg);

double Cp_1(C_type * C, R_type * R, V_type * V, float * Cp, int N)
{
	int i = 0, j = 0;
	float * deg = new float [N];
	memset(deg, 0, sizeof(float) * N);
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
		{	//G[i][C[j]] = V[j];
			deg[i] += V[j];
			//cout<<G[i][C[j]]<<' ';
			//cout<<V[j]<<' ';
	}

	SYSTEM_INFO siSysInfo;
	GetSystemInfo(&siSysInfo); 
	THREADNUM = (int)siSysInfo.dwNumberOfProcessors; 
	//THREADNUM = sysconf(_SC_NPROCESSORS_CONF);	

	Cp_ARG * Properties_arg = new Cp_ARG [THREADNUM];	
	double * Cpsum = new double[THREADNUM];
	for (i = 0; i < THREADNUM; i++)
	{
		Properties_arg[i].id = i;
		Properties_arg[i].C = C;
		Properties_arg[i].R = R;
		Properties_arg[i].V = V;
		Properties_arg[i].N = N;
		Properties_arg[i].Cpsum = Cpsum + i;
		Properties_arg[i].Cp = Cp;
		Properties_arg[i].deg = deg;
	}
	//pthread_t *t = new pthread_t[THREADNUM];
	HANDLE *tHandle = new HANDLE[THREADNUM];
	for (i = 0; i < THREADNUM; i++)
	{
		Cp_ARG *temp = Properties_arg + i;
		//pthread_create(&t[i], NULL, All_Properties_single_thread, (void*)temp);
		tHandle[i] = (HANDLE) _beginthread(Cp_single_thread, 0, (char *)temp);
	}
	for (i = 0; i < THREADNUM; i++)
		WaitForSingleObject(tHandle[i], INFINITE);
		//pthread_join(t[i], NULL);

	double mean_Cp = 0;
	for (i = 0; i < THREADNUM; i++)
	{
		mean_Cp += Cpsum[i];		
	}
	cout<<"mean:"<<mean_Cp/N<<endl;
	
	delete []Properties_arg;
	delete []Cpsum;
	delete []tHandle;
	return mean_Cp/N;
}

void Cp_single_thread(void * voidarg)
{
	Cp_ARG *arg = (Cp_ARG *)voidarg;
	int id = arg->id;
	C_type * C = arg->C;
	R_type * R = arg->R;
	float * V = arg->V;
	int N = arg->N;
	double * Cpsum = arg->Cpsum;
	float * Cp = arg->Cp;
	float * deg = arg->deg;
	
	int i = 0, j = 0, k = 0, ii = 0, jj = 0;
	double Lsum = 0;
	
	int Neighbour_Num = 0;
	double mean_Cp = 0;
	double Ctemp = 0;
	
	for (i = id; i < N; i += THREADNUM)
	{
		
		Neighbour_Num = R[i+1] - R[i];
		if (Neighbour_Num < 2)
		{
			Cp[i] = 0;
			continue;
		}			
		Ctemp = 0;
		for (k = 0; k < Neighbour_Num; k++)
		{
			j = C[R[i]+k];
			int Neighbour_Num_j =  R[j+1] - R[j];
			for (ii = k+1, jj = 0; ii < Neighbour_Num && jj <Neighbour_Num_j; )
			{
				if(C[R[i]+ii] == C[R[j]+jj])
				{	
					Ctemp += (double) V[R[i]+k] + V[R[i]+ii]; 
					++ii; 
					++jj;
				}
				else if(C[R[i]+ii] < C[R[j]+jj]) ++ii;
				else ++jj;				
			}
		}
        Cp[i] = Ctemp/ (Neighbour_Num - 1)/deg[i];
		mean_Cp += Cp[i];
	}
	*Cpsum = (float)mean_Cp;
	return;	
}


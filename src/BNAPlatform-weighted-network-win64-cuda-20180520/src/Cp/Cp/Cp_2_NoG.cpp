# include <windows.h>
# include <iostream>
# include <process.h>
# include <fstream>
# include <math.h>
# include "data_type.h"
using namespace std;
// This struct is used for passing parameters to different threads
struct Cp_ARG_2
{
	int id;
	//float ** G;			// The network
	C_type * C;			
	R_type * R;			// R and C represents the network in CSR format
	V_type * V;
	int N;				// The size of the network
	double *	Cpsum;		
	float * Cp;
};

int THREADNUM_2 = 4;

void Cp_single_thread_2(void * arg);

double Cp_2(C_type * C, R_type * R, V_type * V, float * Cp, int N)
{
	int i = 0, j = 0;
	//float ** G = new float * [N];
    V_type V_max = V[0];
    for (i = 1; i < N; i++)
		V_max = max(V_max,V[i]);

	/*for (i = 0; i < N; i++)
	{
		G[i] = new float [N];
		memset(G[i], 0, sizeof(float) * (N));		
	}*/
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
        {	//G[i][C[j]] = V[j]/V_max;
			V[j] /= V_max; 
			//cout<<G[i][C[j]]<<' ';
			//cout<<V[j]<<' ';
		}

	SYSTEM_INFO siSysInfo;
	GetSystemInfo(&siSysInfo); 
	THREADNUM_2 = (int)siSysInfo.dwNumberOfProcessors; 
	//THREADNUM = sysconf(_SC_NPROCESSORS_CONF);	

	Cp_ARG_2 * Properties_arg = new Cp_ARG_2 [THREADNUM_2];	
	double * Cpsum = new double[THREADNUM_2];
	for (i = 0; i < THREADNUM_2; i++)
	{
		Properties_arg[i].id = i;
		//Properties_arg[i].G = G;
		Properties_arg[i].C = C;
		Properties_arg[i].R = R;
		Properties_arg[i].V = V;
		Properties_arg[i].N = N;
		Properties_arg[i].Cpsum = Cpsum + i;
		Properties_arg[i].Cp = Cp;
	}
	//pthread_t *t = new pthread_t[THREADNUM];
	HANDLE *tHandle = new HANDLE[THREADNUM_2];
	for (i = 0; i < THREADNUM_2; i++)
	{
		Cp_ARG_2 *temp = Properties_arg + i;
		//pthread_create(&t[i], NULL, All_Properties_single_thread, (void*)temp);
		tHandle[i] = (HANDLE) _beginthread(Cp_single_thread_2, 0, (char *)temp);

	}
	for (i = 0; i < THREADNUM_2; i++)
		WaitForSingleObject(tHandle[i], INFINITE);
		//pthread_join(t[i], NULL);

	double mean_Cp = 0;
	for (i = 0; i < THREADNUM_2; i++)
	{
		mean_Cp += Cpsum[i];		
	}
	cout<<"mean:"<<mean_Cp/N<<endl;
	/*for (i = N-1; i >= 0; i--)
		delete G[i];
	delete []G;*/
	delete []Properties_arg;
	delete []Cpsum;
	delete []tHandle;
	return mean_Cp/N;
}

void Cp_single_thread_2(void * voidarg)
{
	Cp_ARG_2 *arg = (Cp_ARG_2 *)voidarg;
	int id = arg->id;
	//float ** G = arg->G;
	C_type * C = arg->C;
	R_type * R = arg->R;
	V_type * V  = arg->V;
	int N = arg->N;
	double * Cpsum = arg->Cpsum;
	float * Cp = arg->Cp;
	
	int i = 0, j = 0, k = 0, ii = 0, jj = 0;
	double Lsum = 0;
	
	int Neighbour_Num = 0;
	double mean_Cp = 0;
	float Ctemp = 0;
	
	/*for (i = id; i < N; i += THREADNUM)
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
	}*/

	for (i = id; i < N; i += THREADNUM_2)
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
			int Neighbour_Num_j = R[j+1] - R[j];
			for (ii = k+1, jj = 0; ii < Neighbour_Num && jj < Neighbour_Num_j;)
			{
				if(C[R[i]+ii] == C[R[j]+jj])
				{	
					Ctemp += pow( (double) V[R[i]+k]*V[R[i]+ii]*V[R[j]+jj],1.0/3 ); 
					++ii; 
					++jj;
				}
				else if(C[R[i]+ii] < C[R[j]+jj]) ++ii;
				else ++jj;
				/*ii = C[R[i]+k];
				jj = C[R[i]+j];
				if (G[ii][jj]!=0)
				Ctemp += pow( (double)G[ii][jj]*G[i][ii]*G[i][jj],1.0/3 );*/			
			}
		}
		Cp[i] = 2*Ctemp / Neighbour_Num / (Neighbour_Num - 1);
		mean_Cp += Cp[i];
	}
	*Cpsum = (float)mean_Cp;
	return;
}

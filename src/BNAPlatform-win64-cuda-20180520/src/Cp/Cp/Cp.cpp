# include <windows.h>
# include <process.h>
# include <fstream>
# include <algorithm>
# include "data_type.h"

// This struct is used for passing parameters to different threads
struct Cp_ARG
{
	int id;
	//bool ** G;			// The network
	C_type * C;			
	R_type * R;			// R and C represents the network in CSR format
	int N;				// The size of the network
	double *	Cpsum;		
	float * Cp;
};

int THREADNUM = 4;

void Cp_single_thread(void * arg);

double Cp(C_type * C, R_type * R, float * Cp, int N)
{
	int i = 0, j = 0;
	/*bool ** G = new bool * [N];
	for (i = 0; i < N; i++)
	{
		G[i] = new bool [N];
		memset(G[i], 0, sizeof(bool) * (N));
	}
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
			G[i][C[j]] = true;*/

	SYSTEM_INFO siSysInfo;
	GetSystemInfo(&siSysInfo); 
	THREADNUM = (int)siSysInfo.dwNumberOfProcessors; 
	//THREADNUM = sysconf(_SC_NPROCESSORS_CONF);	

	Cp_ARG * Properties_arg = new Cp_ARG [THREADNUM];	
	double * Cpsum = new double[THREADNUM];
	for (i = 0; i < THREADNUM; i++)
	{
		Properties_arg[i].id = i;
		//Properties_arg[i].G = G;
		Properties_arg[i].C = C;
		Properties_arg[i].R = R;
		Properties_arg[i].N = N;
		Properties_arg[i].Cpsum = Cpsum + i;
		Properties_arg[i].Cp = Cp;
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
	//for (i = N-1; i >= 0; i--)
	//	delete G[i];
	//delete []G;
	delete []Properties_arg;
	delete []Cpsum;
	delete []tHandle;
	return mean_Cp/N;
}

void Cp_single_thread(void * voidarg)
{
	Cp_ARG *arg = (Cp_ARG *)voidarg;
	int id = arg->id;
	//bool ** G = arg->G;
	C_type * C = arg->C;
	R_type * R = arg->R;
	int N = arg->N;
	double * Cpsum = arg->Cpsum;
	float * Cp = arg->Cp;
	
	int i = 0, j = 0, k = 0, ii = 0, jj = 0;
	double Lsum = 0;
	
	int Neighbour_Num_i = 0;
	int Neighbour_Num_j = 0;
	double mean_Cp = 0;
	int Ctemp = 0;
	for (i = id; i < N; i += THREADNUM)
	{
		Neighbour_Num_i = R[i+1] - R[i];			
		if (Neighbour_Num_i < 2)
		{
			Cp[i] = 0;
			continue;
		}			
		Ctemp = 0;
		for (k = 0; k < Neighbour_Num_i; k++)
		{			
			j = C[R[i]+k];
			int Neighbour_Num_j = R[j+1] - R[j]; 
			for (ii = k+1, jj = 0; ii < Neighbour_Num_i && jj < Neighbour_Num_j; )
			{
				//linear time complexity, but regular operation
				if(C[R[i]+ii] == C[R[j]+jj]) { ++Ctemp; ++ii; ++jj;}
				else if(C[R[i]+ii] < C[R[j]+jj]) ++ii;
				else ++jj;
				
				//log time complexity, but not faster than the previous one
				/*ii = std::lower_bound(C+R[i]+ii,C+R[i]+Neighbour_Num_i,C[R[j]+jj]) - (C+R[i]);
				if (ii == Neighbour_Num_i) break;
				if(C[R[i]+ii] == C[R[j]+jj]) { ++Ctemp; ++ii; ++jj;}
				
				jj = std::lower_bound(C+R[j]+jj,C+R[j]+Neighbour_Num_j,C[R[i]+ii]) - (C+R[j]);
				if (jj == Neighbour_Num_j) break;
				if(C[R[i]+ii] == C[R[j]+jj]) { ++Ctemp; ++ii; ++jj;}*/
			}
		}
		Cp[i] = (double) Ctemp / Neighbour_Num_i / (Neighbour_Num_i - 1) * 2;
		mean_Cp += Cp[i];
	}
	*Cpsum = (float)mean_Cp;
	return;
}

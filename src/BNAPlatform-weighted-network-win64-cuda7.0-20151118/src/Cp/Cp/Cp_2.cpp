# include <windows.h>
# include <iostream>
# include <process.h>
# include <fstream>
# include <math.h>
using namespace std;
// This struct is used for passing parameters to different threads
struct Cp_ARG_2
{
	int id;
	float ** G;			// The network
	int * C;			
	int * R;			// R and C represents the network in CSR format
	int N;				// The size of the network
	double *	Cpsum;		
	float * Cp;
};

int THREADNUM_2 = 4;

void Cp_single_thread_2(void * arg);

double Cp_2(int * C, int * R,float * V, float * Cp, int N)
{
	int i = 0, j = 0;
	float ** G = new float * [N];
    float V_max = V[0];
    for (i = 1; i < N; i++)
		if(V_max < V[i])
			V_max = V[i];

	for (i = 0; i < N; i++)
	{
		G[i] = new float [N];
		memset(G[i], 0, sizeof(float) * (N));		
	}
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
        {	G[i][C[j]] = V[j]/V_max;
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
		Properties_arg[i].G = G;
		Properties_arg[i].C = C;
		Properties_arg[i].R = R;
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
	for (i = N-1; i >= 0; i--)
		delete G[i];
	delete []G;
	delete []Properties_arg;
	delete []Cpsum;
	delete []tHandle;
	return mean_Cp/N;
}

void Cp_single_thread_2(void * voidarg)
{
	Cp_ARG_2 *arg = (Cp_ARG_2 *)voidarg;
	int id = arg->id;
	float ** G = arg->G;
	int * C = arg->C;
	int * R = arg->R;
	int N = arg->N;
	double * Cpsum = arg->Cpsum;
	float * Cp = arg->Cp;
	
	int i = 0, j = 0, k = 0, ii = 0, jj = 0;
	double Lsum = 0;
	
	int Neighbour_Num = 0;
	double mean_Cp = 0;
	float Ctemp = 0;
	
	for (k = id; k < N; k += THREADNUM_2)
	{
		
		Neighbour_Num = R[k+1] - R[k];
		if (Neighbour_Num < 2)
		{
			Cp[k] = 0;
			continue;
		}			
		Ctemp = 0;
		for (i = 0; i < Neighbour_Num; i++)
			for (j = i+1; j < Neighbour_Num; j++)
			{
				ii = C[R[k]+i];
				jj = C[R[k]+j];
				if (G[ii][jj]!=0)
				Ctemp += pow( (double)G[ii][jj]*G[k][ii]*G[k][jj],1.0/3 );			
			}
		Cp[k] = 2*Ctemp / Neighbour_Num / (Neighbour_Num - 1);//onnela
		mean_Cp += Cp[k];
	}
	*Cpsum = (float)mean_Cp;
	return;
}

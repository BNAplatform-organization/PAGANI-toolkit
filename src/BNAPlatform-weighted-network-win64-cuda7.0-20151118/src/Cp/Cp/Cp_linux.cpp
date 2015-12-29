//# include <windows.h>
//# include <process.h>
# include <pthread.h>
# include <fstream>

// This struct is used for passing parameters to different threads
struct Cp_ARG
{
	int id;
	bool ** G;			// The network
	int * C;			
	int * R;			// R and C represents the network in CSR format
	int N;				// The size of the network
	double *	Cpsum;		
	float * Cp;
};

int THREADNUM = 4;

void Cp_single_thread(void * arg);

double Cp_linux(int * C, int * R, float * Cp, int N)
{
	int i = 0, j = 0;
	bool ** G = new bool * [N];
	for (i = 0; i < N; i++)
	{
		G[i] = new bool [N];
		memset(G[i], 0, sizeof(bool) * (N));
	}
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
			G[i][C[j]] = true;

	//SYSTEM_INFO siSysInfo;
	//GetSystemInfo(&siSysInfo); 
	//THREADNUM = (int)siSysInfo.dwNumberOfProcessors; 
	THREADNUM = sysconf(_SC_NPROCESSORS_CONF);	

	Cp_ARG * Properties_arg = new Cp_ARG [THREADNUM];	
	double * Cpsum = new double[THREADNUM];
	for (i = 0; i < THREADNUM; i++)
	{
		Properties_arg[i].id = i;
		Properties_arg[i].G = G;
		Properties_arg[i].C = C;
		Properties_arg[i].R = R;
		Properties_arg[i].N = N;
		Properties_arg[i].Cpsum = Cpsum + i;
		Properties_arg[i].Cp = Cp;
	}
	pthread_t *t = new pthread_t[THREADNUM];
	//HANDLE *tHandle = new HANDLE[THREADNUM];
	for (i = 0; i < THREADNUM; i++)
	{
		Cp_ARG *temp = Properties_arg + i;
		pthread_create(&t[i], NULL, All_Properties_single_thread, (void*)temp);
		//tHandle[i] = (HANDLE) _beginthread(Cp_single_thread, 0, (char *)temp);

	}
	for (i = 0; i < THREADNUM; i++)
		//WaitForSingleObject(tHandle[i], INFINITE);
		pthread_join(t[i], NULL);

	double mean_Cp = 0;
	for (i = 0; i < THREADNUM; i++)
	{
		mean_Cp += Cpsum[i];		
	}
	for (i = N-1; i >= 0; i--)
		delete G[i];
	delete []G;
	return mean_Cp/N;
}

void Cp_single_thread(void * voidarg)
{
	Cp_ARG *arg = (Cp_ARG *)voidarg;
	int id = arg->id;
	bool ** G = arg->G;
	int * C = arg->C;
	int * R = arg->R;
	int N = arg->N;
	double * Cpsum = arg->Cpsum;
	float * Cp = arg->Cp;
	
	int i = 0, j = 0, k = 0, ii = 0, jj = 0;
	double Lsum = 0;
	
	int Neighbour_Num = 0;
	double mean_Cp = 0;
	int Ctemp = 0;
	for (k = id; k < N; k += THREADNUM)
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
				Ctemp += G[ii][jj];			
			}
		Cp[k] = (double) Ctemp / Neighbour_Num / (Neighbour_Num - 1) * 2;
		mean_Cp += Cp[k];
	}
	*Cpsum = (float)mean_Cp;
	return;
}

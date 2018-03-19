# include <iostream>
# include <fstream>
# include <memory.h>
# include <cmath>
# include <cstring>
# include <ctime>
# include "Timer.h"
# include "dirent.h"
# include "Maslov.h"
using namespace std;

// This is an option for the initial value in Leading_Vector()
//# define RANDOM_V0

// set some parameters

int N;								// The number of voxels
//double * v, * v0, * verr, * sumBG;	// some buffers used in the iteration
ofstream fout;						// logfile

// Forward Declaration

void Partition(long long N_C, int * R, int * C, float * V, int * Result);

//template <class Type> double VectorNorm(Type * x, bool * G, long long N);
//void Maslov_weighted_1(int * R_dst, int * C_dst, float * V_dst, int * R_src, int * C_src, float * V_src, int Rlength, int Clength);



int main(int argc, char * argv[])
{
	ofstream flog("BNA_time_log", ios::app);
	clock_t total_time = clock();
	if (argc > 3 && argc <2) 
	{
		cerr<<"Input format: .\\Modularity.exe dir_for_csr num_of_random_networks \nFor example: .\\Modularity_CPU.exe d:\\data 10"<<endl;
		exit(1);	
	}

	/////////////////////////////////////////////////////////////
	/***************** directory traversal *********************/
	DIR *dp;
	struct dirent *dirp;
	if (NULL == (dp = opendir(argv[1])))
	{
		printf("can't open %s", argv[1]);
		exit (1);
	}
	int FileNumber = 0;
	string filenametmp;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);

		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			FileNumber++;
		}
	}
	std::cout<<FileNumber<<" files to be processed."<<endl;
	closedir(dp);
	string *filename = new string[FileNumber];
	dp = opendir(argv[1]);
	long long i = 0;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			filename[i++] = filenametmp;
		}
	}
	/////////////////////////////////////////////////////////////

	for (i = 0; i < FileNumber; i++)
	{
		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\nModular analysis for "<<a.c_str()<<" ..."<<endl;
		
		/////////////////////////////////////////////////
		/*************** read x.csr file ***************/
		ifstream fin(a.c_str(), ios_base::binary);
		
		
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}
				
		int Rlength = 0, Clength = 0, Vlength = 0;
		
		fin.read((char*)&Rlength, sizeof(int));
		int * R = new int [Rlength];
		fin.read((char*)R, sizeof(int) * Rlength);

		fin.read((char*)&Clength, sizeof(int));
		int * C = new int [Clength];
		fin.read((char*)C, sizeof(int) * Clength);
		
		fin.read((char*)&Vlength, sizeof(int));
		if (Vlength!=Clength)
		{	cout << "csr file is damaged!"<<endl; exit(-1); }
		float * V = new float [Vlength];
		fin.read((char*)V, sizeof(float) * Vlength);
		
		fin.close();
		
		long long N;
		N = Rlength - 1;
		cout<<"Number of voxel = "<<N<<endl;
		int *index = new int [N];
		memset(index, 0 , sizeof(int)*N);
		long long i,j;
		int idx = 0;
		long long N_C;
		//int Nm;
		for (i = 0; i < N; i++)
			if(R[i+1]-R[i] != 0)
				index[i] = idx++;
		N_C = idx;	
		cout<<"Number of connected voxel = "<<N_C<<endl;
		
		
		/************** remove isolated voxels ***********/
		int *R_new = new int [N_C+1];
		int *C_new = new int [Clength];
				
		for (i = 0; i < N; i++)
		{
			if(R[i] == R[i+1]) continue;
			R_new[index[i]] = R[i];
			for (j = R[i]; j < R[i+1] && C[j] < N; j++)
			{
				C_new[j] = index[C[j]];
			}
		}
		R_new [N_C] = R[N];
		//////////////////////////////////////////////////
		

		int *Result = new int [N_C];
		
		string X_modu = a.substr(0, a.find_last_of('.')).append("_modu.nm");
		string X_cp_mas = a.substr(0, a.find_last_of('.')).append("_modu.txt");
		fout.open(X_cp_mas.c_str(), ios::out);	// Open the log file

		Setup(0);
		Start(0);

		// allocate buffers used in the iteration
		/*v = new double [N];
		v0 = new double [N];
		verr = new double [N];
		sumBG = new double [N];
		int * Result = new int [N];	*/
		
		/********************* Partition ******************/
		Partition(N_C, R_new, C_new, V, Result);	
		/**************************************************/

		Stop(0);

		flog<<"Modularity\t"<<a.c_str()<<"CPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;
		cout<<"Modularity\t"<<a.c_str()<<"CPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;

		cout<<"Save partition results as "<<X_modu.c_str()<<endl;
		ofstream fresult;
		fresult.open(X_modu.c_str(), ios::binary|ios::out);
		fresult.write((char*)&N, sizeof(int));

		float * fModu_Result = new float [N];
		idx = 0;
		for (i = 0; i < N; i++)
		{
			if (R[i+1] == R[i])
				fModu_Result[i] = 0;
			else 
				fModu_Result[i] = (float) Result[idx++]+1;
		}

		fresult.write((char*)fModu_Result, sizeof(float) * N);
		fresult.close();
		delete [] fModu_Result;
		

		/////////////////////////////////////////////////////////////
		/*************** Analyze random networks *******************/
				
		int Maslov_num = 0;
		if (argc == 3)
			Maslov_num = atoi(argv[2]);
		
		cout<<"Modular analysis for random networks..."<<endl;
		
		int * R_dst = new int [N_C+1];
		int * C_dst = new int [Clength];
		float * V_dst = new float [Clength];
		Setup(0);
		Start(0);
		for (long long l = 0; l < Maslov_num; l++)
		{
			Maslov_weighted_1(R_dst, C_dst, V_dst, R_new, C_new, V, N_C+1, Clength);
			Partition(N_C, R_dst, C_dst, V_dst, Result);
		}
		Stop(0);
		flog<<"Modularity\tRandom"<<"CPU\t(Maslov+kernel) time\t"<<GetElapsedTime(0)<<"s"<<endl;
		// Clean up
		fout.close();
		
		//delete []sumBG;
		//delete []verr;
		//delete []v0;
		//delete []v;
		delete []R;
		delete []C;	
		delete []R_new;
		delete []C_new;		
		delete []V;
		delete []Result;
		delete []R_dst;
		delete []C_dst;
		delete []V_dst;
	}
	cout<<"==========================================================="<<endl;
	total_time = clock() - total_time;
	flog<<"Modularity\tCPU\ttotal time\t"<<1.0*total_time/1000<<"s"<<endl;
	flog<<endl;
	flog.close();
	delete[]filename;
	return 0;
}




# include <iostream>
# include <iomanip>
# include <fstream>
# include <memory.h>
# include "Timer.h"
# include <stdlib.h>
# include <cstring>
# include "dirent.h"
# include "Maslov.h"
using namespace std;

// This struct is used for passing parameters to different threads
struct Cp_ARG
{
	int id;
	bool ** G;			// The network
	int * C;			
	int * R;			// R and C represents the network in CSR format
	int N;				// The size of the network
	double * Cpsum;		
	float * Cp;

};


// Forward Declaration
double Cp_1(int * C, int * R, float * V ,float * Cp, int N);
double Cp_2(int * C, int * R, float * V ,float * Cp, int N);
//void Maslov_weighted_1(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);

int main(int argc, char * argv[])
{
	clock_t total_time = clock();
	ofstream flog("BNA_time_log", ios::app);
	if (argc < 2) 
	{
		cerr<<"Input format: .\\Cp.exe dir_for_csr num_of_random_networks \nFor example: .\\Cp.exe d:\\data 10(option)"<<endl;
		exit(1);	
	}

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
		if(filenametmp.length()>4 && 
			filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0
			&& filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			FileNumber++;
		}
	}
	cout<<FileNumber<<" files to be processed."<<endl;

	closedir(dp);
	string *filename = new string[FileNumber];
	dp = opendir(argv[1]);
	int i = 0;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && 
			filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 
			&& filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			filename[i++] = filenametmp;
		}
	}

	for (int i = 0; i < FileNumber; i++)
	{
		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\ncalculating Cp for "<<a.c_str()<<" ..."<<endl;
		ifstream fin(a.c_str(), ios_base::binary);
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}

		// Read x.csr
		int Rlength = 0, Clength = 0, Clength1=0;
		fin.read((char*)&Rlength, sizeof(int));
		int * R = new int [Rlength];
		fin.read((char*)R, sizeof(int) * Rlength);
		fin.read((char*)&Clength, sizeof(int));
		int * C = new int [Clength];
		fin.read((char*)C, sizeof(int) * Clength);
		fin.read((char*)&Clength1, sizeof(int));
		float * V = new float [Clength];
		fin.read((char*)V, sizeof(float) * Clength);
		fin.close();
		int N = Rlength - 1;

		// Generate the bool network using CSR
		int l = 0;
		Setup(0);
		Start(0);

		double mean_Cp;
		float * Cp_result = new float [N];
		mean_Cp = Cp_2(C, R, V, Cp_result, N);	
		Stop(0);
		cout<<"average Cp:\t"<<setprecision(6)<<mean_Cp<<endl;
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;


		flog<<"Cp\t"<<a.c_str()<<"CPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;


		// Parse file name
		string X_cp = a.substr(0, a.find_last_of('.') + 1).append("cp");
		string X_cp_mas = a.substr(0, a.find_last_of('.')).append("_cp.txt");
		cout<<"Save Clustering Coefficient for each node as "<<X_cp.c_str()<<endl;
		ofstream fout;
		fout.open(X_cp.c_str(), ios::binary|ios::out);
		fout.write((char*)&N, sizeof(int));
		fout.write((char*)Cp_result, sizeof(float) * N);
		fout.close();

		// Analysis for random networks
		
		fout.open(X_cp_mas.c_str(), ios::out);
		fout<<setprecision(6)<<mean_Cp<<endl;
		
		int Maslov_num = atoi(argv[2]);
		cout<<"Calculating Cp for random networks..."<<endl;
		Reset(0);
		Start(0);
		double * Maslov_mean_Cp = new double [Maslov_num + 1];
		int * R_dst = new int [Rlength];
		int * C_dst = new int [Clength];
		float * V_dst = new float [Clength];
		for (l = 0; l < Maslov_num; l++)
		{
			Maslov_weighted_1(R_dst, C_dst,V_dst, R, C, V, Rlength, Clength);
			Maslov_mean_Cp[l] = Cp_1(C_dst, R_dst, V_dst, Cp_result, N);
			fout<<Maslov_mean_Cp[l]<<endl;
		}
		fout.close();
		Stop(0);
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;
		
		flog<<"Cp\tRandom"<<"CPU\t(Maslov+kernel) time\t"<<GetElapsedTime(0)<<"s"<<endl;

		
		cout<<"Save average Cp for the brain network and random networks as "<<X_cp_mas.c_str()<<endl;
		
		
		delete []R_dst;
		delete []C_dst;
		delete []V_dst;
		delete []Cp_result;
		delete []R;
		delete []C;
		delete []V;
	}
	delete []filename;
	total_time = clock() - total_time;
	cout<<"total elapsed time: "<<1.0*total_time/1000<<" s."<<endl;
	cout<<"==========================================================="<<endl;
	
	flog<<"Cp\tCPU\ttotal time\t"<<1.0*total_time/1000<<"s"<<endl;
	flog<<endl;
	flog.close();
	return 0;
}

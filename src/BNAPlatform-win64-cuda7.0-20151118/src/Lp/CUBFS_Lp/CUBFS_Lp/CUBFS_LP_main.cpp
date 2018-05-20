#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <iomanip>
#include "dirent.h"
# include "Timer.h"

using namespace std;
void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);
double CUBFS_Lp(int *row, int *col, int numVertices, int numEdges);

int main(int argc, char** argv) 
{
	// input 1: file name
	// input 2 & 3£º for tuning performance
	//int grid = atoi(argv[2]);
	//int thread = atoi(argv[3]);
	ofstream flog("BNA_time_log", ios::app);
	clock_t total_time = clock();
	if (argc < 3) 
	{
		cerr<<"Input format: .\\CUBFS_Lp.exe dir_for_csr num_of_random_networks \nFor example: .\\CUBFS_Lp.exe d:\\data 10"<<endl;
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
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
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
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			filename[i++] = filenametmp;
		}
	}
	int Nrandom = atoi(argv[2]);
	double *Lp = new double[1+Nrandom];

	for (int i = 0; i < FileNumber; i++)
	{
		
		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\ncalculating Lp for "<<a.c_str()<<" ..."<<endl;

		int numVertices;	// Network vertices number
		int numEdges;	//Network egdes number
		int *row = NULL;	// Cormat CSR row
		int *col = NULL;	// Cormat CSR col

		ifstream fCSR(a.c_str(), ios::binary);
		if(!fCSR.good())
		{
			cerr<<"cannot open file!"<<endl;
			exit(1);
		}
		fCSR.read((char *)&numVertices, sizeof(int));
		numVertices = numVertices - 1;	
		row = new int[numVertices+1];	// Cormat CSR row
		fCSR.read((char*)row, sizeof(int)*(numVertices + 1));
		fCSR.read((char *)&numEdges, sizeof(int));
		col = new int[numEdges];	// Cormat CSR col
		fCSR.read((char*)col, sizeof(int)*numEdges);
		fCSR.close();
		Setup(0);
		Start(0);
		Lp[0] = CUBFS_Lp(row, col, numVertices, numEdges);
		Stop(0);
		cout<<"average Lp:\t"<<setprecision(6)<<Lp[0]<<endl;
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;

		flog<<"BFS_Lp\t"<<a.c_str()<<"GPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;

		cout<<"Calculating Lp for random networks..."<<endl;
		int Rlength = numVertices + 1;
		int Clength = numEdges; 
		int * R_dst = new int [Rlength];
		int * C_dst = new int [Clength];
		Reset(0);
		Start(0);
		for (int l = 0; l < Nrandom; l++)
		{
			Maslov(R_dst, C_dst, row, col, Rlength, Clength);
			Lp[l+1] = CUBFS_Lp(R_dst, C_dst, Rlength - 1, Clength);
		}
		Stop(0);
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;

		flog<<"BFS_Lp\tRandom"<<"GPU\t(Maslov+kernel) time\t"<<GetElapsedTime(0)<<"s"<<endl;

		string Lpoutfile = a.substr(0, a.find_last_of('.')).append("_Lp.txt");
		ofstream fLp(Lpoutfile.c_str());
		for (int i = 0; i < Nrandom + 1; i ++)
			fLp<<setprecision(6)<<Lp[i]<<endl;
		fLp.close();
		cout<<"Save Lp for the brain network and random networks as "<<Lpoutfile.c_str()<<endl;
		delete []R_dst;
		delete []C_dst;
		delete []row;
		delete []col;
	}
	delete []Lp;
	delete []filename;
	total_time = clock() - total_time;
	cout<<"total elapsed time: "<<1.0*total_time/1000<<" s."<<endl;
	cout<<"==========================================================="<<endl;

	flog<<"BFS_Lp\tGPU\ttotal time\t"<<1.0*total_time/1000<<"s"<<endl;
	flog<<endl;
	flog.close();
	return 0;
}


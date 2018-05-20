# include <iostream>
# include <iomanip>
# include <fstream>
# include <memory.h>
# include <stdlib.h>
# include <cstring>
# include <ctime>
# include "dirent.h"
# include "data_type.h"
using namespace std;

typedef float real__t;

int main(int argc, char * argv[])
{
	ofstream flog("BNA_time_log", ios::app);
	clock_t total_time = clock();
	if (argc < 2) 
	{
		cerr<<"Input format: .\\Deg.exe dir_for_csr \nFor example: .\\Deg.exe d:\\data "<<endl;
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

	for (int i = 0; i < FileNumber; i++)
	{
		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\ncalculating Degree for "<<a.c_str()<<" ..."<<endl;
		ifstream fin(a.c_str(), ios_base::binary);
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}

		// Read x.csr
		unsigned int Rlength = 0;
		R_type Clength = 0, Vlength = 0;
		fin.read((char*)&Rlength, sizeof(int));
		R_type * R = new R_type[Rlength];
		fin.read((char*)R, sizeof(R_type) * Rlength);
		fin.read((char*)&Clength, sizeof(R_type));
		C_type * C = new C_type[Clength];
		fin.read((char*)C, sizeof(C_type) * Clength);
		fin.read((char*)&Vlength, sizeof(R_type));
		V_type * V = new V_type[Vlength];
		fin.read((char*)V, sizeof(V_type) * Clength);
		fin.close();
		int N = Rlength - 1;
				
		clock_t time = clock();
		float *Degree = new float [N];
		memset(Degree,0,sizeof(float)*N);
		for (C_type i = 0; i < N; i++) {
			double deg = 0;
			for (size_t j = R[i]; j < R[i + 1]; j++)
				deg += V[j];
			Degree[i] = deg;
		}
		
		time = clock() - time;
		cout<<"Elapsed time: "<<time<<" ms. "<<endl;
		flog<<"Degree\t"<<a.c_str()<<"CPU\tkernel time\t"<<1.0*time/1000<<"s"<<endl;
		// Parse file name
		string X_deg = a.substr(0, a.find_last_of('.')).append("_deg.nm");
		cout<<"Save degree data for each node as "<<X_deg.c_str()<<endl;
		ofstream fout;
		fout.open(X_deg.c_str(), ios::binary|ios::out);
		fout.write((char*)&N, sizeof(int));
		fout.write((char*)Degree, sizeof(float) * N);
		fout.close();

		delete []Degree;
		delete []R;
		delete []C;
		delete []V;
	}
	total_time = clock() - total_time;
	cout<<"total elapsed time: "<<1.0*total_time/1000<<" s."<<endl;
	cout<<"==========================================================="<<endl;
	flog<<"Degree\tCPU\ttotal time\t"<<1.0*total_time/1000<<"s"<<endl;
	flog<<endl;
	flog.close();
	delete []filename;
	return 0;
}

#include "Lp_CPU.h"
#include <iomanip>
#include "dirent.h"
void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);
int main( int argc, char *argv[])
{
	clock_t total_time = clock();
	if (argc < 3) 
	{
		cerr<<"Input format: .\\Lp.exe dir_for_csr num_of_random_networks \nFor example: .\\Lp_CPU.exe d:\\data 10"<<endl;
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
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0)
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
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0)
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
		ifstream fin(a.c_str(), ios_base::binary);
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}

		// Read x.csr
		int Rlength = 0, Clength = 0;
		fin.read((char*)&Rlength, sizeof(int));
		int * R = new int [Rlength];
		fin.read((char*)R, sizeof(int) * Rlength);
		fin.read((char*)&Clength, sizeof(int));
		int * C = new int [Clength];
		fin.read((char*)C, sizeof(int) * Clength);
		fin.close();
		int N = Rlength - 1;
	
		Lp[0] = Lp_CPU(R, C, N, Clength);
		cout<<"Lp:\t"<<setprecision(6)<<Lp[0]<<endl;

		//break;

		cout<<"Calculating Lp for random networks..."<<endl;
		Reset(0);
		Start(0);
		int * R_dst = new int [Rlength];
		int * C_dst = new int [Clength];
		for (int l = 0; l < Nrandom; l++)
		{
			Maslov(R_dst, C_dst, R, C, Rlength, Clength);
			Lp[l+1] = Lp_CPU(R, C, N, Clength);
		}

		string inputfile = string(argv[1]);
		string Lpoutfile = inputfile.append("\\").append(filename[i]).append("_Lp.txt"); 
		ofstream fLp(Lpoutfile.c_str());
		for (int i = 0; i < Nrandom + 1; i ++)
			fLp<<setprecision(6)<<Lp[i]<<endl;
		fLp.close();
		cout<<"Save average Lp for the brain network and random networks as "<<Lpoutfile.c_str()<<endl;

		delete []R_dst;
		delete []C_dst;
		delete []R;
		delete []C;
	}
	delete []Lp;
	delete []filename;
	total_time = clock() - total_time;
	cout<<"total elapsed time: "<<1.0*total_time/1000<<" s."<<endl;
	cout<<"==========================================================="<<endl;
	return 0;
}
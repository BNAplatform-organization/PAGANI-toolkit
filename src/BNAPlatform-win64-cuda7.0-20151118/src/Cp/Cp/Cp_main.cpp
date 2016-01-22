# include <iostream>
# include <iomanip>
# include <fstream>
# include <memory.h>
# include "Timer.h"
# include <stdlib.h>
# include <cstring>
# include "dirent.h"
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
double Cp(int * C, int * R, float * Cp, int N);
void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);

int main(int argc, char * argv[])
{
	clock_t total_time = clock();
	ofstream flog("BNA_time_log", ios::app);
	if (argc < 4) //amend
	{
		cerr<<"Input format: .\\Cp.exe dir_for_csr num_of_random_networks parameter_type \nFor example: .\\Cp.exe d:\\data 10 n"<<endl; //amend
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

	bool global=false,nodal=false,gamma=false;
		//string argv3(argv[3]);
		char* S[]={"g","n","k","gn","gk","nk","gnk"};
       for(int i=0;i<7;i++)
{
  if(strcmp(argv[3],S[i])==0)
  {
     switch(i)
     {
        case 0: global=true;break;
		case 1: nodal=true;break;
		case 2: gamma=true;break;
		case 3:nodal=true;global=true;break;
		case 4:global=true;gamma=true;break;
		case 5:nodal=true;gamma=true;break;
		case 6:global=true;nodal=true;gamma=true;break;
		default:global=true;nodal=true;gamma=true;break;
        
     }
	 break;
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
		int Rlength = 0, Clength = 0;
		fin.read((char*)&Rlength, sizeof(int));
		int * R = new int [Rlength];
		fin.read((char*)R, sizeof(int) * Rlength);
		fin.read((char*)&Clength, sizeof(int));
		int * C = new int [Clength];
		fin.read((char*)C, sizeof(int) * Clength);
		fin.close();
		int N = Rlength - 1;

		// Generate the bool network using CSR
		int l = 0;
		Setup(0);
		Start(0);

		double mean_Cp;
		float * Cp_result = new float [N];
		mean_Cp = Cp(C, R, Cp_result, N);	
		Stop(0);
		cout<<"average Cp:\t"<<setprecision(6)<<mean_Cp<<endl;
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;


		flog<<"Cp\t"<<a.c_str()<<"CPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;

		
		// Parse file name£»output when parameter type is 'n'
/*		bool both=false;
		if(*argv[3] == 'b')
			both=true;
		if(*argv[3] == 'n'||both)  {*/
		if(nodal==true)
		{
		string X_cp = a.substr(0, a.find_last_of('.') ).append("_cp.nm");
		cout<<"Save Clustering Coefficient for each node as "<<X_cp.c_str()<<endl;
		ofstream fout;
		fout.open(X_cp.c_str(), ios::binary|ios::out);
		fout.write((char*)&N, sizeof(int));
		fout.write((char*)Cp_result, sizeof(float) * N);
		fout.close();
		}
	//	}
		// Analysis for random networks
		int * R_dst = new int [Rlength];
		int * C_dst = new int [Clength];
	//	if(*argv[3] == 'g'||both)  {
		int Maslov_num = atoi(argv[2]);
		ofstream fout;
		
		cout<<"Calculating Cp for random networks..."<<endl;
		Reset(0);
		Start(0);
		double * Maslov_mean_Cp = new double [Maslov_num + 1];
		Maslov_mean_Cp[0]=mean_Cp;
		double  small_world_Cp=0;
		for (l = 0; l < Maslov_num; l++)
		{
			Maslov(R_dst, C_dst, R, C, Rlength, Clength);
			Maslov_mean_Cp[l+1] = Cp(C_dst, R_dst, Cp_result, N);
		//	fout<<Maslov_mean_Cp[l+1]<<endl;
			small_world_Cp+=Maslov_mean_Cp[l+1];
		}
		
		Stop(0);
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;
		small_world_Cp/=Maslov_num;
		small_world_Cp=mean_Cp/small_world_Cp;
		if(global==true)
		{
		string X_cp_mas = a.substr(0, a.find_last_of('.')).append("_cp.txt");
		cout<<"Save average Cp for the brain network and random networks as "<<X_cp_mas.c_str()<<endl;
		fout.open(X_cp_mas.c_str(), ios::out);
		for (int i = 0; i < Maslov_num + 1; i ++)
		fout<<setprecision(6)<<Maslov_mean_Cp[i]<<endl;
		fout.close();
		}
		if(gamma==true)
		{
		string X_cp_mas = a.substr(0, a.find_last_of('.')).append("_gamma.txt");
		cout<<"Save corresponding gamma for the brain network and random networks as "<<X_cp_mas.c_str()<<endl;
		fout.open(X_cp_mas.c_str(), ios::binary|ios::out);
		fout<<setprecision(6)<<small_world_Cp<<endl;
		fout.close();
		}
#ifdef Debug
		X_cp_mas = a.substr(0, a.find_last_of('.')).append("_small_world_cp.txt");
		fout.open(X_cp_mas.c_str(),ios::out);
		fout<<setprecision(6)<<small_world_Cp<<endl;
		fout.close();
#endif
		flog<<"Cp\tRandom"<<"CPU\t(Maslov+kernel) time\t"<<GetElapsedTime(0)<<"s"<<endl;

		
	//	}
		delete []R_dst;
		delete []C_dst;
		delete []Cp_result;
		delete []R;
		delete []C;
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

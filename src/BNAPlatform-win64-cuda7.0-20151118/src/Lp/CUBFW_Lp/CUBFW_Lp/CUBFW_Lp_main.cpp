 # include <iostream>
# include <fstream>
# include <ctime>
#include <iomanip>
#include "dirent.h"
# include "Timer.h"
# include "APSP_BFS.h"
using namespace std;
//amend:auto scheduling. method:1.obtain the sparsity of each files and store in filesparsity[];2.judge befpre execute BFS or BFW
#define sparsity  1
float *filesparsity;
void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);
double CUBFW_Lp(int *row, int *col, int numVertices, int numEdges);

float *Li_result;



int main(int argc, char *argv[]){
	ofstream flog("BNA_time_log", ios::app);
	clock_t total_time = clock();
	if (argc < 4) 
	{
		cerr<<"Input format: .\\CUBFW_Lp.exe dir_for_csr num_of_random_networks parameter_type \nFor example: .\\CUBFW_Lp.exe d:\\data 10 g"<<endl;
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
	filesparsity = new float[FileNumber];
	dp = opendir(argv[1]);
	int i = 0;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".csr") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			filename[i] = filenametmp;
			//amend_flag_1
			filesparsity[i]=atof(filenametmp.substr((filenametmp.find("spa")+3),5).c_str());
		//	cout<<filesparsity[i]<<"%"<<endl;
			i++;
		}
	}
	int Nrandom = atoi(argv[2]);
	double *Lp = new double[1+Nrandom];

		bool global=false,nodal=false,lambda=false;
		//string argv3(argv[3]);
		char* S[]={"g","n","l","gn","gl","nl","gnl"};
       for(int i=0;i<7;i++)
{
  if(strcmp(argv[3],S[i])==0)
  {
     switch(i)
     {
        case 0: global=true;break;
		case 1: nodal=true;break;
		case 2: lambda=true;break;
		case 3:nodal=true;global=true;break;
		case 4:global=true;lambda=true;break;
		case 5:nodal=true;lambda=true;break;
		case 6:global=true;nodal=true;lambda=true;break;
		default:global=true;nodal=true;lambda=true;break;
        
     }
	 break;
  } 
   
}

	for (int i = 0; i < FileNumber; i++)
	{
		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\ncalculating Lp for "<<a.c_str()<<" ..."<<endl;

		ifstream fCSR(a.c_str(), ios_base::binary);
		if (!fCSR.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}

		int Nrandom = 0;
		Nrandom = atoi(argv[2]);

		int 	numVertices;	// Network vertices number
		fCSR.read((char *)&numVertices, sizeof(int));
		numVertices = numVertices - 1;	
		int 	*row = new int[numVertices+1];	// Cormat CSR row
		fCSR.read((char*)row, sizeof(int)*(numVertices + 1));
		int 	numEdges;	//Network egdes number
		fCSR.read((char *)&numEdges, sizeof(int));
		int 	*col = new int[numEdges];	// Cormat CSR col
		fCSR.read((char*)col, sizeof(int)*numEdges);
		fCSR.close();

		////generate a symmetric sparse graph
		//int numVertices = 15000;
		//int numEdges = 0;
		//float	*matrix = new float[squareVertices];
		//if (matrix == NULL)
		//{
		//cout<<"Allocation failure"<<endl;
		//return 0;
		//}
		//memset((void*)matrix, 0, sizeof(float) * squareVertices);

		//double sparsity = 0.001;
		//float threshold  = 0.95;
		//
		//for (int i = 0; i < numVertices; i++)
		//{
		//	for (int j = i + 1; j < numVertices; j++)
		//	{
		//		double r = 1.0 * rand()/RAND_MAX;
		//		//float r = matrix[i * numVertices + j];
		//		if (r < sparsity)
		//		//if (r > threshold)
		//		{
		//			matrix[i * numVertices + j] = 1;
		//			matrix[j * numVertices + i] = 1;
		//			numEdges += 2;
		//		}
		//	}
		//}
		//cout<<"sparseness: "<<100.0*numEdges/numVertices/numVertices<<" %"<<endl;

		//int *row = new int[numVertices+1];
		//int *col = new int[numEdges];
		//int count = 0;
		//int colIndex = 0;
		//row[0] = 0;
		//for(int i = 0; i < numVertices; i++)
		//{
		//	for(int j = 0; j < numVertices; j++)
		//	{
		//		if(matrix[i * numVertices + j] ==1)
		//		{
		//			count ++;
		//			col[colIndex++] = j;
		//		}
		//	}
		//	row[i+1] = count;
		//}

		//delete []matrix;
		Li_result = new float[numVertices];
		//amend_flag_2
		if(filesparsity[i]>=1)  {
		cout<<"The file "<<i<<" should use BFW algorithm for sparsity: "<<filesparsity[i]<<'%'<<endl;
		Setup(0);
		Start(0);
		Lp[0] = CUBFW_Lp(row, col, numVertices, numEdges);
		//float *APSP_output = new float[(long long)numVertices * numVertices];
		//Lp[0] =APSP_BFS(APSP_output, row, col, numVertices);
		//delete []APSP_output;


		Stop(0);
		}
		else  {
		float *APSP_output=new float[(long long)numVertices*numVertices];
		cout<<"The file "<<i<<" should use BFS algorithm for sparsity: "<<filesparsity[i]<<'%'<<endl;
		Setup(0);
		Start(0);
		Lp[0]=APSP_BFS(APSP_output,row,col,numVertices);
		//Lp[0] = CUBFS_Lp(row, col, numVertices, numEdges);
		Stop(0);
		
		delete []APSP_output;
		}

		

		string X_eff = a.substr(0, a.find_last_of('.') ).append("_eff.nm");
		ofstream fout;
		int Rlength = numVertices + 1;
		int Clength = numEdges; 
		int * R_dst = new int [Rlength];
		int * C_dst = new int [Clength];
	
	
	/*	bool both=false;
		if(*argv[3]=='b')
			both=true;
		if(*argv[3]=='n'||both)  {*/
		if(nodal==true)
		{
		cout<<"Save nodal efficiency for each node as "<<X_eff.c_str()<<endl;
		
		fout.open(X_eff.c_str(), ios::binary|ios::out);
		fout.write((char*)&numVertices, sizeof(int));
		fout.write((char*)Li_result, sizeof(float) * numVertices);
		fout.close();
		}
	//	}
		//delete [] Li_result;
//		if(*argv[3]=='g'||both)  {
		cout<<"average Lp:\t"<<setprecision(6)<<Lp[0]<<endl;
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;
		flog<<"BFW_Lp\t"<<a.c_str()<<"GPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;
	
		cout<<"Calculating Lp for random networks..."<<endl;

		
		Setup(0);
		Start(0);
		double small_word_Lp=0;
		for (int l = 0; l < Nrandom; l++)
		{
			Maslov(R_dst, C_dst, row, col, Rlength, Clength);
			Lp[l+1] = CUBFW_Lp(R_dst, C_dst, Rlength - 1, Clength);
			//float *APSP_output = new float[(long long)numVertices * numVertices];
			//Lp[l+1] =APSP_BFS(APSP_output, R_dst, C_dst, numVertices);
			//delete []APSP_output;
		   small_word_Lp+=Lp[l+1];
		}
		   small_word_Lp/=Nrandom;
		   small_word_Lp=Lp[0]/small_word_Lp;
		Stop(0);
		cout<<"Elapsed time: "<<GetElapsedTime(0)<<" s."<<endl;

		flog<<"BFW_Lp\tRandom"<<"GPU\t(Maslov+kernel) time\t"<<GetElapsedTime(0)<<"s"<<endl;

		if(global==true)
		{
		string Lpoutfile = a.substr(0, a.find_last_of('.')).append("_Lp.txt");
		ofstream fLp(Lpoutfile.c_str());
		for (int i = 0; i < Nrandom + 1; i ++)
			fLp<<setprecision(6)<<Lp[i]<<endl;
		fLp.close();
		cout<<"Save Lp for the brain network and random networks as "<<Lpoutfile.c_str()<<endl;
		}
		if(lambda==true)
		{
		  string  Lpoutfile = a.substr(0, a.find_last_of('.')).append("_lambda.txt");
		  cout<<"Save corresponding lambda for the brain network and random networks as "<<Lpoutfile.c_str()<<endl;
		 ofstream fLp;
		  fLp.open(Lpoutfile.c_str(),ios::binary|ios::out);
		 fLp<<setprecision(6)<<small_word_Lp<<endl;
		 fLp.close();
		 }
#ifdef Debug
		 Lpoutfile = a.substr(0, a.find_last_of('.')).append("_small_world_Lp.txt");
		 fLp.open(Lpoutfile.c_str(),ios::out);
		 fLp<<setprecision(6)<<small_word_Lp<<endl;
		 fLp.close();
#endif
		
		//}
		delete []row;
		delete []col;
		delete []R_dst;
		delete []C_dst;
		delete []Li_result;
	
	}
	delete []filesparsity;
	delete []Lp;
	delete []filename;
	total_time = clock() - total_time;
	cout<<"Total elapsed time: "<<1.0*total_time/1000<<" s."<<endl
		<<"==========================================================="<<endl;
	flog<<"BFW_Lp\tGPU\ttotal time\t"<<1.0*total_time/1000<<"s"<<endl;
	flog<<endl;
	flog.close();
	return 0;
}


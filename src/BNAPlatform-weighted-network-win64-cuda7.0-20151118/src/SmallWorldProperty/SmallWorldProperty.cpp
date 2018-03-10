// SmallWorldProperty.cpp : 定义控制台应用程序的入口点。
#include "stdafx.h"
# include <iostream>
# include <iomanip>
# include <fstream>
# include <memory.h>
# include "Timer.h"
# include <stdlib.h>
# include <cstring>
# include "dirent.h"
using namespace std;
//#define Debug

int main(int argc, char * argv[])
{
	clock_t total_time = clock();
	ofstream flog("BNA_time_log", ios::app);
	/*if (argc < 4) //amend
	{
		cerr<<"Input format: .\\Cp.exe dir_for_csr num_of_random_networks parameter_type \nFor example: .\\Cp.exe d:\\data 10 n"<<endl; //amend
		exit(1);	
	}*/
	DIR *dp;
	struct dirent *dirp;
	//int i = 0;
	if (NULL == (dp = opendir(argv[1])))
	{
		printf("can't open %s", argv[1]);
		exit (1);
	}
	int swlpFileNumber = 0;
	int swcpFileNumber = 0;
	string filenametmp;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		
		if (filenametmp.find_last_of('.') == -1)
			continue;       //4+7
		if(filenametmp.length()>11 &&filenametmp.substr(filenametmp.find_last_of('.')-7,7).compare("_lambda") == 0&& filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".txt") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			swlpFileNumber++;

		}                    //4+6
		if(filenametmp.length()>10 &&filenametmp.substr(filenametmp.find_last_of('.')-6,6).compare("_gamma") == 0&& filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".txt") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			swcpFileNumber++;

		}
	}
	cout<<swlpFileNumber<<" files within lambda value to be processed."<<endl;
	cout<<swcpFileNumber<<" files within gamma value to be processed."<<endl;

	closedir(dp);
	if(swlpFileNumber==0)
	{
		cerr<<"Lack lambda files. Please run lambda algorithm before executing small world algorithm."<<endl;
	}
	if(swcpFileNumber==0)
	{
	   cerr<<"Lack gamma files. Please run gamma algorithm before executing small world algorithm."<<endl;
	}
	cout<<"Matching...."<<endl;
	string *swlpfilename = new string[swlpFileNumber];
	string *swcpfilename = new string[swcpFileNumber];
	dp = opendir(argv[1]);
	int i = 0;
	int j = 0;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>11 &&filenametmp.substr(filenametmp.find_last_of('.')-7,7).compare("_lambda") == 0&& filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".txt") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			swlpfilename[i++] = filenametmp;
		}
		if(filenametmp.length()>10 &&filenametmp.substr(filenametmp.find_last_of('.')-6,6).compare("_gamma") == 0&& filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".txt") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			swcpfilename[j++] = filenametmp;
		}
	}
#ifdef  Debug
	for(j=0;j<swcpFileNumber;j++)
		cout<<swcpfilename[j].c_str()<<endl;
#endif
	//the data structure of file matching
	char *swcp_match=new char[swcpFileNumber];
	memset(swcp_match, false, sizeof(char)*swcpFileNumber);
	
	char * sucessfulmatch[2];
	sucessfulmatch[0]=new char[swlpFileNumber];
	sucessfulmatch[1]=new char[swlpFileNumber];
	memset(sucessfulmatch[0], 0, sizeof(char)*swlpFileNumber);
	memset(sucessfulmatch[1], 0, sizeof(char)*swlpFileNumber);
	int  sucessfulmatch_count=0;

	char *swlp_single=new char[swlpFileNumber];
	char *swcp_single=new char[swcpFileNumber];
	memset(swlp_single, 0, sizeof(char)*swlpFileNumber);
	memset(swcp_single, 0, sizeof(char)*swcpFileNumber);
	int swlpsingle_count=0,swcpsingle_count=0;
	
	string feature[2];
	//the algorithm of file matching
	for (int i = 0; i <swlpFileNumber; i++)
	{
		//input:swlpfilename and swcpfilename;output:sucessfulmatch, swlp_single,swcp_single
		//1 feature extracting
		feature[0]=swlpfilename[i].substr(1,4);
		feature[1]=swlpfilename[i].substr((swlpfilename[i].find("spa")+3),5);
		//2 matching
		bool successfulmatchflag=false;
		for(j = 0; j <swcpFileNumber; j++)
		{  
			if(swcp_match[j]==true)
				continue;
			else if(feature[0].compare(swcpfilename[j].substr(1,4))==0&&feature[1].compare(swcpfilename[j].substr((swcpfilename[j].find("spa")+3),5))==0)  
			{ 
			  successfulmatchflag=true;
			  sucessfulmatch[0][sucessfulmatch_count]=i;
			  sucessfulmatch[1][sucessfulmatch_count++]=j;
			  swcp_match[j]=true;
			  break;
			}
			else
			  continue;
		}
		if(successfulmatchflag==false)
			  	swlp_single[swlpsingle_count++]=i;
	
	}
	//3 output swcp_single
	for(int j = 0; j <swcpFileNumber; j++)
	{
		if (swcp_match[j]==false)
		{
          swcp_single[swcpsingle_count++]=j;
 
		}
	}
    //output files that single
	if(swlpsingle_count!=0)//||swcpsingle_count!=0)
	{
		cout<<swlpsingle_count<<" files within lambda value can't be used to calculate the small world property because of default corresponding clustering coefficient ratio(gamma value) files"<<endl;
		cout<<" These files are as follows:"<<endl;
		for (int i = 0; i <swlpsingle_count; i++)
		{
		  cout<<"No."<<i+1<<"  "<<"name: "<<swlpfilename[swlp_single[i]].c_str()<<"  "<<endl;
		}
	}
	 //output files that single
	if(swcpsingle_count!=0)
	{
		cout<<swcpsingle_count<<" files within gamma value can't be used to calculate the small world property because of default corresponding  characteristic length path(lambda value) files"<<endl;
		cout<<" These files are as follows:"<<endl;
		for (int i = 0; i <swcpsingle_count; i++)
		{
		 cout<<"No."<<i+1<<"  "<<"name:  "<<swcpfilename[swcp_single[i]].c_str()<<"  "<<endl;
		}
	}
	//calculate the gama
	double *sigma=new double[sucessfulmatch_count];
	double lambda=0;
	double gamma=0;
	for (int i = 0; i <sucessfulmatch_count; i++)
	{
		string temp= string(argv[1]).append("\\").append("test-retest.txt");
		//temp="E:\BNAlauncher\BNU_data_S2_csr\test-retest.txt";
		string a = string(argv[1]).append("\\").append(swlpfilename[sucessfulmatch[0][i]]);
		string b = string(argv[1]).append("\\").append(swcpfilename[sucessfulmatch[1][i]]);
		cout<<"\ncalculating small world property for "<<a.c_str()<<" and "<<b.c_str()<<" ..."<<endl;
		
		ifstream finlp(a.c_str(),ios_base::binary);
		if (!finlp.good())
		{	cout<<"Can't open\t"<<temp.c_str()<<endl;	return 0;}
		finlp>>lambda;
	    finlp.close();
		ifstream fincp(b.c_str(), ios_base::binary);
		if (!fincp.good())
		{	cout<<"Can't open\t"<<b.c_str()<<endl;	return 0;}
		fincp>>gamma;
		fincp.close();
#ifdef Debug		
		cout<<gamma<<endl;
		cout<<lambda<<endl;
#endif		
		sigma[i]=gamma/lambda;
		cout<<"Small world property equals "<<setprecision(6)<<sigma[i]<<endl;
	}
	printf("\n");
	cout<<"Conclude each result as the table view:"<<endl;
	cout<<" File   sparsity   small world property "<<endl;
	string *filename=new string[sucessfulmatch_count];
	for (int i = 0; i <sucessfulmatch_count; i++)
	{ 
	  //get the original file name
	  filename[i]=swlpfilename[sucessfulmatch[0][i]].substr(0,swlpfilename[sucessfulmatch[0][i]].find_last_of('_'));
	  //table output
	  feature[0]=swlpfilename[sucessfulmatch[0][i]].substr(1,4);
	  feature[1]=swlpfilename[sucessfulmatch[0][i]].substr((swlpfilename[sucessfulmatch[0][i]].find("spa")+3),5);
	  cout<<" "<<feature[0].c_str()<<"    "<<feature[1].c_str()<<"%          "<<setprecision(6)<<sigma[i]<<endl;
	}
	cout<<"Already save each result to the corresponding file with .txt format..."<<endl;
	for (int i = 0; i <sucessfulmatch_count; i++)
	{
	ofstream fout;
#ifdef Debug
	cout<<filename[i].c_str()<<endl;
#endif
	string X = string(argv[1]).append("\\").append(filename[i]);
	string X_small_world_property =  X.append("_smallWorld.txt");
	fout.open(X_small_world_property.c_str(), ios::binary|ios::out);
	fout<<setprecision(6)<<sigma[i]<<endl;
	fout.close();
	}
	
	delete[] swlpfilename;
	delete[] swcpfilename;
	delete[] filename;
	for(i=0;i<2;i++)
		delete[] sucessfulmatch[i];
//	delete sucessfulmatch;
	delete[] swlp_single;
	delete[] swcp_single;
	delete[] swcp_match;
	return 0;
}

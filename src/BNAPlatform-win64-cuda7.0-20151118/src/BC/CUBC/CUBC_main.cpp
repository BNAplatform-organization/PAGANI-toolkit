#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <iomanip>
#include "dirent.h"

using namespace std;

void CUBC(float *BC, int * row, int * col, int numVertices, int numEdges);

int main(int argc, char** argv) 
{
	// input 1: file name
	// input 2 & 3£º for tuning performance
	//int grid = atoi(argv[2]);
	//int thread = atoi(argv[3]);

	clock_t total_time = clock();
	if (argc < 2) 
	{
		cerr<<"Input format: .\\CUBC.exe dir_for_csr \nFor example: .\\CUBC.exe d:\\data "<<endl;
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
		// read input file

		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\ncalculating Betweenness for "<<a.c_str()<<" ..."<<endl;
		ifstream fCSR(a.c_str(), ios_base::binary);
		if (!fCSR.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}

		int numVertices;	// Network vertices number
		int numEdges;	//Network egdes number
		int *row = NULL;	// Cormat CSR row
		int *col = NULL;	// Cormat CSR col

		fCSR.read((char *)&numVertices, sizeof(int));
		numVertices = numVertices - 1;	
		row = new int[numVertices+1];	// Cormat CSR row
		fCSR.read((char*)row, sizeof(int)*(numVertices + 1));
		fCSR.read((char *)&numEdges, sizeof(int));
		col = new int[numEdges];	// Cormat CSR col
		fCSR.read((char*)col, sizeof(int)*numEdges);
		fCSR.close();

		float *BC = new float[numVertices];
		CUBC(BC, row, col, numVertices, numEdges);

		// Parse file name
		string X_bc = a.substr(0, a.find_last_of('.') ).append("_bc.nm");
		cout<<"Save betweenness data for each node as "<<X_bc.c_str()<<endl;
		ofstream fout;
		fout.open(X_bc.c_str(), ios::binary|ios::out);
		fout.write((char*)&numVertices, sizeof(int));
		fout.write((char*)BC, sizeof(float) * numVertices);
		fout.close();

		delete []BC;
		delete []row;
		delete []col;
	}
	delete []filename;
	total_time = clock() - total_time;
	cout<<"total elapsed time: "<<1.0*total_time/1000<<" s."<<endl;
	cout<<"==========================================================="<<endl;
	return 0;
}
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <ctime>

#include "BC_CPU.h"

using namespace std;

void read_csr(char * fileCSR, int &numVertices, int &numEdges, int *&row, int *&col)
{
	string csrfilename = string(fileCSR);
	ifstream fCSR(csrfilename.c_str(), ios::binary);
	if(!fCSR.good())
	{
		cout<<"cannot open .csr file!"<<endl<<"note: please add \".csr\" at the end of the file name "<<endl;
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
	return;
}

int main( int argc, char** argv) 
{
	// input 1: file name

	// read input .csr file 
	int numVertices;	// Network vertices number
	int numEdges;	//Network egdes number
	int *row = NULL;	// Cormat CSR row
	int *col = NULL;	// Cormat CSR col
	read_csr(argv[1], numVertices, numEdges, row, col);

	float *BC_GOLDEN = new float[numVertices];

	// computing Betweenness
	cout<<"Start BC computing on CPU. n = "<<numVertices<<", m = "<<numEdges<<endl;	
	clock_t cpu_time = clock();
	Betweenness_CPU(row, col, numVertices, numEdges, BC_GOLDEN);
	cpu_time = clock() - cpu_time;
	cout<<cpu_time<<" ms for CPU BC"<<endl;
	cout<<(double)numEdges*numVertices*1000/cpu_time/1024/1024<<" MTEPS"<<endl;
	
	// write to file
	ofstream fout_golden("BC_gloden.txt");
	for (int i = 0; i < numVertices; i++)
		fout_golden<<BC_GOLDEN[i]<<endl;
	fout_golden.close();
	
	// record running time
	ofstream cputimelog("cputimelog.txt", ios::app);
	cputimelog<<cpu_time<<" ms."<<endl;
	cputimelog.close();
	
	delete []BC_GOLDEN;
	delete []row;
	delete []col;
}
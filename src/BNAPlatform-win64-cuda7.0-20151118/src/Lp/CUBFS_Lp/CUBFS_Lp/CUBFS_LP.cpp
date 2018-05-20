#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <ctime>

#include "BFS_APSP_GPU.h"
using namespace std;

double CUBFS_Lp(int *row, int *col, int numVertices, int numEdges)
{
	int grid = 128;
	int thread = 64;

	long long size = (long long)numVertices*numVertices;
	int *dist = new int[size];

	// computing APSP
	cout<<"Start BFS-based APSP computing on GPU. n = "<<numVertices<<", m = "<<numEdges<<endl;
	clock_t gpu_time = clock();

	APSP_GPU(dist, row, col, numVertices, numEdges, grid, thread);

	gpu_time = clock() - gpu_time;
	cout<<gpu_time<<" ms for GPU APSP"<<endl;
	cout<<(double)numEdges*numVertices/gpu_time*1000/1024/1024<<" MTEPS"<<endl;

	// computing Lp
	double Lp_result = 0.0;
	long long  i;
	for (i = 0; i < size; i++)
	{
		if (dist[i] > 0)
		{
			//APSP_output[i * numVertices + j] = INFINITE;//__FLT_MAX__;
			Lp_result += 1.0f / dist[i];
		}
	}
	Lp_result /= numVertices;
	Lp_result /= numVertices - 1;
	Lp_result = 1 / Lp_result;

	// write to file
	//ofstream fout("Lp.txt");
	//fout<<Lp_result<<endl;
	//fout.close();

	// record running time
	//ofstream timelog("gputimelog.txt", ios::app);
	//timelog<<"grid = "<<grid<<", thread = "<<thread<<", Time = "<<gpu_time<<" ms."<<endl;
	//timelog.close();

	delete []dist;
	return Lp_result;
}

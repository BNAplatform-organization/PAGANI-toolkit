# include <iostream>
# include <fstream>
# include <ctime>
# include "APSP_BFS.h"
# include "FormBlock.h"
# include "BFW.h"

using namespace std;

void cuAPSP(float *costmat, const int numVertices, const int block_size);
extern float *Li_result;

double CUBFW_Lp(int *row, int *col, float *power, int numVertices, int numEdges)
{
	const int sizeBlock = 96*32;
	long long squareVertices = numVertices;
	squareVertices *= numVertices;
	// Padding for block-Floyd-Warshell (non-sparse, non-symmetric)
	int numVerticesPaded = (numVertices + sizeBlock - 1);
	numVerticesPaded /= sizeBlock;
	numVerticesPaded *= sizeBlock;

	int *rowPaded = new int[numVerticesPaded + 1];
	memcpy((void*)rowPaded, (void*)row, sizeof(int) * (numVertices + 1));
	for (int k = numVertices + 1; k < numVerticesPaded + 1; k++)
		rowPaded[k] = row[numVertices];

	int	cntBlock = numVerticesPaded / sizeBlock;
	long long squareVerticesPaded = cntBlock * cntBlock;
	squareVerticesPaded *= sizeBlock;
	squareVerticesPaded *= sizeBlock;

	float *costmatPaded;
	costmatPaded = new float[squareVerticesPaded];
	//	cudaHostAlloc((void **)&costmatPaded, sizeof(float)* squareVerticesPaded, cudaHostAllocDefault);
	if (costmatPaded == NULL)
	{
		cout<<"Allocation failure"<<endl;
		return 0;
	}


	printf( "Blocked FW algorithm on GPU with %d block(s)\n", cntBlock);
	printf( "Start calculating APSP...\n");

	init_block(costmatPaded, rowPaded, col, power, numVerticesPaded, sizeBlock);//(symmetric and) sparse graph -> full dense blocked graph  numVertices? or numVerticesPaded?
	cuAPSP(costmatPaded, numVerticesPaded, sizeBlock);
	//BFW_C	(costmatPaded, (long long) numVerticesPaded, (long long) sizeBlock);
	printf( "End calculating APSP...\n");

	//BFW_C(costmatPaded, numVerticesPaded, sizeBlock);

	//float *dist = new float[squareVertices];
	//if (dist == NULL)
	//{
	//	cout<<"Allocation failure"<<endl;
	//	return 1;
	//}
	//printf( "Start postblocking...\n");
	//post_block(dist, costmatPaded, numVertices, sizeBlock);
	//printf( "End postblocking...\n");

	// computing Lp
	double Lp_result = 0.0;
	long long index = 0;
	memset(Li_result, 0, sizeof(float)*numVertices);
	int  isolated_vertices = 0;
	bool unconnected_mark = false;

	for (int i = 0; i < cntBlock; i++)
	{
		for (int j = 0; j < cntBlock; j++)
		{
			for (int ii = 0; ii < sizeBlock; ii++)
			{
				for (int jj = 0; jj < sizeBlock; jj++)
				{
					int u = i * sizeBlock + ii;
					int v = j * sizeBlock + jj;
					if (!(u == v) && (u < numVertices) && (v < numVertices))
					{
						if (costmatPaded[index]<1e-1)
							;
							//cout<<"ii: "<<ii<<"\tjj: "<<jj<<"\ti: "<<i<<"\tj: "<<j<<"\tindex: "<<index<<"\tcostmatPaded[index]:"<<costmatPaded[index]<<endl;
						else
						{
							Lp_result += 1.0f /costmatPaded[index];
							Li_result[u] += 1.0f /costmatPaded[index];
						}
					}
					index++;
				}
			}
		}
	}

	for (int i = 0; i < numVertices; i++)
	{			
		//if (i > 10000 && i < 10020)
		//cout<<Li_result[i]<<endl;
		Li_result[i] /= (numVertices - 1); 

		if(Li_result[i] < 1.0f/numVertices)
		{
			Li_result[i] = 0;
			isolated_vertices++;
			//cout<<i<<"th node is isolated;\n";
		}
	}

	Lp_result /= numVertices;
	Lp_result /= (numVertices - 1);
	Lp_result = 1 / Lp_result;

	if (isolated_vertices > 0 )
		cout<<"\nisolated vertices number : "<<isolated_vertices<<endl;
	else if (unconnected_mark)
		cout<<"\nThe network is unconnected. "<<endl;

	//// write to file
	//ofstream fout("Lp.txt");
	//fout<<Lp_result<<endl;
	//cout<<"Lp: "<<Lp_result<<endl;
	//fout.close();

	//// verify
	//// CPU
	//float *output_C = new float[squareVertices];
	//if (output_C == NULL)
	//{
	//	cout<<"Allocation failure"<<endl;
	//	return 1;
	//}
	//APSP_BFS(output_C, row, col, numVertices);

	//cout<<"comparing..."<<endl;
	//for (int i = 0; i < numVertices; i++)
	//	for (int j = 0; j < numVertices; j++)
	//	{
	//		if (!(fabs(output_C[i * numVertices + j])<1e-12))
	//		{
	//			if ((fabs(output_C[i * numVertices + j]-dist[i * numVertices + j]) > 1e-12))
	//			{
	//				cout<<"cuole\t"<<"\t"
	//					<<"i\t"<<i<<'\t'
	//					<<"j\t"<<j<<"\t"
	//					<<output_C[i * numVertices + j]<<"\t"
	//					<<dist[i * numVerticesPaded + j]
	//				<<endl;

	//				getchar();
	//			}
	//		}
	//	}
	//	cout<<"niule_gpu"<<endl;

	//	delete []output_C;

	delete []costmatPaded;
	delete []rowPaded;
	return Lp_result;
}


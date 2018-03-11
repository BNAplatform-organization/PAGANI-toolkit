#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "dirent.h"
#include <list>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc != 2 && argc != 3)
	{
		cout << "No Directory Input" << endl;
		return -1;
	}
	

	DIR *dp;
	Dirent *dirp;	
	list<string> fileNameList;
	dp = opendir(argv[1]);
	if(!dp)
	{
		cout << "can't open " << argv[1] << endl;
		return -1;
	}
	
	if (argv[1][strlen(argv[1]) - 1] == '\\')
		argv[1][strlen(argv[1]) - 1] = 0;
	
	string fileNametemp;
	while((dirp = readdir(dp)) != NULL)
	{
		fileNametemp = string(dirp->d_name);
		fileNameList.push_back(fileNametemp);
	}

	closedir(dp);
	
	list<string> fileOnProcess;

	for (list<string>::iterator iter = fileNameList.begin(); iter != fileNameList.end(); ++iter)
	{
		int pos = iter->find(".csr");
		if (pos != string::npos)
		{
			iter->replace(pos, 4, "");
			fileOnProcess.push_back(*iter);
		}
	}
	

	for (list<string>::iterator iter = fileOnProcess.begin(); iter != fileOnProcess.end(); ++iter)
	{
		string fileName = string(argv[1]) + string("\\")+(*iter);
		string fileName_csr = fileName + string(".csr");
		string fileName_modu = fileName + string("_modu.nm");
		string fileName_pc = fileName + string("_pc.nm");
		string fileName_pc_Normalized = fileName_pc + string("_Normalized.nm");
		cout<<"calculate participant coefficient of "<<fileName_csr<<endl<<endl;
		//data input
		ifstream infile;
		infile.open(fileName_csr.c_str(),ios::in|ios::binary);
		if (!infile)
		{
			cout << fileName_csr << " not found!" << endl;
			return 0;
		}

		int NumVertex;
		infile.read((char*)&NumVertex, sizeof(int));
		int* EdgeIndex = new int[NumVertex];
		infile.read((char*)EdgeIndex, sizeof(int) * NumVertex);
		int NumEdges;
		infile.read((char*)&NumEdges, sizeof(int));
		int* EdgeTargets = new int[NumEdges];
		infile.read((char*)EdgeTargets, sizeof(int) * NumEdges);
		int Vlength;
		infile.read((char*)&Vlength, sizeof(int));
		float * V = new float [Vlength];
		infile.read((char*)V, sizeof(float) * Vlength);
		infile.close();


		infile.open(fileName_modu.c_str(), ios::in|ios::binary);
		if (!infile)
		{
			cout << fileName_modu << " not found!" << endl;
			continue;
		}
		infile.read((char*)&NumVertex, sizeof(int));
		float* module = new float[NumVertex];
		infile.read((char*)module, sizeof(float) * NumVertex);
		infile.close();

		//data output
		float* original, *normalized;
		original = new float[NumVertex];

		int NumModu = 0;
		for (int i = 0; i < NumVertex; ++i)
		{
			if((int) module[i] > NumModu)
			{
				NumModu = (int) module[i];
			}
		}
        //++NumModu;

		float* VertexDegree = new float[NumVertex];
		memset(VertexDegree, 0, NumVertex*sizeof(float));
		for (int i = 0; i < NumVertex; ++i)
		{
			for(int j = EdgeIndex[i]; j < EdgeIndex[i+1]; j++)
			VertexDegree[i] += V[j];
		}

		float* WithinModule = new float [NumModu];
		for (int i = 0; i < NumVertex; ++i)
		{
			if (VertexDegree[i] == 0)
			{
				original[i] = 0.0f;
			}
			else
			{
                /*for (int j = 1; j <= NumModu; ++j)
				{
					int sum = 0;
					for (int k = EdgeIndex[i]; k < EdgeIndex[i+1]; ++k)
					{
						if (module[EdgeTargets[k]] == j)
						{
							++sum;
						}
					}
					WithinModule[j] = sum;
				}
				*/
				memset(WithinModule,0,NumModu*sizeof(int));
				
				int m;
				for (int k = EdgeIndex[i]; k < EdgeIndex[i+1]; ++k)
					{
						m = (int) module[EdgeTargets[k]];
						WithinModule[m-1] += V[k];
					}
				


				float sum = 0.0f;
				for (int j = 0; j < NumModu; ++j)
				{
					float temp = float(WithinModule[j]) / float(VertexDegree[i]);
					sum += temp * temp;
				}
				original[i] = 1.0f - sum;
			}
		}

		if (NumModu != 1)
		{
			float max = float(NumModu - 1) / float(NumModu);
			normalized = new float[NumVertex];
			for (int i = 0; i < NumVertex; ++i)
			{
				normalized[i] = original[i] / max;
			}
		}
		else
		{
			normalized = original;
		}
		ofstream outfile;
		outfile.open(fileName_pc.c_str(),ios::out|ios::binary);
		outfile.write((char*)&NumVertex, sizeof(int));
		
		//if (argc == 3 && argv[2][0] == 'n' )
		//	outfile.write((char*)normalized, sizeof(float) * NumVertex);
		//else 
			outfile.write((char*)original, sizeof(float) * NumVertex);
		
		outfile.close();

		outfile.open(fileName_pc_Normalized.c_str(),ios::out|ios::binary);
		outfile.write((char*)&NumVertex, sizeof(int));
		outfile.write((char*)normalized, sizeof(float) * NumVertex);
		
		outfile.close();

		delete [] EdgeIndex;
		delete [] EdgeTargets;
		delete [] V;
		delete [] module;
		delete [] original;
		delete [] VertexDegree;
		delete [] WithinModule;
		if (NumModu != 1)
		{
			delete [] normalized;
		}
	}

	return 1;
}

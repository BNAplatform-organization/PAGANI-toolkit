# include <iostream>
# include <fstream>
# include <ctime>
# include <cstring>
# include <stdlib.h>
# include "Maslov.h"
using namespace std;

//void Maslov_weighted_1(int * R_dst, int * C_dst, float * V_dst, int * R_src, int * C_src, float *V_src, int Rlength, int Clength);

int main(int argc, char * argv[])
{
/*
input:
	argv[1] :	network name, exclude the suffix
*/
	//cout<<RAND_MAX<<endl;
	if ( argc != 2 ) 
	{
		printf("1 argument is required.\nInput format example: .\\Maslov.exe ..\\..\\data\\data_avr_20_0.06\n");
		exit(1);
	}
	int i = 0, j = 0, k = 0;

	// Parse file name
	int len_arg1 = strlen(argv[1]);

	char * OutCSR = new char [len_arg1 + 8]; 
	strcpy(OutCSR, argv[1]);
	strcpy(OutCSR + len_arg1 - 4, "_Maslov.csr");
	OutCSR[len_arg1 + 7] = 0;

	// Read CSR file
	ifstream fin(argv[1], ios::binary);
	if (!fin.good())
	{	cout<<"Can't open\t"<<argv[1]<<endl;	system("pause"); return 0;}
	cout<<"Reading CSR file..."<<endl;
	int Rlength = 0, Clength = 0, Vlength = 0;
	fin.read((char*)&Rlength, sizeof(int));
	int * R = new int [Rlength];
	fin.read((char*)R, sizeof(int) * Rlength);
	fin.read((char*)&Clength, sizeof(int));
	int * C = new int [Clength];
	fin.read((char*)C, sizeof(int) * Clength);
	fin.read((char*)&Vlength, sizeof(int));
	float * V = new float [Vlength];
	fin.read((char*)V, sizeof(int) * Vlength);
	fin.close();
	int N = Rlength - 1;

	clock_t time;
	time = clock();
	// Construct Adjacent Graph
	cout<<"Rewiring..."<<endl;	
	int * R_dst = new int [Rlength];
	int * C_dst = new int [Clength];
	float * V_dst = new float [Vlength];
	Maslov_weighted_1(R_dst, C_dst, V_dst, R, C, V, Rlength, Clength);

	time = clock() - time;
	cout<<"Elapsed time:\t"<<((float)time)/1000<<"s"<<endl;
	// Output to new CSR file
	ofstream fout(OutCSR, ios::binary | ios::out);
	if (!fin.good())
	{	cout<<"Can't open\t"<<OutCSR<<endl;	return 0;}
	fout.write((char*)&Rlength, sizeof(int));
	fout.write((char*)R_dst, sizeof(int) * Rlength);
	fout.write((char*)&Clength, sizeof(int));
	fout.write((char*)C_dst, sizeof(int) * Clength);
	fout.write((char*)&Vlength, sizeof(int));
	fout.write((char*)V_dst, sizeof(int) * Vlength);
	fout.close();	 
	cout<<"==========================================================="<<endl;
	delete []R;
	delete []C;
	delete []V;
	delete []R_dst;
	delete []C_dst;
	delete []V_dst;
	return 0;
}


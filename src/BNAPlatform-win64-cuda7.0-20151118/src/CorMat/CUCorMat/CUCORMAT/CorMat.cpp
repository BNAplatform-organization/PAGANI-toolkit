
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <ctime>
#include <cmath>
#include <iomanip>
#include "dirent.h"
#include <iostream>
#include <fstream>
#include <stack>

using namespace std;

typedef float real__t;
typedef unsigned int uint__t;

const int HdrLen = 352;
double ProbCut = 0.5; 
const int blocksize = 1024*1024*48;

bool rs_flag = true;
bool cormat_flag = false;
char mask_dt;

void CorMat_cpu(real__t * Cormat, real__t * BOLD, int N, int L);
int CorMat_gpu(real__t * Cormat, real__t * BOLD, int N, int L, int Batch_size);
long long FormFullGraph(bool * adjacent, real__t * Cormat, int N, real__t threshold);
real__t FormFullGraph_s(bool * adjacent, real__t * Cormat, int N, long long threshold);
void post_block (real__t * Cormat, real__t * Cormat_blocked, int dim, int block_size,bool fishtran_flag);
void FormCSRGraph(int * R, int * C, bool * adjacent, int N);
long long find_max(real__t *Cormat, long long M1);
real__t select_GPU(real__t *Cormat, long long M1, long long k);

real__t fisher_trans(real__t r)
{
	real__t z;
	if (r==1) r -= 1e-6; 
	z = 0.5*log((1+r)/(1-r));
	return z;	
}

real__t inv_fisher_trans(real__t z)
{
	real__t r;
	r = exp(2*z);
	r = (r-1)/(r+1);
	return r;	
}

int main(int argc, char * argv[])
{
	clock_t total_time = clock();
	if (argc < 6) 
	{
		cerr<<"Input format: .\\CorMat.exe  Dir_for_BOLD  threshold_for_mask(0~1) to_average(yf/yn/bf/bn/n) to_save_cormatrix(y/n) threshold_type(r/s) threshold_for_correletaion_coefficient(s)(0~1)\n"
			<<"For example: .\\CorMat.exe  d:\\BOLD  0.5 y n r 0.2 0.25 0.3\n"<<endl;
		exit(1);	
	}
	int L, N = 0, i = 0, j = 0, k = 0, l = 0, total_size;
	clock_t time;

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
		//cout<<filenametmp.c_str()<<"!"<<endl;
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".nii") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			if (filenametmp.compare("mask.nii")!=0)
				FileNumber++;
		}
	}
	cout<<FileNumber<<" file(s) to be processed."<<endl;

	closedir(dp);
	string *filename = new string[FileNumber];
	dp = opendir(argv[1]);
	i = 0;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		//	cout<<"here";
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && filenametmp.substr(filenametmp.find_last_of('.'),4).compare(".nii") == 0 && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 3)
		{
			if (filenametmp.compare("mask.nii")!=0)
				filename[i++] = filenametmp;
		}
	}

	real__t ProbCut = (real__t)atof(argv[2]);
	int NumS = argc - 6;
	real__t * r_thresh = new real__t [NumS];
	real__t * s_thresh = new real__t [NumS];
	if (argv[5][0] == 'r' || argv[5][0] == 'R' )
		for (i = 0; i < NumS; i++)
			r_thresh[i] = (real__t)atof(argv[6+i]);
	else if (argv[5][0] == 's' || argv[5][0] == 'S' )
	{
		rs_flag = false;
		memset(r_thresh, 0, sizeof(real__t)*NumS);
		for (i = 0; i < NumS; i++)
			s_thresh[i] = (real__t)atof(argv[6+i]);
	}
	else {
		cout << "threshold type error! \nr for correlation threshold\ns for spasity threshold\n";
		exit(1);
	}

	if(argv[4][0]=='y' || argv[4][0]=='Y')
	{
		cormat_flag = true;
	}
	else if (argv[4][0]=='N' || argv[4][0]=='n')
	{
		cormat_flag = false;
	}
	else
	{
		cout << "to_save_cor_matrix type error! \ny to save the whole correlation matrix \nn to save only csr format of adjacency matrix\n";
		exit(1);
	}
	// read input files and parameters
	if (argv[1][strlen(argv[1]) - 1] == '\\')
		argv[1][strlen(argv[1]) - 1] = 0;
	ifstream fin(string(argv[1]).append("\\mask.nii").c_str(), ios::binary);
	if (!fin.good())
	{	cout<<"Can't open\t"<<string(argv[1]).append("\\mask.nii").c_str()<<endl;	return 0;}
	short hdr[HdrLen / 2];
	fin.read((char*)hdr, HdrLen);
	cout<<"mask datatype : "<< hdr[35]<<"  "<<hdr[36]<<endl;
	mask_dt = hdr[35];
	//typedef unsigned char mask_type; 
	total_size = hdr[21] * hdr[22] * hdr[23];	// Total number of voxels
	
	real__t * mask = new float [total_size];

	if (mask_dt==2)
	{
		unsigned char *mask_uc = new unsigned char[total_size];
		fin.read((char *) mask_uc, sizeof(unsigned char) * total_size);
		for (int vm = 0; vm<total_size; vm++)
			mask[vm] = (float) mask_uc[vm];
		delete [] mask_uc;
	}
	else if(mask_dt==16)
	{
		fin.read((char *)mask, sizeof(float) * total_size);
	}
	else
	{
		cout<<"mask data-type error, Only the type of unsigned char and float can be handled.\n";
		//system("pause");
		return -1;
	}

	fin.close();
	// Count the number of the valid voxels	
	for (k = 0; k < total_size; k++)
		N += ( (float) mask[k] >= ProbCut);	
	cout<<"Data size: "<<hdr[21] <<"x"<<hdr[22]<<"x"<<hdr[23]  <<", Grey voxel count: "<<N<<"."<<endl;
	// swap the largest threshold to the beginning
	int min_idx = 0;
	for (i = 0; i < NumS; i++)
		if (r_thresh[i] < r_thresh[min_idx])
			min_idx = i;
	real__t temp = r_thresh[0];
	r_thresh[0] = r_thresh[min_idx];
	r_thresh[min_idx] = temp;
	// CSR format
	int Rlength = N + 1;
	long long Clength;
	int * R = new int[Rlength];
	int * C;
	
	//process, do not average
	if (argv[3][0] == 'n' || argv[3][0] == 'N' || argv[3][0] == 'b' || argv[3][0] == 'B' )
		for (int i = 0; i < FileNumber; i++)
		{
			string a = string(argv[1]).append("\\").append(filename[i]);
			cout<<"\ncalculating correlation for "<<a.c_str()<<" ..."<<endl;
			ifstream fin(a.c_str(), ios_base::binary);

			// Get the length of the time sequence 
			if (!fin.good())
			{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}
			fin.read((char*)hdr, HdrLen);
			L = hdr[24];
			real__t * BOLD = new real__t [L * N];
			if (hdr[36] == 64) // double
			{
				double *  InData = new double [L * total_size];
				fin.read((char*)InData, sizeof(double) * L * total_size);
				fin.close();
				// Get the BOLD signal for all the valid voxels
				for (int i = -1, k = 0; k < total_size; k++)
					if ((float) mask[k] >= ProbCut)
						for (i++, l = 0; l < L; l++)
						{
							BOLD[l*N+i] = InData[l*total_size+k];
						}
						cout<<"BOLD length: "<<L<<", Data type: double."<<endl;
						delete []InData;
			}
			else if (hdr[36] == 32)	   //float
			{
				real__t *InData = new float [L * total_size];
				fin.read((char*)InData, sizeof(float) * L * total_size);
				fin.close();
				// Get the BOLD signal for all the valid voxels
				for (int i = -1, k = 0; k < total_size; k++)
					if ((float) mask[k] >= ProbCut)
						for (i++, l = 0; l < L; l++)
						{
							BOLD[l*N+i] = InData[l*total_size+k];
						}
						cout<<"BOLD length: "<<L<<", Data type: float."<<N<<endl;
						delete []InData;
			}
			else
			{
			  cerr<<"Error: Data type is neither float nor double."<<endl;
			}
			
			// set some parameters 
			int Batch_size = 1024 * 3;				// should not be smaller than 1024 !
			int Num_Blocks = (N + Batch_size - 1) / Batch_size;
			long long M2 = Num_Blocks * (Num_Blocks + 1) / 2;
			M2 *= Batch_size * Batch_size;
			real__t * Cormat_blocked = new real__t [M2];

			long long M1 = (N-1);
			M1 *= N;
			M1 /= 2;
			real__t * Cormat_gpu = new real__t [M1];

			// Begin computing correlation	
			time = clock();
			CorMat_gpu(Cormat_blocked, BOLD, N, L, Batch_size);
			
			time = clock() - time;
			memset((void*)Cormat_gpu, 0, sizeof(real__t) * M1);
			post_block(Cormat_gpu, Cormat_blocked, N, Batch_size,0);		
			cout<<"GPU correlation time: "<<time<<"ms"<<endl;
		
//					CorMat_cpu(Cormat_cpu, BOLD, N, L);

			//CorMat_gpu(Cormat_blocked, BOLD, N, L, Batch_size);
			//time = clock() - time;
			//cout<<"GPU correlation time: "<<time<<"ms"<<endl;
			//time = clock();
			//memset((void*)Cormat_gpu, 0, sizeof(real__t) * M1);
			//cout<<"Postprocessing ..."<<endl;
			//post_block(Cormat_gpu, Cormat_blocked, N, Batch_size);

			/*int kkk = 0;
			for (int i = 0; i < N; i++)
				for (int j = i+1; j < N; j++)
				{
						if (fabs(Cormat_gpu[kkk] - Cormat_cpu[kkk]) > 1e-6)
					{
						cout<<"cuole"<<'\t'<<Cormat_gpu[kkk]<<'\t'<<Cormat_cpu[kkk]<<'\t'<<i<<'\t'<<j<<endl;
						getchar();
					}
					kkk++;
				}*/

			delete []Cormat_blocked;

			char sparsity[30];
			char Graph_size[30];
			string OutCor = a.substr(0, a.find_last_of('.')).append("_").append(string(itoa(N, Graph_size, 10)));
			if(cormat_flag == true)
			{
				string cormat_filename = OutCor;
				cormat_filename.append(".cormat");
				ofstream cormat_file;
				cormat_file.open(cormat_filename.c_str(), ios::binary | ios::out);
				cormat_file.write((char*)&M1, sizeof(int));
				cormat_file.write((char*)Cormat_gpu,sizeof(real__t)*M1);
				cormat_file.close();
			}

			string Outfilename;
			ofstream fout;

			bool * adjacent = new bool [(long long)N*N];
			for (k = 0; k < NumS; k++)
			{
				// Form full graph
				time = clock();
				//cout<<"threshold for r = "<<r_thresh[k]<<endl;
				if (rs_flag){
					Clength = FormFullGraph(adjacent, Cormat_gpu, N, r_thresh[k]);
					s_thresh[k] = 100.0 * Clength / M1 / 2;
				}
				else {
					cout<<"\nthreshold for sparsity = "<<s_thresh[k]<<endl;
					Clength =  (long long) 2*M1*s_thresh[k]/100.0;
					Clength += Clength%2;
					r_thresh[k] = FormFullGraph_s(adjacent, Cormat_gpu, N, Clength/2);
				}
				//cout<<"Clength = "<<Clength<<endl;
			
				time = clock() - time;
						cout<<"time for forming full graph (ms):\t" <<time<<endl;
				//if (k == 0)
					C = new int [Clength];
	
				// form csr graph
				//time = clock();
				FormCSRGraph(R, C, adjacent, N);
				if(R[N]!=Clength) {cout<<"R[N] != Clength"<<endl; system("pause");}
				//time = clock() - time;
				//		cout<<"time for forming csr graph (ms):\t" <<time<<endl;

				// write files

			    //for (int i=0;i<Rlength;i++)
			    //	   cout<<R[i]<<' ';
				sprintf_s(sparsity, "_spa%.3f%%_cor%.3f", s_thresh[k],r_thresh[k]);
				Outfilename = OutCor;
				Outfilename.append(sparsity).append(".csr");
				cout<<"generating "<<Outfilename.c_str()<< "..."<<endl;
				fout.open(Outfilename.c_str(), ios::binary | ios::out);
				fout.write((char*)&Rlength, sizeof(int));
				fout.write((char*)R, sizeof(int) * Rlength);
				fout.write((char*)&Clength, sizeof(int));
				fout.write((char*)C, sizeof(int) * Clength);
				fout.close();

				delete []C;	
				
			}
			/*
			long long hisogram[21] = {0};
			real__t avg=0;
			int illegal=0;
			for (long long z=0; z<M1; z++){
				    if (Cormat_gpu[z]<1000 || Cormat_gpu[z]>-1000)
					{
					  hisogram[(int) ((Cormat_gpu[z]+1.05)*10)]++;
					  avg += Cormat_gpu[z];
					}
					else illegal++;
			}
			avg = avg/M1;
			ofstream fo;
			fo.open("E:\\xumo\\onenii\\hisogram.txt");
			
			fo <<"\n=====================================\n";
			for (int num=0; num<21;num++){
				fo<<(100.0*(double) hisogram[num])/M1<<endl;								
			} 
			fo<<"\n=====================================\n";
			fo<<endl<<avg;
			fo<<endl<<100.0*((double) illegal)/M1;
			fo<<"\n=====================================\n";
			fo.close();
			*/

				
			delete []adjacent;
			delete []Cormat_gpu;
			delete []BOLD;
		}

		if (argv[3][0] == 'y' || argv[3][0] == 'Y' || argv[3][0] == 'b' || argv[3][0] == 'B' )
		{
			// set some parameters 
			int Batch_size = 1024 * 3;				// should not be smaller than 1024 !
			int Num_Blocks = (N + Batch_size - 1) / Batch_size;
			long long M2 = Num_Blocks * (Num_Blocks + 1) / 2;
			M2 *= Batch_size * Batch_size;
			real__t * Cormat_blocked = new real__t [M2];

			long long M1 = (N-1);
			M1 *= N;
			M1 /= 2;
			real__t * Cormat_gpu = new real__t [M1];
			memset((void*)Cormat_gpu, 0, sizeof(real__t) * M1);

			//real__t * Cormat_blocked_ = new real__t [M2];
			//real__t * Cormat_gpu_ = new real__t [M2];

			for (int i = 0; i < FileNumber; i++)
			{
				//check whether the first int are identical!
				string a = string(argv[1]).append("\\").append(filename[i]);
				cout<<"\ncalculating correlation for "<<a.c_str()<<" ..."<<endl;
				ifstream fin(a.c_str(), ios_base::binary);

				// Get the length of the time sequence 
				if (!fin.good())
				{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}
				fin.read((char*)hdr, HdrLen);
				L = hdr[24];
				real__t * BOLD = new real__t [L * N];
				if (hdr[36] == 64) // double
				{
					double *  InData = new double [L * total_size];
					fin.read((char*)InData, sizeof(double) * L * total_size);
					fin.close();
					// Get the BOLD signal for all the valid voxels
					for (int i = -1, k = 0; k < total_size; k++)
						if ((float)mask[k] >= ProbCut)
							for (i++, l = 0; l < L; l++)
							{
								BOLD[l*N+i] = InData[l*total_size+k];
							}
							cout<<"BOLD length: "<<L<<", Data type: double."<<endl;
							delete []InData;
				}
				else if (hdr[36] == 32)	   //float
				{
					real__t *InData = new float [L * total_size];
					fin.read((char*)InData, sizeof(float) * L * total_size);
					fin.close();
					// Get the BOLD signal for all the valid voxels
					for (int i = -1, k = 0; k < total_size; k++)
						if ((float) mask[k] >= ProbCut)
							for (i++, l = 0; l < L; l++)
							{
								BOLD[l*N+i] = InData[l*total_size+k];
							}
							cout<<"BOLD length: "<<L<<", Data type: float."<<N<<endl;
							delete []InData;
				}
				else
				{
					cerr<<"Error: Data type is neither float nor double."<<endl;
				}

				// Begin computing correlation
				time = clock();
			//	if (i == 0)
				{	   CorMat_gpu(Cormat_blocked, BOLD, N, L, Batch_size);
						//CorMat_cpu(Cormat_gpu, BOLD, N, L);
				}
			//	else
				{
					//	CorMat_gpu(Cormat_blocked_, BOLD, N, L, Batch_size);
					//CorMat_cpu(Cormat_gpu_, BOLD, N, L);
				}
				
				time = clock() - time;
				cout<<"GPU correlation time: "<<time<<"ms"<<endl;
				if (argv[3][1] == 'f' || argv[3][1] == 'F')
					post_block(Cormat_gpu, Cormat_blocked, N, Batch_size,1);
				else if (argv[3][1] == 'n' || argv[3][1] == 'N')
					post_block(Cormat_gpu, Cormat_blocked, N, Batch_size,0);
			  /*
				for (int i = 0; i < M1; i++)
				{
					Cormat_gpu[i] = i;
				}
				int id = 0;
				for (int i = 0; i < (N + Batch_size-1)/Batch_size; i++)
				{
					for (int j = i+1; j < (N + Batch_size-1)/Batch_size; j++)
					{
						for (int ii = 0; ii < Batch_size ; ii++)
						{
							for (int jj = 0; jj < Batch_size; jj++)
							{
								Cormat_blocked[id]=(float)((long long)i * );
								id ++;
							}
						}
					}
				}		   */
							
				delete []BOLD;
			}
			/*int count = 0;
			for (int i = 0; i < M2; i++)
			{
				if(fabs(Cormat_blocked[i] - Cormat_blocked_[i])>0.2)
					{
						//cout<<i<<"\t"<<Cormat_blocked[i]<<"\t"<<Cormat_blocked_[i]<<endl;
						count ++;//system("pause");
				}
				Cormat_blocked[i] +=  Cormat_blocked_[i];
			}	
			cout<<"count = "<<count<<endl;
			
			
			for (int i = 0; i < M1; i++)
			{
				Cormat_gpu[i] += Cormat_gpu_[i];
			}		*/
			
			delete []Cormat_blocked;
			
			if (argv[3][1] == 'f' || argv[3][1] == 'F')
				for (int i = 0; i < M1; i++)
				{
					Cormat_gpu[i] /= FileNumber;
					Cormat_gpu[i] = inv_fisher_trans(Cormat_gpu[i]);
				}
			else
				for (int i = 0; i < M1; i++)
				{
					Cormat_gpu[i] /= FileNumber;					
				}

			char sparsity[30];
			char Graph_size[30];
			string b = string(argv[1]);
			string OutCor;
			if (b.find_last_of('\\') ==  -1)
				OutCor = b.append("\\");
			else
				OutCor = b.append(b.substr(b.find_last_of('\\'), b.length() - b.find_last_of('\\')));
			OutCor.append("_").append(string(itoa(FileNumber, Graph_size, 10))).append("_").append(string(itoa(N, Graph_size, 10)));
			
			if(cormat_flag == true)
			{
				string cormat_filename = OutCor;
				cormat_filename.append(".cormat");
				ofstream cormat_file;
				cormat_file.open(cormat_filename.c_str(), ios::binary | ios::out);
				cormat_file.write((char*)&M1, sizeof(int));
				cormat_file.write((char*)Cormat_gpu,sizeof(real__t)*M1);
				cormat_file.close();
			}

			string Outfilename;
			ofstream fout;

			//// output the dense float corr matrix to compute p values 
			//fout.open("..\\..\\..\\..\\..\\corr_mat", ios::binary | ios::out);
			//fout.write((char*)&N, sizeof(int));
			//fout.write((char*)Cormat_gpu, sizeof(float) * M1);
			//fout.close();

			bool * adjacent = new bool [(long long)N*N];
			for (k = 0; k < NumS; k++)
			{
				// Form full graph
				time = clock();
				if (rs_flag){
					Clength = FormFullGraph(adjacent, Cormat_gpu, N, r_thresh[k]);
					s_thresh[k] = 100.0 * Clength / M1 / 2;
				}
				else {
					Clength =  (long long) 2*M1*s_thresh[k]/100.0;
					Clength += Clength%2;
					r_thresh[k] = FormFullGraph_s(adjacent, Cormat_gpu, N, Clength/2);
				}
				time = clock() - time;
				cout<<"time for forming full graph (ms):\t" <<time<<endl;
				
				//if (k == 0)
				C = new int [Clength];

				//form csr graph
				time = clock();
				FormCSRGraph(R, C, adjacent, N);
				time = clock() - time;
						cout<<"time for forming csr graph (ms):\t" <<time<<endl;

				// write files
				sprintf_s(sparsity, "_spa%.3f%%_cor%.3f", s_thresh[k],r_thresh[k]);
				Outfilename = OutCor;
				Outfilename.append(sparsity).append(".csr");
				cout<<"generating "<<Outfilename.c_str()<< "..."<<endl;
				fout.open(Outfilename.c_str(), ios::binary | ios::out);
				fout.write((char*)&Rlength, sizeof(int));
				fout.write((char*)R, sizeof(int) * Rlength);
				fout.write((char*)&Clength, sizeof(int));
				fout.write((char*)C, sizeof(int) * Clength);
				fout.close();

				delete []C;
			}

			/*
			long long hisogram[21] = {0};
			real__t avg=0;
			int illegal=0;
			for (long long z=0; z<M1; z++){
				if (Cormat_gpu[z]<1000 || Cormat_gpu[z]>-1000)
				{
					hisogram[(int) ((Cormat_gpu[z]+1.05)*10)]++;
					avg += Cormat_gpu[z];
				}
				else illegal++;
			}
			avg = avg/M1;
			ofstream fo;
			fo.open("E:\\xumo\\4Dnii\\hisogramavg.txt");

			fo <<"\n=====================================\n";
			for (int num=0; num<21;num++){
				fo<<(100.0*(double) hisogram[num])/M1<<endl;								
			} 
			fo<<"\n=====================================\n";
			fo<<endl<<avg;
			fo<<endl<<100.0*((double) illegal)/M1;
			fo<<"\n=====================================\n";
			fo.close();
			*/

			
			delete []adjacent;
			delete []Cormat_gpu;
		}

		delete []R;
		delete []mask;
		delete []r_thresh;
		delete []filename;
		total_time = clock() - total_time;
		cout<<"total elapsed time: "<<1.0*total_time/1000<<" s."<<endl;
		cout<<"==========================================================="<<endl;
		return 0;
}


long long  FormFullGraph(bool * adjacent, real__t * Cormat, int N, real__t threshold)
{
	long long index = 0;
	long long nonZeroCount = 0;
	memset(adjacent, 0, sizeof(bool) * (long long)N * N);
	for(int i = 0; i < N; i++)
		for(int j = i+1; j < N; j++)
		{   if (Cormat[index] > 1 || Cormat[index] < -1)
			{
				cout<<"FormFullGraph:illegal correlation coefficient\n";
				printf("position: (%d , %d)\n",i, j);
				//cin >> j;
			}
			if (Cormat[index] > threshold)
			{
				nonZeroCount += (adjacent[(long long)i * N + j] = adjacent[(long long)j * N + i] = true);
			}
			index ++;
		}
	return nonZeroCount * 2;
}

long long partition(real__t *A, long long  m,long long  p){
        long long i;
		real__t	tem;
		real__t v;
        v=A[m];i=m+1;
		//cout<<v<<endl;
        while(1){
			while(A[i]>v && i<p)
				i++;
            while(A[p]<=v && i<=p)
				p--;
			if(i<p){
                tem=A[i];A[i]=A[p];A[p]=tem;
            }else break;
        }
        A[m]=A[p];A[p]=v;return (p);
}

void select(real__t *A,long long n,long long k){
	long long j,m,r;
	m=0;r=n-1;
	cout<<"k = "<<k<<endl;
	while(1){		
		j=r;
		//cout<<"m = "<<m<<" ; r = "<< r<< endl;
		//cout<<"A[m] = "<<A[m]<<" ; A[r] = "<< A[r]<< endl;
		//clock_t time = clock();
		j=partition(A, m,j);
		//time = clock() - time;
		//cout<<"partition time : "<<time<<"; j ="<<j<<endl;
		//cout<<"j = "<<j<<endl;
		if(k-1==j)break;
		else if(k-1<j) r=j-1;
		else m=j+1;
	}
}


real__t  FormFullGraph_s(bool * adjacent, real__t * Cormat, int N, long long threshold)
{
	

	long long M1 = (N-1);
	M1 *= N;
	M1 /= 2;

	//long long M = ((M1-2+blocksize)/blocksize)*blocksize+1;
	
	//  clock_t time = clock();
	//  cout<<"max_correlation : "<<Cormat[find_max(Cormat,M1)]<<endl;
	//  time = clock()-time;
	//  cout<<"find max time : "<<time<<endl;

	real__t * Cormat_copy = new real__t[M1];
	//real__t * Cormat_copy = new real__t[M];
	memcpy (Cormat_copy, Cormat, (long long) M1*sizeof(real__t));
	//memcpy (Cormat_copy, Cormat, (long long) M1*sizeof(real__t));
	//memset (Cormat_copy+M1, 0, (long long) (M-M1)*sizeof(real__t));

	/*for (long long i=0 ; i < M1; i++)
	{
		if (Cormat[i] >0.94) cout<< Cormat_copy[i]<<endl;
		if (Cormat_copy[i] != Cormat[i]) cout<<"i = "<<i<<endl;
	}*/
	//select_GPU(Cormat_copy, M1, threshold);
	select(Cormat_copy, M1, threshold);
	real__t r_threshold = Cormat_copy[threshold-1];
	cout <<"corresponding r_threshold = "<<r_threshold<<endl;
	delete []Cormat_copy;
	
	//long long n_edge = FormFullGraph(adjacent,  Cormat, N, r_threshold);
	long long index = 0;
	long long nonZeroCount = 0;
	stack<int> s_i;
	stack<int> s_j;
	int i_for_del;
	int j_for_del;

	memset(adjacent, 0,  sizeof(bool) *(long long) N * N);
	for(int i = 0; i < N; i++)
	{
		//if (nonZeroCount >= threshold)
		//		break;
		for(int j = i+1; j < N; j++)
		{   if (Cormat[index] > 1 || Cormat[index] < -1)
			{
		        cout<<"FormFullGraph:illegal correlation coefficient\n";
				printf("position: (%d , %d)\n",i, j);
				//cin >> j;
				//exit(1);
			}
			
		
			if (Cormat[index] >= r_threshold )
				nonZeroCount += (adjacent[(long long)i * N + j] = adjacent[(long long)j * N + i] = true);
			if (Cormat[index] == r_threshold)
			{	s_i.push(i);s_j.push(j);	}
			
			if (nonZeroCount > threshold)
			{					
				if (s_i.empty()) cout<<"nonZeroCount error!";
				adjacent[(long long)(s_i.top()) * N + s_j.top()] = false;
				adjacent[(long long)(s_j.top()) * N + s_i.top()] = false;
				s_i.pop();   s_j.pop();  nonZeroCount--;
			}

			index ++;
		}
	}

	if (nonZeroCount != threshold) 
	{
		cout<<"expected edge number = "<<threshold<<"\nreal edge number = "<<nonZeroCount<<endl; 
	    cin >> threshold;
		exit(1);
	}
	return  r_threshold;
	
}

void post_block (real__t * Cormat, real__t * Cormat_blocked, int N, int block_size,bool fishtran_flag)
{
/*	int i = 0;
	int block_cnt = (N + block_size - 1) / block_size;
	real__t * src, * dst;
	src = Cormat_blocked;
	dst = Cormat;
	for (i = 0; i < (block_cnt - 1) * block_size ; i++)
	{
		int offset = i % block_size + 1;
		memcpy(dst, src + offset, sizeof(real__t) * (block_size - offset));
		dst += block_size - offset;
		src += block_size * block_size;
		for (int j = i / block_size + 1; j < block_cnt - 1; j++)
		{
			memcpy(dst, src, sizeof(real__t) * block_size);
			dst += block_size;
			src += block_size * block_size;
		}
		memcpy(dst, src, sizeof(real__t) * (N % block_size));
		dst += N % block_size;
		src += block_size - ((i+1) % block_size != 0) * block_size * block_size * (block_cnt - 1 - i / block_size);
	}
	for (; i < N; i++)
	{
		int offset = i % block_size + 1;
		memcpy(dst, src + offset, sizeof(real__t) * (N % block_size - offset));
		dst += (N % block_size - offset);
		src += block_size;
	}
	*/

	int block_cnt = (N + block_size - 1)/ block_size;
	long long index = 0;
	int nonzeros = 0;
	if (fishtran_flag)
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
			{
				int block_row = i / block_size;
				int block_col = j / block_size;
				int block_i   = i % block_size;
				int block_j   = j % block_size;
				long long offset    =block_row * (2*block_cnt - block_row + 1) / 2 + block_col - block_row;
				offset *= block_size * block_size;
				if (Cormat_blocked[(long long)block_i*block_size + block_j + offset]>=-1 && Cormat_blocked[(long long)block_i*block_size + block_j + offset]<=1)
						Cormat[index] += fisher_trans(Cormat_blocked[(long long)block_i*block_size + block_j + offset]);
				index ++;
			}
	else
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
			{
				int block_row = i / block_size;
				int block_col = j / block_size;
				int block_i   = i % block_size;
				int block_j   = j % block_size;
				long long offset    =block_row * (2*block_cnt - block_row + 1) / 2 + block_col - block_row;
				offset *= block_size * block_size;
				if (Cormat_blocked[(long long)block_i*block_size + block_j + offset]>=-1 && Cormat_blocked[(long long)block_i*block_size + block_j + offset]<=1)
						Cormat[index] += Cormat_blocked[(long long)block_i*block_size + block_j + offset];
				index ++;
			}
}
void FormCSRGraph(int * R, int * C, bool * adjacent, int N)
{

	int count = 0;
	long long Cindex = 0;
	R[0] = 0;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			if(adjacent[(long long)i*N+j])
			{
				count ++;
				C[Cindex++] = j;
			}
		}
		R[i+1] = count;
	} 
}


void CorMat_cpu(real__t * Cormat, real__t * BOLD, int N, int L)
{	
	// transposing the BOLD signal
	real__t * BOLD_t = new real__t [L * N];
	memset(BOLD_t, 0, sizeof(real__t) * L * N);
	for (int i = 0; i < L; i ++)
		for (int j = 0; j < N; j++)
		{
			BOLD_t[j * L + i] = BOLD[i * N + j];
		}

	// Normalize
	for (int i = 0; i < N; i++)
	{
		real__t * row = BOLD_t + i * L;
		double sum1 = 0, sum2 = 0;
		for (int l = 0; l < L; l++)
		{
			sum1 += row[l];
		}
		sum1 /= L;
		for (int l = 0; l < L; l++)
		{
			sum2 += (row[l] - sum1) * (row[l] - sum1);
		}
		sum2 = sqrt(sum2);
		for (int l = 0; l < L; l++)
		{
			row[l] = (row[l] - sum1) / sum2;;
		}
	}


	// GEMM
	double sum3;
	for (int k = 0, i = 0; i < N; i++)
		for (int j = i+1; j < N; j++)
		{
			sum3 = 0;
			for (int l = 0; l < L; l++)
			{
				sum3 += BOLD_t[i*L+l] * BOLD_t[j*L+l];
			}
			Cormat[k++] = sum3;
		}
	delete []BOLD_t;
}
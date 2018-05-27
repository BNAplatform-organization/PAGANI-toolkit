#include <stdlib.h>
#include <memory.h>
#include <string>
#include <ctime>
#include <cmath>
#include <iomanip>
#include "dirent.h"
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <algorithm>
#include "help_func.cuh"
#include "data_type.h"
#include "shlwapi.h"       
//used for creating folder
#pragma comment(lib,"shlwapi.lib")
using namespace std;

#define argnum 7
typedef float real__t;
typedef unsigned int uint__t;

const int HdrLen = 352;
double ProbCut = 0.5; 
const int blocksize = 1024*1024*48;

bool rs_flag = false;
bool cormat_flag = false;

void CorMat_cpu(real__t * Cormat, real__t * BOLD, int N, int L);

real__t CorMat_spa2rth_blocking(string OutCor, real__t * BOLD_t, int N,  int L, int  Batch_size,real__t *s_thresh, real__t* result, int num_spa, const int blocknum);
int CorMat_gpu_blocking(string OutCor, real__t * BOLD_t, const int &N, const int &N0, const int &Num_Blocks, const int &L, const int &Batch_size, V_type *r_thresh, const int &NumS, const int blocknum);

real__t CorMat_spa2rth(string OutCor, real__t * BOLD_t, int N,  int L, int  Batch_size,real__t *s_thresh, real__t* result, int num_spa);
int CorMat_gpu(string OutCor, real__t * BOLD_t, const int &N, const int &N0, const int &Num_Blocks, const int &L, const int &Batch_size, V_type *r_thresh, const int &NumS);

void SpearmanAssignmentAndNormalization(real__t* BOLD_t, int N, int L);
void SpearmanAssignmentBroadcast(real__t* BOLD_t, int N, int L, bool blocking);
void output(real__t* BOLD_t, int N, int L,string str, bool);

//long long FormFullGraph(bool * adjacent, real__t * Cormat, int N, real__t threshold);
//real__t FormFullGraph_s(bool * adjacent, real__t * Cormat, int N, long long threshold);
//void post_block (real__t * Cormat, real__t * Cormat_blocked, int dim, int block_size,bool fishtran_flag);
//void FormCSRGraph(int * R, int * C, real__t *V, bool * adjacent, int N , real__t *Cormat);
//long long find_max(real__t *Cormat, long long M1);
//real__t select_GPU(real__t *Cormat, long long M1, long long k);

real__t fishertrans(real__t r)
{
	real__t z;
	if (r==1) r -= 1e-6; 
	z = 0.5*log((1+r)/(1-r));
	return z;	
}

real__t inv_fishertrans(real__t z)
{
	real__t r;
	r = exp(2*z);
	r = (r-1)/(r+1);
	return r;	
}

enum string_code {
    spearman,
    pearson,
    spearmanAndpearson
};
string_code hashit (std::string const& inString) {
	if (inString == "s") return spearman;
	if (inString == "p") return pearson;
	if (inString == "ps" || inString == "ps") return spearmanAndpearson;
}

/************** specialized for assignment rank on CPU *******************/
//struct myPair{
//real__t value ;
//int addr;
//} *Pair;
//bool functorByValue (myPair i,myPair j) 
//{ return (i.value<j.value); }
/*************************************************************************/

int main(int argc, char * argv[])
{
	clock_t total_time = clock();
	if (argc < 7) 
	{
		cerr<<"Input format: .\\CorMat.exe  Dir_for_BOLD Path_for_mask threshold_for_mask(0~1) to_average(yf/yn/bf/bn/n) threshold_type(r/s) threshold_for_correletaion_coefficient(s)(0~1)\n"
			<<"For example: .\\CorMat.exe  d:\\BOLD d:\\MASK\\mask.nii 0.5 y n r 0.2 0.25 0.3\n"<<endl;
		exit(1);	
	}
	int L, N = 0, i = 0, j = 0, k = 0, l = 0, total_size;
	clock_t time;

	/*************************************************************************************************/
	/*                                 Search fMRI (nifti) files                                     */
	/*************************************************************************************************/
	DIR *dp;
	struct dirent *dirp;
	if (NULL == (dp = opendir(argv[1])))
	{
		printf("can't open %s", argv[1]);
		system("pause");
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
			if (filenametmp.compare("mask.nii")!=0&&filenametmp.compare("grey_mask.nii")!=0&&filenametmp.compare("1mm_grey.nii")!=0&&filenametmp.compare("1mm_grey_mask.nii")!=0)
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
			if (filenametmp.compare("mask.nii")!=0&&filenametmp.compare("grey_mask.nii")!=0)
				filename[i++] = filenametmp;
		}
	}

	/*************************************************************************************************/
	/*                                 threshold strategy setup                                      */
	/*************************************************************************************************/
	int NumS = argc - argnum;
	V_type * r_thresh = new V_type [NumS];
	float * s_thresh = new float [NumS];
	if (argv[5][0] == 'r' || argv[5][0] == 'R' )
		for (i = 0; i < NumS; i++)
			r_thresh[i] = (real__t)atof(argv[argnum+i]);
	else if (argv[5][0] == 's' || argv[5][0] == 'S' )
	{
		rs_flag = true;
		memset(r_thresh, 0, sizeof(real__t)*NumS);
		for (i = 0; i < NumS; i++)
			s_thresh[i] = (real__t)atof(argv[argnum+i]);
	}
	else {
		cout << "threshold type error! \nr for correlation threshold that is sole currently.\n";
		system("pause");
		exit(1);
	}
	
	/*************************************************************************************************/
	/*								Read Mask file, obtain mask_idx                                  */
	/*************************************************************************************************/
	real__t ProbCut = (real__t)atof(argv[3]);
	string mask_file = string(argv[2]); 
	ifstream fin(mask_file.c_str(), ios_base::binary);
	if (!fin.good())
	{	cout<<"Can't open\t"<<mask_file.c_str()<<endl;	return 0; }
	short hdr[HdrLen / 2];
	fin.read((char*)hdr, HdrLen);
	cout<<"mask datatype : "<< hdr[35]<<"  "<<hdr[36]<<endl;
	char mask_dt = hdr[35];
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
		system("pause");
	    return -1;
	}
	fin.close();
	
	// Count the number of the valid voxels, and obtain mask_idx vector	
	vector<int> mask_idx;
	for (k = 0; k < total_size; k++)
		if (mask[k] >= ProbCut)
			mask_idx.push_back(k);
	N = mask_idx.size();
	cout<<"Data size: "<<hdr[21] <<"x"<<hdr[22]<<"x"<<hdr[23]  <<", Grey voxel count: "<<N<<"."<<endl;
	//cant figure out logic. so reserve mask
	
	/*************************************************************************************************/
	/*								Read fMRI files and process                                      */
	/*************************************************************************************************/
	if (argv[1][strlen(argv[1]) - 1] == '\\')
		argv[1][strlen(argv[1]) - 1] = 0;
	string str = string(argv[1]).append("\\").append("unweighted");
	if (!PathIsDirectory(str.c_str()))
	{
		::CreateDirectory(str.c_str(), NULL);
	}

//#ifdef figure
//	for (int fi = 0; fi < 1; fi++)//FileNumber; fi++)
//	{
//		
//		string a = string(argv[1]).append("\\").append(filename[fi]);
//		cout<<"\ncalculating correlation for "<<a.c_str()<<" ..."<<endl;
//		ifstream fin(a.c_str(), ios_base::binary);
//
//		// Get the length of the time sequence 
//		if (!fin.good())
//		{	cout<<"Can't open\t"<<a.c_str()<<endl;	system("pause"); return 0;  }
//		fin.read((char*)hdr, HdrLen);
//		L = hdr[24];
//		L = 128;
//		streampos sp = fin.tellg();
//		for(N = 200000;N < 1200000;N = N + 200000)
//		{
//			cout<<"calculate case for N="<<N<<" .."<<endl;
//			//L = 120;
//			fin.seekg(sp);
//			if (L==1)
//			{
//				cout<<a.c_str()<<"is not a 4D nifti file. Continue to the next file."<<endl;
//				fin.close();
//				continue;
//			}
//
//			/*********************   Get BOLD signal data   ********************/
//			int Batch_size = 1024 * 3;				
//			Batch_size *= 2; //6 fail
//			//Batch_size = 1024;
//			const int Num_Blocks = (N + Batch_size - 1) / Batch_size;
//			uint__t N0 = Num_Blocks * Batch_size;
//			real__t * BOLD_t = new real__t [L * N0];
//			memset(BOLD_t, 0, sizeof(real__t) * L * N0);
//			//real__t * BOLD = new real__t [L * N];
//			if (hdr[36] == 64) // double
//			{
//				double * InData = new double [total_size];
//				for (l = 0; l < L; l++) 
//				{
//					fin.read((char*)InData, sizeof(double) * total_size);
//					// Get the BOLD signal for all the valid voxels
//					for (int i = 0; i < N; i++)
//					{
//						BOLD_t[i*L+l] = InData[mask_idx[i]];
//					}
//				}
//				cout<<"BOLD length: "<<L<<", Data type: double. "<<N<<endl;
//				delete []InData;
//				//fin.close();
//			}
//			else if (hdr[36] == 32)	   //float
//			{
//#ifndef figure3
//				try{
//					real__t *InData = new float [L * total_size];
//					fin.read((char*)InData, sizeof(float) * L * total_size); 
//					//L = 120;
//					//fin.close();
//					// Get the BOLD signal for all the valid voxels
//					for (int i = -1, k = 0; k < total_size; k++)
//						if (mask[k] >= ProbCut)
//							for (i++, l = 0; l < L; l++)
//							{
//								BOLD_t[l*N+i] = InData[l*total_size+k];
//							}
//					delete []InData;
//				}
//				catch(...){
//					cout<<"total size: "<<total_size<<"*"<<L<<endl;
//					cout<<"The file will be masked in parts."<<endl;
//					//real__t * BOLD_com = new real__t [L * N];
//					uint__t piece = 100;
//					uint__t tile = total_size*piece;
//					int sheet = L / piece + 1;
//					real__t *InData = new float [tile];//[L * total_size];
//					for (uint__t z = 0; z < sheet; z++)
//					{
//						cout<<"Complete the "<<z<<"th part."<<endl;
//						uint__t bound = piece*(z+1)< L?piece*(z+1):L;
//						fin.read((char*)InData, sizeof(float) * tile);
//						for (int i = -1, k = 0; k < total_size; k++)
//							if (mask[k] >= ProbCut)
//								for (i++, l = 0+piece*z; l <bound; l++)
//								{
//									BOLD_t[l*N+i] = InData[(l-piece*z)*total_size+k];
//								}
//					}
//					//fin.close();
//					delete []InData;
//				}  
//#else
//				string ssss = "E:\\project\\CuCorMat-spearman-pearson-pan\\Data.datt";
//				ifstream fdata(ssss.c_str(), ios_base::binary);
//				if (!fdata.good())
//				{	cout<<"Can't open\t"<<ssss.c_str()<<endl;	system("pause"); }
//				fdata.read((char*)BOLD_t, sizeof(float) * N * L);
//#endif
//				cout<<"BOLD length: "<<L<<", Data type: float. "<<N<<endl;
//			}
//			else
//			{
//				cerr<<"Error: Data type is neither float nor double."<<endl;
//				fin.close();
//				continue;
//			}	
//		
//			/*****************   CPU programme for comparison   *****************/	
//			/*
//			Pair = new myPair[L];
//		
//			real__t *temp = new real__t[L];
//
//			for (int i = 0; i < N; i++)
//			{
//				real__t * row = BOLD_t + i * L;
//				for (int j = 0; j < L; j++)
//				{
//					Pair[j].addr = j;
//					Pair[j].value = row[j];
//				}
//				sort(Pair,Pair + L, functorByValue);
//			   未完待续！！
//
//
//
//			}
//			delete[] Pair;*/
//
//			/******************  whether or not start blocked transmission  ***************/
//			bool blocking = false;
//			cudaSetDevice(0);
//			cudaDeviceProp deviceProp;
//			cudaGetDeviceProperties(&deviceProp, 0);
//			double requiredMem = N0 * L * sizeof(real__t) ;
//			blocking =  requiredMem > (deviceProp.totalGlobalMem * 0.750f) ;
//#ifdef myDebug
//			cout<<"blocking:"<<blocking<<endl;  
//			cout<<"deviceProp.totalGlobalMem:"<<deviceProp.totalGlobalMem<<endl;
//#endif
//			/********************   correlation type    ************************/
//			//execution order: spearman\pearson and kendall?
//			unsigned int state = 0x00;
//			switch ( hashit ( string(argv[6]) ) )
//			{
//			case spearman: state = 0x10; break;
//			case pearson:  state = 0x01; break;
//			case spearmanAndpearson:state = 0x11; break;
//			default:state = 0x01;break;
//			}
//			for (int e = 0; e < 2; e++)
//			{
//				if (e == 0 && state == 0x01 )
//					continue;       
//				if (e == 1 && state == 0x10 )
//					continue;       
//			
//				/**********************   Spearman preprocessing   ******************/
//				//an auto judgement function for whether or not blocking ought to be writen in near future. difficulty: dont,t sure how to simply get gpu mem size.
//				if ( e==0 )
//				{   //need to modify resultant file name!!!!!
//					cout<<"computing spearman coefficient:"<<endl;
//					SpearmanAssignmentBroadcast(BOLD_t, N, L, blocking);
//					//SpearmanAssignmentAndNormalization(BOLD_t, N, L);
//				}else
//				{
//					cout<<"Computing pearson coefficient:"<<endl;
//				}
//
//				/**********************   Normalize BOLD_t   **********************/	
//				for (int i = 0; i < N; i++)
//				{
//					real__t * row = BOLD_t + i * L;
//					double sum1 = 0, sum2 = 0;
//					for (int l = 0; l < L; l++)
//					{
//						sum1 += row[l];
//					}
//					sum1 /= L; //average mean
//					for (int l = 0; l < L; l++)
//					{
//						sum2 += (row[l] - sum1) * (row[l] - sum1);
//					}
//					sum2 = sqrt(sum2); //standardized variance
//					for (int l = 0; l < L; l++)
//					{
//						row[l] = (row[l] - sum1) / sum2;;
//					}
//				}                     //normalized
//		
//				/**********************   Cormat Computation   **********************/
//				char Graph_size[20];
//		
//				a=string(argv[1]).append("\\weightedSpearman\\");
//				a.append(filename[fi]);
//				string OutCor;
//				if ( e == 0)
//					OutCor = a.substr(0, a.find_last_of('.')).append("_").append(string(itoa(N, Graph_size, 10))).append("_spearman");//cascade N to file name;
//				else
//					OutCor = a.substr(0, a.find_last_of('.')).append("_").append(string(itoa(N, Graph_size, 10))).append("_pearson");//cascade N to file name;
//				//calculate r_thresh for s_thresh
//				if (blocking == false)
//				{
//					if (rs_flag)
//					CorMat_spa2rth(OutCor, BOLD_t, N, L, Batch_size,s_thresh,  r_thresh,NumS);
//					//cout<<"*r_thresh:"<<*r_thresh<<endl;
//					// sort r_thresh (increase)
//					sort(r_thresh, r_thresh+NumS); //r_th_min == r_thresh[0]
//					/*************  Construct CSR network for each r_thresh  ************/
//					CorMat_gpu(OutCor, BOLD_t, N, N0, Num_Blocks, L, Batch_size, r_thresh, NumS); 	
//				}else
//				{
//					const int transferblocknum = 102;
//					if (rs_flag)
//						CorMat_spa2rth_blocking(OutCor, BOLD_t, N, L, Batch_size,s_thresh,  r_thresh,NumS, transferblocknum);
//						//cout<<"*r_thresh:"<<*r_thresh<<endl;
//						// sort r_thresh (increase)
//					sort(r_thresh, r_thresh+NumS); //r_th_min == r_thresh[0]
//					/*************  Construct CSR network for each r_thresh  ************/
//					CorMat_gpu_blocking(OutCor, BOLD_t, N, N0, Num_Blocks, L, Batch_size, r_thresh, NumS, transferblocknum); 
//				}
////#ifdef multiblock
////				const int transferblocknum = 152;
////				if (rs_flag)
////					CorMat_spa2rth_blocking(OutCor, BOLD_t, N, L, Batch_size,s_thresh,  r_thresh,NumS, transferblocknum);
////					//cout<<"*r_thresh:"<<*r_thresh<<endl;
////					// sort r_thresh (increase)
////				sort(r_thresh, r_thresh+NumS); //r_th_min == r_thresh[0]
////				/*************  Construct CSR network for each r_thresh  ************/
////				CorMat_gpu_blocking(OutCor, BOLD_t, N, N0, Num_Blocks, L, Batch_size, r_thresh, NumS, transferblocknum); 	
////#else
////				if (rs_flag)
////					CorMat_spa2rth(OutCor, BOLD_t, N, L, Batch_size,s_thresh,  r_thresh,NumS);
////					//cout<<"*r_thresh:"<<*r_thresh<<endl;
////					// sort r_thresh (increase)
////				sort(r_thresh, r_thresh+NumS); //r_th_min == r_thresh[0]
////				/*************  Construct CSR network for each r_thresh  ************/
////				CorMat_gpu(OutCor, BOLD_t, N, N0, Num_Blocks, L, Batch_size, r_thresh, NumS); 	
////#endif
//
//			}
//		
//			delete []BOLD_t;			
//		}
//		fin.close();
//			
//	}	
//
//#else
	for (int fi = 0; fi < FileNumber; fi++)
	{
		
		string a = string(argv[1]).append("\\").append(filename[fi]);
		cout<<"\ncalculating correlation for "<<a.c_str()<<" ..."<<endl;
		ifstream fin(a.c_str(), ios_base::binary);

		// Get the length of the time sequence 
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	system("pause"); return 0;  }
		fin.read((char*)hdr, HdrLen);
		L = hdr[24];
		//L = 120;
		if (L==1)
		{
			cout<<a.c_str()<<"is not a 4D nifti file. Continue to the next file."<<endl;
			fin.close();
			continue;
		}

		/*********************   Get BOLD signal data   ********************/
		int Batch_size = 1024 * 3;				
		Batch_size *= 2; //6 fail
		//Batch_size = 1024;
		const int Num_Blocks = (N + Batch_size - 1) / Batch_size;
		uint__t N0 = Num_Blocks * Batch_size;
		real__t * BOLD_t = new real__t [L * N0];
		memset(BOLD_t, 0, sizeof(real__t) * L * N0);
		//real__t * BOLD = new real__t [L * N];
		if (hdr[36] == 64) // double
		{
			double * InData = new double [total_size];
			for (l = 0; l < L; l++) 
			{
				fin.read((char*)InData, sizeof(double) * total_size);
				// Get the BOLD signal for all the valid voxels
				for (int i = 0; i < N; i++)
				{
					BOLD_t[i*L+l] = InData[mask_idx[i]];
				}
			}
			cout<<"BOLD length: "<<L<<", Data type: double. "<<N<<endl;
			delete []InData;
			fin.close();
		}
		else if (hdr[36] == 32)	   //float
		{

			real__t *InData = new float [total_size];
			for (l = 0; l < L; l++) 
			{
				fin.read((char*)InData, sizeof(float) * total_size);
				/*      Read BOLD signal for all valid voxels     */
				for (int i = 0; i < N; i++)
				{
					BOLD_t[i*L+l] = InData[mask_idx[i]];
				}
			}
			cout<<"BOLD length: "<<L<<", Data type: float. "<<N<<endl;
			delete []InData;
			fin.close();

		}
		else
		{
			cerr<<"Error: Data type is neither float nor double."<<endl;
			fin.close();
			continue;
		}	
		
		/*****************   CPU programme for comparison   *****************/	
		/*
		Pair = new myPair[L];
		
		real__t *temp = new real__t[L];

		for (int i = 0; i < N; i++)
		{
			real__t * row = BOLD_t + i * L;
			for (int j = 0; j < L; j++)
			{
				Pair[j].addr = j;
				Pair[j].value = row[j];
			}
			sort(Pair,Pair + L, functorByValue);
           未完待续！！



		}
		delete[] Pair;*/

		/******************  whether or not start blocked transmission  ***************/
		bool blocking = false;
		cudaSetDevice(0);
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, 0);
		double requiredMem = N0 * L * sizeof(real__t) ;
		blocking =  requiredMem > (deviceProp.totalGlobalMem * 0.750f) ;
		/********************   correlation type    ************************/
		//execution order: spearman\pearson and kendall?
		unsigned int state = 0x00;
		switch ( hashit ( string(argv[6]) ) )
		{
		case spearman: state = 0x10; break;
		case pearson:  state = 0x01; break;
		case spearmanAndpearson:state = 0x11; break;
		default:state = 0x01;break;
		}
		for (int e = 0; e < 2; e++)
		{
			if (e == 0 && state == 0x01 )
				continue;       
			if (e == 1 && state == 0x10 )
				continue;       
			
		
			/************************   debug mode    **************************/
#ifdef myDebug
			string plus = "\\";
			string partialStr = str + plus + filename[fi];
			output(BOLD_t,N,L,partialStr,false);
#endif
			/**********************   Spearman preprocessing   ******************/
			//an auto judgement function for whether or not blocking ought to be writen in near future. difficulty: dont,t sure how to simply get gpu mem size.
			if ( e==0 )
			{   //need to modify resultant file name!!!!!
				cout<<"computing spearman coefficient:"<<endl;
				SpearmanAssignmentBroadcast(BOLD_t, N, L,blocking);
			}else
			{
				cout<<"Computing pearson coefficient:"<<endl;
			}
			

			/************************   debug mode    **************************/
#ifdef myDebug
			string Plus = "\\";
			string partialStrRanked = str + Plus + filename[fi];
			output(BOLD_t,N,L,partialStrRanked,false);
#endif

			/**********************   Normalize BOLD_t   **********************/	
#ifdef myDebug
			clock_t normTime = clock();
#endif
			for (int i = 0; i < N; i++)
			{
				real__t * row = BOLD_t + i * L;
				double sum1 = 0, sum2 = 0;
				for (int l = 0; l < L; l++)
				{
					sum1 += row[l];
				}
				sum1 /= L; //average mean
				for (int l = 0; l < L; l++)
				{
					sum2 += (row[l] - sum1) * (row[l] - sum1);
				}
				sum2 = sqrt(sum2); //standardized variance
				for (int l = 0; l < L; l++)
				{
					row[l] = (row[l] - sum1) / sum2;;
				}
			}                     //normalized
#ifdef myDebug
			normTime = clock() - normTime;
			cout<<"elapsed time for cpu normalization:"<<normTime/1000.0<<"s."<<endl;
#endif
			/*for (int i = 0; i < N; i++)
			{
				real__t * row = BOLD_t + i * L;
				for (int l = 0; l < L; l++)
				{
					sum1 += row[l];
				}
			}*/
		
			/**********************   Cormat Computation   **********************/
			char Graph_size[20];
		
			a=string(argv[1]).append("\\unweighted\\");
			a.append(filename[fi]);
			string OutCor;
			if ( e == 0)
				OutCor = a.substr(0, a.find_last_of('.')).append("_").append(string(itoa(N, Graph_size, 10))).append("_spearman");//cascade N to file name;
			else
				OutCor = a.substr(0, a.find_last_of('.')).append("_").append(string(itoa(N, Graph_size, 10))).append("_pearson");//cascade N to file name;
			//calculate r_thresh for s_thresh

			if (blocking == false)
			{
				if (rs_flag)
				CorMat_spa2rth(OutCor, BOLD_t, N, L, Batch_size,s_thresh,  r_thresh,NumS);
				//cout<<"*r_thresh:"<<*r_thresh<<endl;
				// sort r_thresh (increase)
				sort(r_thresh, r_thresh+NumS); //r_th_min == r_thresh[0]
				/*************  Construct CSR network for each r_thresh  ************/
				CorMat_gpu(OutCor, BOLD_t, N, N0, Num_Blocks, L, Batch_size, r_thresh, NumS); 	
			}else
			{
				const int transferblocknum = 102;
				if (rs_flag)
					CorMat_spa2rth_blocking(OutCor, BOLD_t, N, L, Batch_size,s_thresh,  r_thresh, NumS, transferblocknum);
					//cout<<"*r_thresh:"<<*r_thresh<<endl;
					// sort r_thresh (increase)
				sort(r_thresh, r_thresh+NumS); //r_th_min == r_thresh[0]
				/*************  Construct CSR network for each r_thresh  ************/
				CorMat_gpu_blocking(OutCor, BOLD_t, N, N0, Num_Blocks, L, Batch_size, r_thresh, NumS, transferblocknum); 
			}

		}
		
		delete []BOLD_t;			
		
			
	}	

	delete []mask;
	mask_idx.clear();
	mask_idx.swap(vector<int>()) ; 
	delete []r_thresh;
	delete []filename;
	total_time = clock() - total_time;
	cout<<"total elapsed time: "<<1.0*total_time/1000<<" s."<<endl;
	cout<<"==========================================================="<<endl;
	
	return 0;	
}


/*long long  FormFullGraph(bool * adjacent, real__t * Cormat, int N, real__t threshold)
{
	long long index = 0;
	long long nonZeroCount = 0;
	memset(adjacent, 0, sizeof(bool) * (long long)N * N);
	for(int i = 0; i < N; i++)
		for(int j = i+1; j < N; j++)
		{   if (Cormat[index] > 1)
			{
		        cout<< Cormat[index]<<endl;
				cout<<"FormFullGraph:illegal correlation coefficient\n";				
				exit(1);
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
}*/


//real__t  FormFullGraph_s(bool * adjacent, real__t * Cormat, int N, long long threshold)
//{
//	
//
//	long long M1 = (N-1);
//	M1 *= N;
//	M1 /= 2;
//
//	//long long M = ((M1-2+blocksize)/blocksize)*blocksize+1;
//	
//	//  clock_t time = clock();
//	//  cout<<"max_correlation : "<<Cormat[find_max(Cormat,M1)]<<endl;
//	//  time = clock()-time;
//	//  cout<<"find max time : "<<time<<endl;
//
//	real__t * Cormat_copy = new real__t[M1];
//	//real__t * Cormat_copy = new real__t[M];
//	memcpy (Cormat_copy, Cormat, (long long) M1*sizeof(real__t));
//	//memcpy (Cormat_copy, Cormat, (long long) M1*sizeof(real__t));
//	//memset (Cormat_copy+M1, 0, (long long) (M-M1)*sizeof(real__t));
//
//	/*for (long long i=0 ; i < M1; i++)
//	{
//		if (Cormat[i] >0.94) cout<< Cormat_copy[i]<<endl;
//		if (Cormat_copy[i] != Cormat[i]) cout<<"i = "<<i<<endl;
//	}*/
//	//select_GPU(Cormat_copy, M1, threshold);
//	select(Cormat_copy, M1, threshold);
//	real__t r_threshold = Cormat_copy[threshold-1];
//	cout <<"corresponding r_threshold = "<<r_threshold<<endl;
//	delete []Cormat_copy;
//	
//	//long long n_edge = FormFullGraph(adjacent,  Cormat, N, r_threshold);
//	long long index = 0;
//	long long nonZeroCount = 0;
//	stack<int> s_i;
//	stack<int> s_j;
//	int i_for_del;
//	int j_for_del;
//
//	memset(adjacent, 0,  sizeof(bool) *(long long) N * N);
//	for(int i = 0; i < N; i++)
//	{
//		//if (nonZeroCount >= threshold)
//		//		break;
//		for(int j = i+1; j < N; j++)
//		{   if (Cormat[index] > 1 || Cormat[index] < -1)
//			{
//		        cout<<"FormFullGraph:illegal correlation coefficient\n";
//				cin >> j;
//				exit(1);
//			}
//			
//		
//			if (Cormat[index] >= r_threshold )
//				nonZeroCount += (adjacent[(long long)i * N + j] = adjacent[(long long)j * N + i] = true);
//			if (Cormat[index] == r_threshold)
//			{	s_i.push(i);s_j.push(j);	}
//			
//			if (nonZeroCount > threshold)
//			{					
//				if (s_i.empty()) cout<<"nonZeroCount number error!";
//				adjacent[(long long)(s_i.top()) * N + s_j.top()] = false;	
//				adjacent[(long long)(s_j.top()) * N + s_i.top()] = false;
//				s_i.pop();   s_j.pop();  nonZeroCount--;
//			}
//
//			index ++;
//		}
//	}
//
//	if (nonZeroCount != threshold) 
//	{
//		cout<<"expected edge number = "<<threshold<<"\nreal edge number = "<<nonZeroCount<<endl; 
//	    cin >> threshold;
//		exit(1);
//	}
//	return  r_threshold;
//	
//}

void output(real__t* BOLD_t, int N, int L, string path, bool afterRank)
{
	ofstream fout;
	string Outfilename;
	if (afterRank==false)
		Outfilename = path.append("_BOLD_t.matrix");
	else
		Outfilename = path.append("_BOLD_t_AfterRank.matrix");
	fout.open(Outfilename.c_str(), ios::binary | ios::out);
	if (!fout)
	{
		cout<<"create outfile unsuccessfully. error code:  "<<GetLastError()<<endl;		exit(false);
	}	
	fout.write((const char*)&BOLD_t[0],sizeof(real__t) * L * N );
	fout.close();
	if (afterRank==false)
		cout<<"Complete BOLD_t matrix transfering."<<endl;
	else
		cout<<"Complete ranked BOLD_t matrix transfering."<<endl;
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

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <process.h>
#include <windows.h>

# include <iostream>
# include <fstream>
using namespace std;

typedef float real__t;
typedef unsigned int uint__t;

const int HdrLen = 352;
double ProbCut = 0.5; 
int numThread;


struct Corr_ARG
{
	int id;
	real__t * Cormat;
	real__t * Bold_t;
	real__t * moment2;
	int N;
	int L;
};


void CorMat_cpu(real__t * Cormat, real__t * BOLD, int N, int L);
void CorMat_cpu_MT(real__t * Cormat, real__t * BOLD, int N, int L);
int FormFullGraph(bool * adjacent, real__t * Cormat, int N, real__t threshold);
void FormCSRGraph(int * R, int * C, bool * adjacent, int N);
void Cross_term(void *input_arg);

int main(int argc, char * argv[])
{
	if (argc < 4) 
	{
		cerr<<"Input format: .\\CorMat.exe  niiFile_for_BOLD  niiFile_for_mask  threshold_for_mask  threshold_for_corre_coeff(multiple)\n"
			<<"For example: .\\CorMat.exe  X.BOLD.nii  mask.nii  0.5  0.25  0.3  0.4\n"<<endl;
		exit(1);	
	}
	int L, N = 0, i = 0, j = 0, k = 0, l = 0, total_size;
	clock_t time;

	// read input files and parameters
	real__t ProbCut = (real__t)atof(argv[3]);
	int NumS = argc - 4;
	real__t * S_thresh = new real__t [NumS];
	for (i = 0; i < NumS; i++)
		S_thresh[i] = (real__t)atof(argv[4+i]);

	char * input = argv[1];
	char * mask_file = argv[2];
	ifstream fin(mask_file, ios::binary);
	if (!fin.good())
	{	cout<<"Can't open\t"<<mask_file<<endl;	return 0;}
	short hdr[HdrLen / 2];
	fin.read((char*)hdr, HdrLen);
	total_size = hdr[21] * hdr[22] * hdr[23];	// Total number of voxels
	real__t * mask = new real__t [total_size];
	fin.read((char*)mask, sizeof(real__t) * total_size);
	fin.close();

	// Count the number of the valid voxels	
	for (k = 0; k < total_size; k++)
		N += (mask[k] >= ProbCut);				

	// Get the length of the time sequence 
	fin.open(input, ios::binary);
	if (!fin.good())
	{	cout<<"Can't open\t"<<input<<endl;	return 0;}
	fin.read((char*)hdr, HdrLen);
	L = hdr[24];
	real__t * InData = new real__t [L * total_size];
	fin.read((char*)InData, sizeof(real__t) * L * total_size);
	fin.close();
	real__t * BOLD = new real__t [L * N];

	// Get the BOLD signal for all the valid voxels
	for (i = -1, k = 0; k < total_size; k++)
		if (mask[k] >= ProbCut)
			for (i++, l = 0; l < L; l++)
			{
				BOLD[l*N+i] = InData[l*total_size+k];
			}
	cout<<"Network size\t\t"<<N<<endl;

	long long M1 = (N-1);
	M1 *= N;
	M1 /= 2;

	// calculate on CPU 
	SYSTEM_INFO siSysInfo;
	GetSystemInfo(&siSysInfo); 
	numThread = (int)siSysInfo.dwNumberOfProcessors; 

	real__t * Cormat_cpu = new real__t [M1];
	time = clock();
	CorMat_cpu_MT(Cormat_cpu, BOLD, N, L);
	time = clock() - time;
	cout<<"Multi thread time (ms):\t"<<time<<endl;

/*
	// compute on single core for testing
	real__t * Cormat_cpu_st = new real__t [M1];
	time = clock();
	CorMat_cpu(Cormat_cpu_st, BOLD, N, L);
	time = clock() - time;
	cout<<"Single thread time (ms):\t"<<time<<endl;
	
	for (long long i = 0; i < M1; i++)
		if (fabs(Cormat_cpu[i] - Cormat_cpu_st[i]) > 1e-6)
		{
			cout<<"coule"<<endl;
			getchar();
		}
		cout<<"niule"<<endl;
*/

	// swap the largest threshold to the beginning
	int min_idx = 0;
	for (i = 0; i < NumS; i++)
		if (S_thresh[i] < S_thresh[min_idx])
			min_idx = i;
	real__t temp = S_thresh[0];
	S_thresh[0] = S_thresh[min_idx];
	S_thresh[min_idx] = temp;

	// CSR format
	int Rlength = N + 1;
	int Clength;
	int * R = new int[Rlength];
	int * C;
	
	cout<<"generating CSR files..."<<endl;
	argv[1][strlen(argv[1])-4] = 0;
	char sparsity[30];
	char Graph_size[30];
	string OutCor = string(argv[1]).append("_").append(string(itoa(N, Graph_size, 10)));
	string Outfilename;
	ofstream fout;
	
	bool * adjacent = new bool [(long long)N*N];
	for (k = 0; k < NumS; k++)
	{
		// Form full graph
		time = clock();
		Clength = FormFullGraph(adjacent, Cormat_cpu, N, S_thresh[k]);
		time = clock() - time;
//		cout<<"time for forming full graph (ms):\t" <<time<<endl;
		if (k == 0)
			C = new int [Clength];

		// form csr graph
		time = clock();
		FormCSRGraph(R, C, adjacent, N);
		time = clock() - time;
//		cout<<"time for forming csr graph (ms):\t" <<time<<endl;

		// write files
		sprintf_s(sparsity, "_%.2f%%", 100.0 * Clength / M1 / 2);
		Outfilename = OutCor;
		Outfilename.append(sparsity).append(".csr");
		fout.open(Outfilename.c_str(), ios::binary | ios::out);
		fout.write((char*)&Rlength, sizeof(int));
		fout.write((char*)R, sizeof(int) * Rlength);
		fout.write((char*)&Clength, sizeof(int));
		fout.write((char*)C, sizeof(int) * Clength);
		fout.close();
	}
/*
	delete []R;
	delete []C;
	delete []adjacent;
	delete []Cormat_cpu;
	delete []mask;
	delete []BOLD;
*/
	cout<<"==========================================================="<<endl;
	return 0;
}


int FormFullGraph(bool * adjacent, real__t * Cormat, int N, real__t threshold)
{
	int index = 0;
	int nonZeroCount = 0;
	memset(adjacent, 0, sizeof(bool) * N * N);
	for(int i = 0; i < N; i++)
		for(int j = i+1; j < N; j++)
		{
			if (fabs(Cormat[index]) > threshold)
			{
				nonZeroCount += (adjacent[i * N + j] = adjacent[j * N + i] = true);
			}
			index ++;
		}
	return nonZeroCount * 2;
}



void FormCSRGraph(int * R, int * C, bool * adjacent, int N)
{
	int count = 0;
	int Cindex = 0;
	R[0] = 0;
	for(int i = 0; i < N; i++)
 	{
		for(int j = 0; j < N; j++)
		{
			if(adjacent[i*N+j])
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
	// calculate moments
	int i = 0, j = 0, k = 0, l = 0;
	// transposing the BOLD signal
	real__t * BOLD_t = new real__t [L*N];
	for (i = 0; i < L; i ++)
		for (j = 0; j < N; j++)
		{
			BOLD_t[j*L+i] = BOLD[i*N+j];
		}

	// remove DC
	real__t * moment1 = new real__t [N];
	real__t * moment2 = new real__t [N];
	real__t * row;
	double sum1;
	double sum2;
	for (i = 0; i < N; i ++)
	{
		sum1 = sum2 = 0;
		row = BOLD_t + i*L;
		for (j = 0; j < L; j++)
		{
			sum1 += row[j];
			sum2 += row[j]*row[j];
		}
		moment1[i] = (real__t)sum1;
		moment2[i] = (real__t)sqrt(sum2 - sum1*sum1 / L);
		sum1 /= L;
		for (j = 0; j < L; j++)
			row[j] -= (real__t)sum1;
	}

	// compute the cross terms
	float sum3;
	for (k = 0, i = 0; i < N; i++)
		for (j = i+1; j < N; j++)
		{
			for (l = 0, sum3 = 0; l < L; l++)
			{
				sum3 += BOLD_t[i*L+l] * BOLD_t[j*L+l];
			}
			Cormat[k] = sum3 / moment2[i] / moment2[j];
			if (moment2[i] == 0 || moment2[j] == 0)
				Cormat[k] = 0;
			k++;
		}
	
	delete []BOLD_t;
	delete []moment1;
	delete []moment2;
}

void CorMat_cpu_MT(real__t * Cormat, real__t * BOLD, int N, int L)
{	
	// calculate moments
	int i = 0, j = 0, k = 0, l = 0;
	// transposing the BOLD signal
	real__t * BOLD_t = new real__t [L*N];
	for (i = 0; i < L; i ++)
		for (j = 0; j < N; j++)
		{
			BOLD_t[j*L+i] = BOLD[i*N+j];
		}

	// remove DC
	real__t * moment1 = new real__t [N];
	real__t * moment2 = new real__t [N];
	real__t * row;
	double sum1;
	double sum2;
	for (i = 0; i < N; i ++)
	{
		sum1 = sum2 = 0;
		row = BOLD_t + i*L;
		for (j = 0; j < L; j++)
		{
			sum1 += row[j];
			sum2 += row[j]*row[j];
		}
		moment1[i] = (real__t)sum1;
		moment2[i] = (real__t)sqrt(sum2 - sum1*sum1 / L);
		sum1 /= L;
		for (j = 0; j < L; j++)
			row[j] -= (real__t)sum1;
	}

	// compute the cross terms
	HANDLE *tHandle = new HANDLE[numThread];
	Corr_ARG *arg = new Corr_ARG[numThread];
	for (int i = 0; i < numThread; i++)
	{
		arg[i].id = i;
		arg[i].Cormat = Cormat;
		arg[i].Bold_t = BOLD_t;
		arg[i].moment2 = moment2;
		arg[i].N = N;
		arg[i].L = L;
	}
	for (int i = 0; i < numThread; i++)
	{
		Corr_ARG *temp = arg + i;
		//pthread_create(&thread[i], NULL, partial_CorMat, (voi*)(temp));
		tHandle[i] = (HANDLE) _beginthread(Cross_term, 0, (char *)temp);
	}	
	for (int i = 0; i < numThread; i++)
		WaitForSingleObject(tHandle[i], INFINITE);

	/*float sum3;
	for (k = 0, i = 0; i < N; i++)
		for (j = i+1; j < N; j++)
		{
			for (l = 0, sum3 = 0; l < L; l++)
			{
				sum3 += BOLD_t[i*L+l] * BOLD_t[j*L+l];
			}
			Cormat[k] = sum3 / moment2[i] / moment2[j];
			if (moment2[i] == 0 || moment2[j] == 0)
				Cormat[k] = 0;
			k++;
		}*/
	
	delete []BOLD_t;
	delete []moment1;
	delete []moment2;
}


void Cross_term(void *input_arg) 
{	
	Corr_ARG * temp = (Corr_ARG * )input_arg;
	int id = temp->id;
	real__t * Cormat = temp->Cormat;
	real__t * BOLD_t = temp->Bold_t;
	real__t * moment2 = temp->moment2;
	int N = temp->N;
	int L = temp->L;

	float sum3;
	int i, j, k, l;
	for (i = 0; i < N; i++)
		for (j = i+1+id; j < N; j += numThread)
		{
			for (l = 0, sum3 = 0; l < L; l++)
			{
				sum3 += BOLD_t[i*L+l] * BOLD_t[j*L+l];
			}
			k = (2*N-i-1)*i/2 + j - (i+1);
			Cormat[k] = sum3 / moment2[i] / moment2[j];
			if (moment2[i] == 0 || moment2[j] == 0)
				Cormat[k] = 0;
		}
	return;
}
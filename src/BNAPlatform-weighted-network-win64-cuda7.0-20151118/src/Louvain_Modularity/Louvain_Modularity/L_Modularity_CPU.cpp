# include <iostream>
# include <fstream>
# include <cmath>
# include <memory.h>
# include <cstring>
# include "Timer.h"
# include "dirent.h"   
# include <time.h> 
#include <Windows.h>
#include <process.h>
#include <ppl.h>
#include "Maslov.h"
using namespace Concurrency;
using namespace std;
//using System;
//using System.Threading;

#define ep  1e-10
typedef float real_t;

int seed = 1;
clock_t total_time = 0;
ofstream fout;

void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);

//int calc_k_GPU(long long n0, real_t *W, real_t *k);
//int calc_dQ_GPU(int n0, real_t *knmi, real_t *km, real_t ki, real_t Wii, int Mi, double s);
//int init_GPU_memory(int N);

typedef struct arg
{
	int id;
	int n0;
	int numThread;
	double *tQ;
	real_t *tknm_i;
	real_t *tkm;
	real_t k_i;	
	real_t Wii;
	int node_i;
	int Mi;
	double s0;
	int *max_index;
	double *max_Q;

} ARG;

int Compare(const void *elem1,const void *elem2)
{
	return*((int*)(elem1))-*((int*)(elem2));
}

void randperm(int *a, int n)
{	
	srand(seed++);
	for (int i=n-1; i>0; --i)
	{  
		int idx = (rand() * (0x7fff+1) + rand()) % (i+1);  
		if (idx == i) continue;  
		int tmp=a[idx];  
		a[idx] = a[i];  
		a[i] = tmp;
	}
}

double maxq(double *x, int n, int *idx)
{
	double max_x = x[0];
	*idx = 0;
	for (int i = 0; i<n; i++)
		if (max_x < x[i])
		{
			max_x = x[i];
			*idx = i;
		}
	return max_x;
}

int unique(int *M, int N)
{	
	int *M1 = new int [N];
	memcpy(M1,M,sizeof(int)*N);
	qsort(M1, N, sizeof(int),Compare);
	int idx = 0;
	for (int i=0; i<N; i++)
	{
		if (i>0 && M1[i]==M1[i-1])
				continue;
		
		for (int j = 0; j < N; j++)
			if (M[j]==M1[i])
				M[j] = idx;
		
		idx++;
	}
	return idx;
}

double Louvain_modularity_sub_parr( real_t **W0, long long *n1, int *Ci, long long N, double s)
{	
	real_t *W = *W0;
	int i=0,j=0;
	long long n0=*n1;
	double q_new = 0;
	real_t *k = new real_t [n0];
	memset(k,0,sizeof(real_t)*n0);
	//clock_t time = clock();
	parallel_for(long long(0), n0,  [&](long long ii) 
	{	
		double tmp_k = 0;
		for(int jj=0; jj<n0; jj++)
			tmp_k += *(W+ii*n0+jj);
		k[ii] = tmp_k;
	});	
	//time = clock() - time;
	//cout<<"time for calc_k = "<<time<<endl;
		
	real_t *km = new real_t [n0];
	memcpy(km, k, sizeof(real_t)*n0);
	/*memset(km,0,sizeof(real_t)*n0);
	//for(i=0; i<n0; i++) cout<<km[i]<<endl;
	for(i=0; i<n0; i++)
	{
		for(j=0; j<n0; j++)
			km[i] += *(W+i*n0+j);	
		//if (km[i] != k[i])
		printf("km[%d] = %f, k[%d] = %f\n",i, km[i], i, k[i]);
	}*/
	real_t *knm = new real_t [n0*n0];
	memcpy(knm, W, sizeof(real_t)*n0*n0);
	
	int *M = new int [n0];
	for (i = 0; i < n0; i++)
		M[i]=i;
		
	int *Nm = new int[n0];
	for (i = 0; i < n0; i++)
		Nm[i]=i;
	
	double *dQ = new double [n0];
    int cycle = 0;
	bool flag = true;
	while (flag)
	{
		flag=false;
		cout<<"cycle : "<<++cycle<<endl;
		clock_t time = clock();
		int *a = new int[n0];
		for (i=0; i<n0; i++)
			a[i] = i;
		randperm(a, n0);
		//for(int idx=0; idx<n0; idx++)
		//	cout<<a[idx]<<endl;

		
		for(int idx=0; idx<n0; idx++)
		{   i=a[idx];
		   
			
			//cout<<"i = "<<i<<endl;
			double tmp1 = W[i*n0+i] - knm[i*n0+M[i]] ;
			double tmp2 = k[i] - km[M[i]];
			
			//clock_t time = clock();
			/*parallel_for (long long(0), n0, [&](long long jj)
			{						
				dQ[jj]=(*(knm+i*n0+jj)+tmp1)-k[i]*(km[jj]+tmp2)*s;
			});*/
			for (j=0; j<n0; j++)
			{						
				dQ[j]=(*(knm+i*n0+j)+tmp1)-k[i]*(km[j]+tmp2)*s;
			}
			dQ[M[i]]=0;
			//time = clock() - time;
			//cout<<"time for calc_Q = "<<time<<endl;

			//for (j=0; j<n0; j++) cout<<"dQ "<<j<<" : "<<dQ[j]<<endl;

			double max_dQ = maxq(dQ,n0,&j);
			  
			//cout<<"max_dQ and id= "<<max_dQ<<' '<<j<<endl;
			if (max_dQ > ep)
			{				
				int z=0;
				for (z=0; z<n0; z++)
					*(knm+z*n0+j) += *(W+z*n0+i);
				for (z=0; z<n0; z++)
					*(knm+z*n0+M[i]) -= *(W+z*n0+i);

				km[j] += k[i];
				km[M[i]] -= k[i];

				Nm[j]++;
				Nm[M[i]]--;

				M[i]=j;
				flag = true;
			}
		}
		time = clock() - time;
		cout<<"cycle elapsed time : "<<time<<endl;
	}
	free(dQ);
		
	long long n = (long long) unique(M, n0); 
	//if (n==n0)       // No change
	//{cout<<"No change"<<endl;return -1;}
	//cout<<"n = "<<n <<endl;

	for (i=0; i<N; i++)
	{
		Ci[i] = M[Ci[i]];
		//cout<<Ci[i];
	}

	real_t *W1 = new real_t[n*n];
	memset(W1,0,sizeof(real_t)*n*n);

	for (i=0; i<n0; i++)
	{
		*(W1+M[i]*n+M[i]) += *(W+n0*i+i);
		for (j=i+1; j<n0; j++)
		{
			*(W1+M[i]*n+M[j]) += *(W+n0*i+j);
			*(W1+M[j]*n+M[i]) += *(W+n0*i+j);
		}
	}
	/*for (i=0; i<n; i++)
	{
		cout<<endl;
		for (j=0; j<n; j++)
			cout<<*(W1+i*n+j)<<' ';
	}*/

	free (*W0);
	*W0 = new real_t[n*n];
	memcpy(*W0, W1, sizeof(real_t)*n*n);
	free (W1);

	W=*W0;
		
	double s2=s*s;
	
	clock_t time1 = clock();
	
	for (i = 0;  i<n; i++)
	{	q_new +=*(W+n*i+i)*s;
		for (j=0; j<n; j++)
		{				
			for( int z=0; z<n; z++)
				q_new-=(*(W+i*n+z))*(*(W+z*n+j))*s2;
			//cout<<i<<' '<<j<<' '<<q_new<<endl;
		}
	}
	time1 = clock() - time1;
	cout<<"calc_q_new time : "<<time1<<endl;
	//total_time += time;

	*n1 = n;
	free(M);
	free(Nm);
	return q_new;
}

double Louvain_modularity_sub( real_t **W0, long long *n1, int *Ci, long long N, double s)
{	
	real_t *W = *W0;	
	
	int i=0,j=0;
	long long n0=*n1;
	//real_t *tmp = new real_t [n0];
	//memcpy(tmp, W+(n0-2)*n0,sizeof(real_t)*n0);

	double q_new = 0;
	real_t *k = new real_t [n0];
	memset(k,0,sizeof(real_t)*n0);
	
	//clock_t time = clock(); 
	for(i=0; i<n0; i++)
		for(j=0; j<n0; j++)
			k[i] += *(W+i*n0+j);
	//time = clock() - time;
	//cout<<"time for calc_k = "<<time<<endl;

	real_t *km = new real_t [n0];
	memcpy(km, k, sizeof(real_t)*n0);
	//for(i=0; i<n0; i++) cout<<km[i]<<endl;
		
	real_t *knm = new real_t [n0*n0];
	memcpy(knm, W, sizeof(real_t)*n0*n0);
	
	int *M = new int [n0];
	for (i=0; i<n0; i++)
		M[i]=i;
	
	int *Nm = new int[n0];
	for (i=0; i<n0; i++)
		Nm[i]=1;
	
	double *dQ = new double [n0];
        
	bool flag = true;
	int cycle = 0;
	while (flag)
	{
		flag=false;

		cout<<"cycle : "<<++cycle<<endl;
		clock_t time = clock();

		int *a = new int[n0];
		for (i=0; i<n0; i++)
			a[i] = i;
		randperm(a, (int) n0);
		//for(int idx=0; idx<n0; idx++)
		//	cout<<a[idx]<<endl;

		
		for(int idx=0; idx<n0; idx++)
		{   i=a[idx];
		   
			
			//cout<<"i = "<<i<<endl;
			//clock_t time = clock();
			for (j=0; j<n0; j++)
			{						
				dQ[j]=(*(knm+i*n0+j)-*(knm+i*n0+M[i])+*(W+i*n0+i))-k[i]*(km[j]-km[M[i]]+k[i])*s;
			}
			dQ[M[i]]=0;
			//time = clock() - time;
			//cout<<"time for calc_Q = "<<time<<endl;
			//for (j=0; j<n0; j++) cout<<"dQ "<<j<<" : "<<dQ[j]<<endl;

			double max_dQ = maxq(dQ,(int)n0,&j);
			  
			//cout<<"max_dQ and id= "<<max_dQ<<' '<<j<<endl;
			if (max_dQ > ep)
			{				
				int z=0;
				for (z=0; z<n0; z++)
					*(knm+z*n0+j) += *(W+z*n0+i);
				for (z=0; z<n0; z++)
					*(knm+z*n0+M[i]) -= *(W+z*n0+i);

				km[j] += k[i];
				km[M[i]] -= k[i];

				Nm[j]++;
				Nm[M[i]]--;

				M[i]=j;
				flag = true;
			}
		}
		time = clock() - time;
		cout<<"cycle elapsed time : "<<time<<endl;
	}
	free(dQ);
		
	int n = unique(M, (int)n0); 
	//if (n==n0)       // No change
	//{cout<<"No change"<<endl;return -1;}
	//cout<<"n = "<<n <<endl;

	for (i=0; i<N; i++)
	{
		Ci[i] = M[Ci[i]];
		//cout<<Ci[i];
	}

	real_t *W1 = new real_t[n*n];
	memset(W1,0,sizeof(real_t)*n*n);

	for (i=0; i<n0; i++)
	{
		*(W1+M[i]*n+M[i]) += *(W+n0*i+i);
		for (j=i+1; j<n0; j++)
		{
			*(W1+M[i]*n+M[j]) += *(W+n0*i+j);
			*(W1+M[j]*n+M[i]) += *(W+n0*i+j);
		}
	}
	/*for (i=0; i<n; i++)
	{
		cout<<endl;
		for (j=0; j<n; j++)
			cout<<*(W1+i*n+j)<<' ';
	}*/

	free (*W0);
	*W0 = new real_t[n*n];
	memcpy(*W0, W1, sizeof(real_t)*n*n);
	free (W1);

	W=*W0;
		
	double s2=s*s;
		
	for (i = 0;  i<n; i++)
	{	q_new +=*(W+n*i+i)*s;
		for (j=0; j<n; j++)
		{				
			for( int z=0; z<n; z++)
				q_new-=(*(W+i*n+z))*(*(W+z*n+j))*s2;
			//cout<<i<<' '<<j<<' '<<q_new<<endl;
		}
	}
	*n1 = n;
	free(M);
	free(Nm);
	return q_new;
}


double Louvain_modularity_sub1(  real_t **W0, int *R, int *C, float *V, long long *n1, int *Ci, long long N, double s)
{
	real_t *knm = *W0;
	//real_t *W = *W0;
	int i=0,j=0;
	long long n0=*n1;
	double q_new = 0;
	real_t *k = new real_t [n0];
	memset(k,0,sizeof(real_t)*n0);
	//clock_t time = clock();
	parallel_for(long long(0), n0,  [&](long long ii) 
	{
		for (j = R[ii]; j < R[ii+1]; j++)
			k[ii] += V[j];
	});
	/*for(i = 0; i < n0; i++) 
		for (j = R[i]; j < R[i+1]; j++)
		k[i] += V[j];
	*/		
	//time = clock() - time;
	//cout<<"time for calc_k = "<<time<<endl;
		
	real_t *km = new real_t [n0];
	memcpy(km, k, sizeof(real_t)*n0);
	/*memset(km,0,sizeof(real_t)*n0);
	//for(i=0; i<n0; i++) cout<<km[i]<<endl;
	for(i=0; i<n0; i++)
	{
		for(j=0; j<n0; j++)
			km[i] += *(W+i*n0+j);	
		//if (km[i] != k[i])
		printf("km[%d] = %f, k[%d] = %f\n",i, km[i], i, k[i]);
	}*/
	//real_t *knm = new real_t [n0*n0];
	//memcpy(knm, W, sizeof(real_t)*n0*n0);
	
	int *M = new int [n0];
	for (i = 0; i < n0; i++)
		M[i] = i;
	
	
	int *Nm = new int[n0];
	for (i = 0; i < n0; i++)
		Nm[i] = (int) i;
		
	double *dQ = new double [n0];
        
	bool flag = true;
	int cycle = 0;
	while (flag)
	{
		flag=false;
		
		cout<<"cycle : "<<++cycle<<endl;
		Setup(1);
		Start(1);
		int *a = new int[n0];
		for (i=0; i<n0; i++)
			a[i] = i;
		randperm(a, n0);
		//for(int idx=0; idx<n0; idx++)
		//	cout<<a[idx]<<endl;

		
		for(int idx=0; idx<n0; idx++)
		{   i=a[idx];
		   
			
			//cout<<"i = "<<i<<endl;
			double tmp1 = 0 - knm[i*n0+M[i]] ;
			double tmp2 = k[i] - km[M[i]];
			
			//clock_t time = clock();
			/*parallel_for (long long(0), n0, [&](long long jj)
			{						
				dQ[jj]=(*(knm+i*n0+jj)+tmp1)-k[i]*(km[jj]+tmp2)*s;
			});*/
					
			for (j = 0; j < n0; j++)
				dQ[j]=(*(knm+i*n0+j)+tmp1)-k[i]*(km[j]+tmp2)*s;
			
			dQ[M[i]]=0;
			//time = clock() - time;
			//cout<<"time for calc_Q = "<<time<<endl;

			//for (j=0; j<n0; j++) cout<<"dQ "<<j<<" : "<<dQ[j]<<endl;
			//int max_j;
			double max_dQ = maxq(dQ,n0,&j);
			
			//cout<<"max_dQ and id= "<<max_dQ<<' '<<j<<endl;
			
			if (max_dQ > ep)
			{				
				int z=0;
				
				for (z=R[i]; z<R[i+1]; z++)
				{
					*(knm+C[z]*n0+j) += V[z];
					*(knm+C[z]*n0+M[i]) -= V[z];
				}
					//for (z=0; z<n0; z++)
				//	*(knm+z*n0+M[i]) -= *(W+z*n0+i);
				
				km[j] += k[i];
				km[M[i]] -= k[i];

				Nm[j]++;
				Nm[M[i]]--;

				M[i] = j;
				flag = true;
			}
		}
		Stop(1);
		cout<<"Cycle elapsed time : "<<GetElapsedTime(1)<<'s'<<endl;
	}
	free(dQ);
		
	long long n = (long long) unique(M, n0); 
	//if (n==n0)       // No change
	//{cout<<"No change"<<endl;return -1;}
	//cout<<"n = "<<n <<endl;

	for (i=0; i<N; i++)
	{
		Ci[i] = M[Ci[i]];
		//cout<<Ci[i];
	}

	real_t *W1 = new real_t[n*n];
	memset(W1,0,sizeof(real_t)*n*n);

	for (i = 0; i < n0; i++)
		for (j=R[i]; j<R[i+1]; j++)
			*(W1+M[i]*n+M[C[j]]) += V[j];
	
	/*for (i=0; i<n; i++)
	{
		cout<<endl;
		for (j=0; j<n; j++)
			cout<<*(W1+i*n+j)<<' ';
	}*/

	free (*W0);
	*W0 = new real_t[n*n];
	memcpy(*W0, W1, sizeof(real_t)*n*n);
	free (W1);

	real_t *W=*W0;
		
	double s2=s*s;
	q_new = 0;
	
	clock_t time1 = clock();
	
	
	/*
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
			q_new += 1.0 * (Ci[i] == Ci[C[j]]);
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			q_new -= 1.0 * (Ci[i] == Ci[j]) * (R[i+1]-R[i]) * (R[j+1]-R[j]) *s;
	q_new *= s;
	*/
	for (i = 0;  i<n; i++)
	{	q_new +=*(W+n*i+i)*s;
		for (j=0; j<n; j++)
		{				
			for( int z=0; z<n; z++)
				q_new-=(*(W+i*n+z))*(*(W+z*n+j))*s2;
			//cout<<i<<' '<<j<<' '<<q_new<<endl;
		}
	}
	time1 = clock() - time1;
	cout<<"calc_q_new time : "<<time1<<endl;
	//total_time += time;
	*n1 = n;
	free(M);
	free(Nm);
	return q_new;
}

double Louvain_modularity( int *Ci, real_t *W, int *R, int *C, float *V, long long N, int *Num_module)
{   
	double s = 0;
	double q = 0, q_new = 0;
	long long n0 = N;
	long long i=0,j=0;
	
	for (i = 0; i<R[N]; i++)
		s += V[i];
	//cout<<"edge number : "<<s<<endl;
	s = 1/s;
	
	for (i=0; i<N; i++)
		Ci[i]=i;
	
	i=1;
	q_new = Louvain_modularity_sub1( &W, R, C, V, &n0, Ci, N, s);  
	cout <<"Round "<<i<<endl;
	fout <<"Round "<<i++<<endl;
	cout <<"Number of modules : "<<n0<<endl;
	fout <<"Number of modules : "<<n0<<endl;
	cout <<"Modularity Q : "<<q_new<<endl<<endl;
	fout <<"Modularity Q : "<<q_new<<endl<<endl;

	while (true)
	{
		q=q_new;
		q_new=Louvain_modularity_sub_parr(&W, &n0, Ci, N, s);
		//q_new=Louvain_modularity_sub(&W, &n0, Ci, N, s);
		//q_new=Louvain_modularity_GPU_sub(&W, &n0, Ci, N, s);
		
		/*for (i=0; i<n0; i++)
		{
		cout<<endl;
		for (j=0; j<n0; j++)
			cout<<*(W+i*n0+j)<<' ';
		}*/

		//if (q_new < 0)   //No change from the last merging
		//{  q_new = q; break;  }
		cout <<"Round "<<i<<endl;
		fout <<"Round "<<i++<<endl;
		cout <<"Number of modules : "<<n0<<endl;
		fout <<"Number of modules : "<<n0<<endl;
		cout <<"Modularity Q : "<<q_new<<endl<<endl;
		fout <<"Modularity Q : "<<q_new<<endl<<endl;
		if (q_new-q < ep)
			break;
	}
	*Num_module = n0;
	return q_new;	
}


int main(int argc, char * argv[])
{
	ofstream flog("BNA_time_log", ios::app);
	clock_t total_time = clock();
	if (argc < 2) 
	{
		cerr<<"Input format: .\\Modularity.exe dir_for_csr num_of_random_networks \nFor example: .\\Modularity_CPU.exe d:\\data 10"<<endl;
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

	for (i = 0; i < FileNumber; i++)
	{
		string a = string(argv[1]).append("\\").append(filename[i]);
		cout<<"\nModular analysis for "<<a.c_str()<<" ..."<<endl;
		ifstream fin(a.c_str(), ios_base::binary);
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}
		// Read x.csr
		int Rlength = 0, Clength = 0, Vlength = 0;
		
		fin.read((char*)&Rlength, sizeof(int));
		int * R = new int [Rlength];
		fin.read((char*)R, sizeof(int) * Rlength);

		fin.read((char*)&Clength, sizeof(int));
		int * C = new int [Clength];
		fin.read((char*)C, sizeof(int) * Clength);
		
		fin.read((char*)&Vlength, sizeof(int));
		float * V = new float [Vlength];
		fin.read((char*)V, sizeof(float) * Vlength);
		
		fin.close();
		
		long long N;
		N = Rlength - 1;
		cout<<"Number of voxel = "<<N<<endl;
		int *index = new int [N];
		memset(index, 0 , sizeof(int)*N);
		long long i,j;
		int idx = 0;
		long long N_C;
		int Nm;
		for (i = 0; i < N; i++)
			if(R[i+1]-R[i] != 0)
				index[i] = idx++;
		N_C = idx;	
		cout<<"Number of connected voxel = "<<N_C<<endl;
		//N=8000;
		//long long N2 = N * N; 
		long long N2 = N_C * N_C;
		real_t *W = new real_t[N2];
		memset(W, 0, sizeof(real_t)*N2);
		
		int *R_new = new int [N_C+1];
		int *C_new = new int [Clength];
				
		for (i = 0; i < N; i++)
		{
			if(R[i] == R[i+1]) continue;
			R_new[index[i]] = R[i];
			for (j = R[i]; j < R[i+1] && C[j] < N; j++)
			{
				C_new[j] = index[C[j]];
				W[index[i]*N_C+C_new[j]] = V[j];
			}
		}
		R_new [N_C] = R[N];
		/*
		for (i = 0; i < N; i++)
			for (j = R[i]; j < R[i+1] && C[j] < N; j++)
				W[i*N+C[j]] = 1;
		*/		
		int *Ci = new int [N_C];
		//int *Ci = new int [N];

		string X_modu = a.substr(0, a.find_last_of('.')).append("_modu.nm");
		string X_cp_mas = a.substr(0, a.find_last_of('.')).append("_modu.txt");
			
		fout.open(X_cp_mas.c_str(), ios::out);	// Open the log file


		Setup(0);
		Start(0);
		double q = Louvain_modularity(Ci, W, R_new, C_new, V, N_C, &Nm);

		
		//cout<<"total_time for calc_q_new : "<<total_time<<endl;
		cout<<"q = "<<q<<endl;
		
		int *result = new int [N];
		idx = 0;
		for (i = 0; i < N; i++)
		{
			if (R[i+1] == R[i])
				result[i] = ++Nm;
			else 
				result[i] = Ci[idx++]+1;
		}

		fout<<"Final Modularity Q = "<<q<<endl;
		fout<<"Final number of modules : "<<Nm<<endl;
		fout<<"\n************************* Modularity of Random Networks *************************\n\n";
		/*
		int M = (R[N] - R[0])/ 2;	
		double Q = 0;
		for (i = 0; i < N_C; i++)
			for (j = R_new[i]; j < R_new[i+1]; j++)
				Q += 1.0 * (Ci[i] == Ci[C_new[j]]);
		for (i = 0; i < N_C; i++)
			for (j = 0; j < N_C; j++)
				Q -= 1.0 * (Ci[i] == Ci[j]) * (R_new[i+1]-R_new[i]) * (R_new[j+1]-R_new[j]) / 2 / M;
		Q /= 2 * M;
		cout<<"Q = "<<Q<<endl;	
		*/
		
		
		ofstream fresult;
		fresult.open(X_modu.c_str(), ios::binary|ios::out);
		
		float *result_f = new float [N];
		for (int k = 0; k < N; k++)
			result_f[k] = (float) result[k];
		fresult.write((char*)&N, sizeof(int));
		fresult.write((char*)result_f, sizeof(float) * N);

		
		int Maslov_num = atoi(argv[2]);
		cout<<"Modular analysis for random networks..."<<endl;

		int * R_dst = new int [N_C+1];
		int * C_dst = new int [Clength];
		float *V_dst = new float [Vlength];
		//Setup(0);
		//Start(0);
		
		for (int l = 0; l < Maslov_num; l++)
		{
			Maslov_weighted_1(R_dst, C_dst,V_dst, R_new, C_new, V, N_C+1, Clength);
			real_t *Wr = new real_t[N2];
			memset(Wr, 0, sizeof(real_t)*N2);

			for (long long ii = 0; ii < N_C; ii++)
			{
				//R_new[index[i]] = R[i];
				for (long long jj = R_dst[ii]; jj < R_dst[ii+1]; jj++)
				{
					//C_new[j] = index[C[j]];
					Wr[ii*N_C+C_dst[jj]] = V_dst[jj];
				}
			}
			fout<<"Random Network "<<l+1<<":"<<endl;
			q = Louvain_modularity(Ci, Wr, R_dst, C_dst, V_dst, N_C,&Nm);
			fout<<"Final Modularity Q = "<<q<<endl;
			fout<<"\n*****************************************************\n\n";
			//Partition(R_dst, C_dst, Modu_Result);
			//Partition_GPU(R_dst, C_dst, Modu_Result);
		}
		
		Stop(0);
		
		flog<<"Modularity(+maslov)\t"<<a.c_str()<<"CPU\tkernel time\t"<<GetElapsedTime(0)<<"s"<<endl;

		cout<<"total elapsed time : "<<GetElapsedTime(0)<<'s'<<endl;
		fout<<"total elapsed time : "<<GetElapsedTime(0)<<'s'<<endl;
		fresult.close();
		fout.close();
		delete [] R;
		delete [] R_new;
		delete [] C;
		delete [] C_new;
		delete [] V;
		delete [] R_dst;
		delete [] C_dst;
		delete [] V_dst;
		//system("pause");
	}
}
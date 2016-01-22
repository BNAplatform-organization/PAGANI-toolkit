# include <stdlib.h>
# include <stdio.h>
# include <iostream>
# include <fstream>
# include <cmath>
# include <iomanip>
# include <memory.h>
# include <cstring>
# include "Timer.h"
# include "dirent.h"   
# include "modularity_GPU.h" 
# include <time.h> 
using namespace std;

// This is an option for the initial value in Leading_Vector()
# define RANDOM_V0

// set some parameters

const int MAX_ITER = 10000;			// The maximum iteration times in the power method
const int ITERNUMBER =500;
const double BETA_Adjust = 0;		// An optional parameter for quicker convergence. Its effect is uncertain
const double Epsilon = 0.000001;	// If |x - x0| < Epsilon, quit iteraion 
const double LAMBDA = 0.01;		// if labmda > LAMBDA, initiate the division
const int MIN_GROUP = 1;			// The minimum nodes of an allowed module 

long long N;								// The Number of voxels
long long Ntemp;                              // The Number of voxels for each round
long long seed;
double * v, * v0, * verr, *vv;
double * sumBG;	// some buffers used in the iteration


ofstream fout;						// log file

// Forward Declaration
void Partition(int * R, int * C, int * Result);
bool Sub_Partition(int * OriR, int * OriC, int * R, int * C, int M, long long innerM, int * Result, int * Max_Result,int * AD);
double Lead_Vector(int * OriR, int * R, int * C, int M, double * sumBG1, double beta, int *AD);
template <class Type> double VectorNorm(Type * x, long long N);
void Maslov(int * R_dst, int * C_dst, int * R_src, int * C_src, int Rlength, int Clength);

/* 
This function does the partition, no return value.
R and C represent the adjacency matrix in CSR format.
Result stores the partition results.
*/
void Partition(int * R, int * C, int * Result)
{
	Setup(0);
	Start(0);
	int M = (R[N] - R[0]) / 2;			// The total Number connection in the network
	int Round = 0;						// The iteration round
	int Max_Result = 0;					// Maximum index of modules
	int Module_Num = 1;					// Used for adjust module index.
	bool Issub;							// Return by function Sub_Partition()
	bool * G = new bool [N];			// G[i] = 1 if node i is involved in this round of partition
	memset(Result, 0, sizeof(int) * N);
	int * Adjust_Result = new int [N];	// Map the results to consecutive intergers starting from 1 
	int Index = 1;						// Used in the adjusted results, starting from 1, increase by 1 at each successful division
	int * NewRow = new int [N];
	int * NewCol = new int [R[N]];
	int * Index_Result = new int [N];  // Used for matching the order number for every round partition
	long long i = 0, j = 0;
	int NumG = 0;						//
	long long innerM = 2*M;
	int Newtemp1 = 0;
	int	Newtemp2 = 0;
	while (Round <= Max_Result)	   //???????
	{
		innerM =2*M;
		NumG = 0;
		for (int i = 0; i < N; i++)
		{
			G[i] = (Result[i] == Round); 
			if (G[i])
			{
				Index_Result[NumG] = i;
				NumG ++;				// G[i] = 1 if node i is involved in this round
				continue;
			}
			innerM -= R[i+1] - R[i];
		}
		//Select all the involved node to form new row and col
		Newtemp1 = 0;
		Newtemp2 = 0;
		NewRow[Newtemp2] = Newtemp1;
		Newtemp2++;
		for(int i = 0;i < N; i++)
		{
			if(!G[i])
				continue;
			for (int j = R[i];j < R[i+1];j++)
			{
				if(!G[C[j]])
					continue;
				NewCol[Newtemp1] = C[j];
				Newtemp1++;
			}
			NewRow[Newtemp2] = Newtemp1;
			Newtemp2++;
		}
		Ntemp = Newtemp2 - 1;
		//main part of the partition
		if (NumG)						 
		{
			cout<<"\nRound:\t"<<Round<<'\t';
			cout<<"number of nodes:\t"<<NumG<<'\t';
			fout<<"\nRound:\t"<<Round<<'\t';
			fout.flush();
			fout<<"number of nodes:\t"<<NumG<<'\t';	
			Setup(1);
			Start(1);
			Issub = Sub_Partition(R,C,NewRow, NewCol, M, innerM, Result, &Max_Result, Index_Result);
										// call Sub_Partition() for this round of division
			Stop(1);
			cout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;
			fout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;
			if (!Issub)					// If divided, record the adjusted result
				Adjust_Result[Round] = Index++;
			Module_Num += Issub;		// Update the total number modules 
		}
		Round++;
	}
	Stop(0);

	// calculate Q
	double Q = 0;
	for (int i = 0; i < N; i++)
		for (int j = R[i]; j < R[i+1]; j++)
			Q += 1.0 * (Result[i] == Result[C[j]]);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Q -= 1.0 * (Result[i] == Result[j]) * (R[i+1]-R[i]) * (R[j+1]-R[j]) / 2 / M;

	Q /= 2 * M;

	cout<<"\nNumber of Modules: "<<Module_Num<<",\tQ="<<Q<<endl;
	fout<<"\nNumber of Modules: "<<Module_Num<<",\tQ="<<Q<<endl;
	cout<<"Elapsed time:   "<<GetElapsedTime(0)<<"s"<<endl;
	fout<<"Elapsed time:   "<<GetElapsedTime(0)<<"s"<<endl;

	// Adjust the results
	for (int i = 0; i < N; i++)	
		Result[i] = Adjust_Result[Result[i]];
	delete []Adjust_Result;
	delete []G;
	delete []Index_Result;
	return;
}

/* 
This function does the division of each round.
It returns a bool variable, indicating whether this round of partition is successful.
OriR represent the original adjacency matrix including every nodes of the graph.
R and C represent the new adjacency matrix  that only includes the nodes participate in this round
M is the total connection of the network..=
G indicates whether a node is involved in this round.
Result stores the partition results.
Max_Result is updated for outer function to decide whether to terminate partition.
IR is to help to match node order of this round with original node order.
*/
bool Sub_Partition(int * OriR, int * OriC, int * R, int * C, int M, long long innerM, int * Result, int * Max_Result,int * AD)
{
	long long i = 0, j = 0;
	// Initialize sumBG for this round
	for (i = 0; i < Ntemp; i++)
	{
		sumBG[i] = 0;
		sumBG[i] = R[i+1] - R[i] - (OriR[AD[i]+1] - OriR[AD[i]]) * (double)innerM / 2 / M;
	}
	/*for(int ii=0;ii<Ntemp;ii++)
	{
		if (sumBG[ii]!=0)
			cout<<"sumBG["<<ii<<"] = "<<sumBG[ii]<<endl;
    }*/
	// Call Lead_Vector() to calucate the most positive eigenvalue lambda
	double lambda = 0;
	lambda = Lead_Vector(OriR, R, C, M, sumBG, BETA_Adjust,AD);
	lambda -= BETA_Adjust;
	// If lambda < 0, calucate the leading eigenvalue for  B - lambda * I
	if (lambda < 0)
		lambda += Lead_Vector(OriR, R, C, M, sumBG, -lambda,AD);

	cout<<"Eigen Value: "<<lambda<<'\t';
	fout<<"Eigen Value: "<<lambda<<'\t';

	// Decide whether this round of partition is successful 
	long long subN = 0, subP = 0;
	for (i = 0; i < Ntemp; i++)
	{
		subP += (v[i] > 0);
		subN += (v[i] <= 0);
	}
	bool Issub = (lambda > LAMBDA && subP > MIN_GROUP && subN > MIN_GROUP);

	cout<<"Divide?: "<<Issub<<'\t';
	fout<<"Divide?: "<<Issub<<'\t';
	// If not divided, return; otherwise update Result and Max_Result
	if (!Issub)
		return 0;
	for (i = 0; i < Ntemp; i++)
		Result[AD[i]] = *Max_Result + 1 + (v[i] * (subP - subN + 0.5) <= 0);
	// notice: this is wrong  Result[i] = *Max_Result + 1 + (v[i] * (subP - subN) >= 0);
	(*Max_Result) += 2;
	return Issub;
}


/* 
This funtion calculates the most positive eigenvalue for matrix B - beta * I.
It returns the most positive value.

R and C represent the adjacency matrix in CSR format.
M is the total connection of the network.
G indicates whether a node is involved in this round.
For definition of sumBG1, see Sub-Partition() or [1].
For definition of beta, see the first line of this comment.

Global parameters v, v0, verr are also used in this function.
v and v0 are vectors for iteration. verr is their difference. 
*/
double Lead_Vector(int * OriR, int * R, int * C, int M, double * sumBG1, double beta, int *AD)
{
	long long i = 0, j = 0;
	// Initialize v. Two methods are optional. Define RANDOM_V0 if you want to use random starting vector
#ifdef RANDOM_V0
	srand(seed);
	for (i = 0; i < Ntemp; i++)
		v[i] = AD[i];
#else
	for (i = 0; i < N && !G[i]; i++)
		;
	v[i] = 1;
#endif

	double err1 = 1, err2 = 1;
	int ITER = 0;
	double vNorm = 0;
	double v_k;
	double temp1;

	while (err1 > Epsilon && err2 > Epsilon && ITER < MAX_ITER)
	{
		for (i = 0; i < Ntemp; i++)
		{
			v0[i] = v[i];
			vv[AD[i]] = v[i];
		}
		v_k = 0;
		// The dot product of v * k
		for (i = 0; i < Ntemp; i++)
			v_k += v0[i] * (OriR[AD[i]+1] - OriR[AD[i]]); //???
	
		// Do the matrix-vector multiplication
		for (i = 0; i < Ntemp; i++)
		{
			temp1 = 0;
			for (j = R[i]; j < R[i+1]; j++)
				temp1 += vv[C[j]];
			//cout<<"v["<<i<<"] = "<<temp1<<endl;
			temp1 -= v_k / 2 / M * (OriR[AD[i]+1] - OriR[AD[i]]) + (sumBG1[i] - beta) * v0[i];
			v[i] = temp1;
		}
		/* for(int ii=0;ii<Ntemp;ii++)
	   {
			if (v[ii]!=0)
				cout<<"v["<<ii<<"] = "<<v[ii]<<endl;
	   }*/
	  /* int ii=0;
	   while (v[ii]==0)
	   {
		   ii++;
	   }
	   cout<<"v["<<ii<<"] = "<<v[ii];*/
		// Decide whether converge
		vNorm = VectorNorm(v, N);
		for (i = 0; i < Ntemp; i++)
			v[i] = v[i] / vNorm;
		for (i = 0; i < Ntemp; i++)
			verr[i] = v[i] - v0[i];
		err1 = VectorNorm(verr, N);
		for (i = 0; i < Ntemp; i++)
			verr[i] = v[i] + v0[i];
		err2 = VectorNorm(verr, N);
		ITER++;
	}
	cout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';
	fout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<min(err1, err2)<<'\t';

	// return the eigenvalue
	long long max_index = 0;
	for (i = 0; i < Ntemp; i++)
		if (fabs(v[i]) > fabs(v[max_index]))
			max_index = i;
	return (v[max_index] * v0[max_index] > 0) ? vNorm: -vNorm;
}

/* 
This function returns the norm of the input vector x[G].
G is the logic subscriber and N is the matrix dimension.
*/
template <class Type>
double VectorNorm(Type * x, long long N)
{
	long long i = 0;
	double Norm = 0;
	for (i = 0; i < Ntemp; i++)
		Norm += x[i] * x[i];
	return sqrt(Norm);
}






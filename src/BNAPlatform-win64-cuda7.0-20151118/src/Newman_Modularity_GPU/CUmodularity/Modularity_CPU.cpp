# include <iostream>
# include <fstream>
# include <cmath>
# include <memory.h>
# include <cstring>
# include <sstream>
# include "Timer.h"
# include "dirent.h"   
# include <time.h> 
using namespace std;

// This is an option for the initial value in Leading_Vector()
// # define RANDOM_V0
typedef  unsigned int u_int;
// set some parameters

extern const int MAX_ITER = 10000;			// The maximum iteration times in the power method
extern const double BETA_Adjust = 0.1;		// An optional parameter for quicker convergence. Its effect is uncertain
extern const double Epsilon = 1e-5;	// If |x - x0| < Epsilon, quit iteraion 
extern const double LAMBDA = 0.001;		// if labmda > LAMBDA, initiate the division
extern const double DQ_MIN = 1e-10;
extern const int MIN_GROUP = 1;			// The minimum nodes of an allowed module 
extern ofstream fout;						// log file

bool Sub_Partition_GPU_test(int N, int * ind, int * R, int * C, double *K, double *sumBG, double m, int * Result, int * Num_module);

// Forward Declaration

/* 
This function returns the norm of the input vector x[G].
G is the logic subscriber and N is the matrix dimension.
*/
template <class Type>
double VectorNorm2(Type * x, int N)
{
	double Norm = 0;
	for (int i = 0; i < N; i++)
		Norm += (x[i] * x[i]);
	return sqrt(Norm);
}

template <class Type>
double VectorNorm(Type * x, int N)
{
	double Norm = 0;
	for (int i = 0; i < N; i++)
	{
		Norm = ( fabs(x[i]) > fabs(Norm) ) ? x[i] : Norm;
	}
	return (Norm);
}

/*
This function is equivalent to unique in MATLAB
*/
//int unique(int *M, int N)
//{	
//	int *M1 = new int [N];
//	memcpy(M1,M,sizeof(int)*N);
//	qsort(M1, N, sizeof(int),Compare);
//	int idx = 0;
//	for (int i=0; i<N; i++)
//	{
//		if (i>0 && M1[i]==M1[i-1])
//				continue;
//		
//		for (int j = 0; j < N; j++)
//			if (M[j]==M1[i])
//				M[j] = idx;
//		
//		idx++;
//	}
//	return idx;
//}

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
double Lead_Vector(int N, int * R, int * C, float *V, double *K,  double * sumBG, double m, double beta, double *u)
{
	long long i = 0, j = 0;
	//double *u = new double [N];
	
	//double *uerr = new double [N];
	// Initialize u. Two methods are optional. Define RANDOM_V0 if you want to use random starting vector
#ifdef RANDOM_V0
	srand(time(0));
	//srand(2016);
	u[0] = 1;                        //normalized;
	for (i = 1; i < N; i++)
		u[i] = rand()*1.0f/RAND_MAX;
	
#else
	for (i = 0; i < N; i++)
		u[i] = 1.0f*i/(N-1);
#endif
	
	double *u0 = new double [N];
	//double err1 = 1, err2 = 1;
	double err = 1;
	u_int ITER = 0;
	double uNorm = 0, uNorm0 = 0;
	double u_k;
	double temp1;

	while (err > Epsilon && ITER < MAX_ITER)
	{
		//for (i = 0; i < N; i++)
		//	v0[i] =  v[i];
		memcpy(u0,u,sizeof(double)*N);
		u_k = 0;
		
		// The dot product of k'*u, 
		for (i = 0; i < N; i++)
			u_k += u0[i] * (K[i]);
		
		//cout<<"u_k = "<<u_k<<endl;

		// Do the matrix-vector multiplication
		for (i = 0; i < N; i++)
		{
			temp1 = 0;
			for (j = R[i]; j < R[i+1]; j++)
				temp1 += V[j] * u0[C[j]];
			temp1 -= u_k / m * (K[i]) + (sumBG[i] + beta) * u0[i];
			u[i] = temp1;
		}
		//ofstream debugfile;
		//debugfile.open("debug_u",ios::binary|ios::out);
		//debugfile.write((char *)u, N*sizeof(double));
		//debugfile.close();
		
		//minus the second item: k*(k'*u)/m;  
		/*for (i = 0; i < N; i++)
			u[i] -= u_k*K[i]/m ;*/
		//minus the third item: (diag(SumBG)+beta*I)*u		
		/*for (i = 0; i < N; i++)
			u[i] -= (sumBG[i]+beta)*u0[i];*/
	
		uNorm = VectorNorm<double>(u, N);
		

		// Decide whether converge, using infinity norm
		err = fabs(uNorm-uNorm0);
		//cout<<uNorm<<" - "<<uNorm0<<" residual: "<<err<<endl;
		uNorm0 = uNorm;
		for (i = 0; i < N; i++)
			u[i] = u[i] / uNorm ;
		/*for (i = 0; i < N; i++)
			uerr[i] = u[i] - u0[i];
		err1 = VectorNorm(uerr, N);
		for (i = 0; i < N; i++)
			uerr[i] = u[i] + u0[i];
		err2 = VectorNorm(uerr, N);*/

		ITER++;
	}
	cout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<err<<'\t';
	fout<<"Iterations:\t"<<ITER<<'\t'<<"residual:\t"<<err<<'\t';
	delete []u0;
	u0 = NULL;
	// return the eigenvalue
	//u_int max_index = 0;
	//u[0] = G[0] ? v[0] : 0;
	//for (i = 0; i < N; i++)
	//	if (fabs(u[i]) > fabs(u[max_index]))
	//		max_index = i;
	return (uNorm);
}

/* 
This funtion calculates dQ for each round to decide whether to split the current submodule.
It returns a double variable dQ.
*/
double calculate_dQ(signed char *S, int N, int * R, int * C, float *V, double *K, double *sumBG, double m)
{
	double dQ = 0;
	//double * x = new double [N];
	int i,j;
	double temp;
	double k_s = 0;

	// x = S'*Bsub*S = S'*(Bsub*S) = S'*[(Asparse-k*k'/m)*S]
	for (i = 0; i < N; i++)
			k_s += S[i] * K[i];
	
	//cout<<"k_s = "<<k_s<<endl;
	for (i = 0; i < N; i++)
		{
			temp = 0;
			for (j = R[i]; j < R[i+1]; j++)
				temp += V[j] * S[C[j]];			
			//dQ += temp;
			dQ += S[i] * temp - sumBG[i];
		}
	//cout<<"pre dQ = "<<dQ<<endl;
	
	//for (i = 0; i < N; i++)
	//	 dQ += x[i] * S[i];
	dQ -= k_s*k_s/m;

	//delete []x;
	return (dQ);
}

double qmax(int N, double *Qit, char *indSub, int *imax)
{
	double Qmax = Qit[0];
	for(int i = 1; i < N; i++)
		if (Qmax < Qit[i] && indSub[i])
		{
			Qmax = Qit[i];
			*imax = i;
		}
	return Qmax;
}

double fine_tune_S(double dQ, int N, int * R, int * C, float * V, double * K, double m, signed char *S)
{
	signed char *Sit = new signed char [N];
	memcpy(Sit,S,sizeof(signed char)*N);
	signed char *S_Si = new signed char [N];
	double *Qit = new double [N];
	memset(Qit,0,sizeof(double)*N);
	char  * indSub = new char [N];
	fill(indSub,indSub+N,1); 
	double Qmax = dQ;
	double temp;
	double k_s;
	int i,j,k;
	int imax;
	int ITER = 0;
	double Q = dQ;
	//bool flag = TRUE;
	while (ITER<N)
	{
		
		for (i = 0; i < N; i++)
		{
			if (indSub[i] == 0)
				continue;

			for (j = 0; j < N; j++)
				S_Si[j] = Sit[j]*Sit[i];

			//Sit[k] = -Sit[k];
			S_Si[i] = 0;
			//k_s = 0;
			////////////////////////////////////////////////////
			//calculate Qit(k)=(Sit')*( Bsub-diag(Bsub) )*Sit
			//for (i = 0; i < N; i++)
			//	k_s += S[i] * K[i];

			temp = 0;
			for (j = R[i]; j < R[i+1]; j++)
				temp += S_Si[C[j]] * V[j];			
			//	Qit[k] += S[i] * temp + K[i]*K[i]/m;
			//}
			for (j = 0; j < N; j++)
				temp -= S_Si[j]*K[i]*K[j]/m;
			Qit[i] = Qmax - 4*temp;
			
			//cout<<Qit[k]<<'\t';
			////////////////////////////////////////////////////
			//Sit[k] = -Sit[k];
		}

		
		//for (i = 0; i < N; i++)
		//	Qit[i] *= indSub[i];
		
		Qmax = qmax(N, Qit, indSub, &imax);
		//cout<<"Qmax = "<<Qmax<<"\t imax = "<<imax<<endl;
		Sit[imax] = -Sit[imax];
		indSub[imax] = 0;
		if (Qmax > Q)
		{
			Q = Qmax;
			memcpy(S,Sit,sizeof(signed char)*N);			
		}
		else 
			break;
		ITER++;
	}
	delete []Sit;
	delete []Qit;
	delete []indSub;
	cout<<"fine tune Q: "<<Q<<endl;
	return Q;
}

/* 
This funtion does the division of each round.
It returns a bool variable, indicating whether this round of partition is successful.
R and C represent the adjacency matrix in CSR format.
M is the total connection of the network..=
G indicates whether a node is involved in this round.
Result stores the partition results.
Max_Result is updated for outer function to decide whether tto terminate partition.
*/
bool Sub_Partition(int N, int * ind, int * R, int * C, float *V, double *K, double *sumBG, double m, int * Result, int * Num_module)
{
	long long i = 0, j = 0;
	//long long temp1, temp2;
	double * eigv = new double [N];
	signed char *S = new signed char [N];
	// Call Lead_Vector() to calucate the most positive eigenvalue lambda
	double lambda = 0;
	lambda = Lead_Vector(N, R, C, V, K,  sumBG, m, BETA_Adjust, eigv);
	lambda += BETA_Adjust;
	// If lambda < 0, calucate the leading eigenvalue for  B - lambda * I
	if (lambda < 0)
		lambda += Lead_Vector(N, R, C, V, K,  sumBG, m, lambda, eigv);

	cout<<"Eigen Value: "<<lambda<<'\t';
	fout<<"Eigen Value: "<<lambda<<'\t';

	// Decide whether this round of partition is successful 
	long long subN = 0, subP = 0;
	for (i = 0; i < N; i++)
	{
		subP += (eigv[i] >= 0);  		
		subN += (eigv[i] < 0);	 		
		S[i] = ((eigv[i] >= 0) ? 1 : -1);
	}
	cout<<"subP: "<<subP<<endl;
	cout<<"subN: "<<subN<<endl;
	//calculate dQ;
	double dQ = 0;
	dQ = calculate_dQ(S, N, R, C, V, K, sumBG, m);
	cout<< "dQ = "<<dQ<<endl;
	
	//                fine tune results           //
	if (dQ>DQ_MIN)  
	{
		dQ = fine_tune_S(dQ, N, R, C, V, K, m, S);	
	}
	subP = 0;
	subN = 0;
	for (i = 0; i < N; i++)
	{
		subP += (S[i] >= 0);  		
		subN += (S[i] < 0);	 		
	}
	//dQ = calculate_dQ(S, N, R, C, V, K, sumBG, m);
	cout<< "after fine tune, dQ = "<<dQ<<endl;
	///////////////////////////////////////////////

	bool isSplit = (dQ > DQ_MIN && subP > MIN_GROUP && subN > MIN_GROUP);

	cout<<"Divide?: "<<isSplit<<'\t';
	fout<<"Divide?: "<<isSplit<<'\t';
	// If not divided, return; otherwise update Result and Max_Result
	if (isSplit)
	{
	//	for (i = 0; i < N; i++)
	//		Result[ind[i]] = *Max_Result + 1 + (eigv[i] * (subP - subN + 0.5) <= 0);
	//// notice: this is wrong  Result[i] = *Max_Result + 1 + (v[i] * (subP - subN) >= 0);
	//	(*Max_Result) += 2;
		(*Num_module) += 1;
		if (subP>subN)
			for (i = 0; i < N; i++)
				Result[ind[i]] = ((S[i] >= 0 ) ? Result[ind[i]] : (*Num_module));
		else
			for (i = 0; i < N; i++)
				Result[ind[i]] = ((S[i] < 0 ) ? Result[ind[i]] : (*Num_module));
	}
	delete []S;
	delete []eigv;
	return isSplit;
}


/* 
This function does the partition, no return value.
R and C represent the adjacency matrix in CSR format.
Result stores the partition results.
*/
void Partition(long long N, int * R, int * C, int * Result)
{
	long long M = R[N];
	long long i = 0,j = 0;
	
	double * K;
	K = new double [N];
	double m = 0;
	memset(K, 0, sizeof(double)*N);
	for (i = 0; i < N; i++){
		K[i] = R[i + 1] - R[i];
	}
	//double m_inv = 1.0/m;
	cout<<"check K: "<<K[0]<<'\t'<<K[N/2]<<'\t'<<K[N-1]<<endl; //check
	cout<<"check m: "<<m<<endl;
	//memset(Result, 0, sizeof(int) * N);
	fill(Result, Result+N, 1); //initialize Result
	//int * Adjust_Result = new int [N];   //can be optimized,old version bfs style
	int result_idx = 1;
	int Num_module = 1;
	int NumG = N;
	int * index = new int [N];
	
	int * ind = new int [N];            //later try to put in the loop, using NumG
	for (i = 0; i < N; i++)
		ind[i] = i;
	//int NumG = 0;
	int Round = 1;						// The iteration round
	//int Max_Result = 1;					// Maximum index of modules
	int * R_new = new int [N+1];
	memcpy(R_new, R, sizeof(int)*(N+1));
	int * C_new = new int [M];
	memcpy(C_new, C, sizeof(int)*(M));
	double *K_new = new double [N];
	memcpy(K_new, K, sizeof(double)*N);
	double *sumBG = new double [N];
	memset(sumBG, 0, sizeof(double)*N);
	
	bool isSplit;
	int ITER = 0;
	while (Round <= Num_module)
	{		
		/**********************************************************************************/
		/************************************ Partition ***********************************/
		if (NumG>1)
		{
			cout<<"\nRound:\t"<<Round<<'\t'<<"Iter:\t"<<ITER<<endl;
			cout<<"number of nodes:\t"<<NumG<<'\t';
			cout<<"number of non-zero elements:\t"<<R_new[NumG]<<'\t';
			cout<<"density of this submodule:\t"<<R_new[NumG]*1.0/( (double) NumG * NumG)<<'\t';
			fout<<"\nRound:\t"<<Round<<'\t';
			fout<<"number of nodes:\t"<<NumG<<'\t';
			fout.flush();
			
			Setup(1);
			Start(1);
			//isSplit = Sub_Partition(NumG, ind, R_new, C_new, V_new, K_new, sumBG, m, Result, &Num_module); //if split, Num_module+1;
			isSplit = Sub_Partition_GPU_test(NumG, ind, R_new, C_new, K_new, sumBG, m, Result, &Num_module); //if split, Num_module+1;
			Stop(1);
			
			cout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;
			fout<<"sub_partition time:   "<<GetElapsedTime(1)<<"s"<<endl;
			//if (!isSplit)					// If divided, record the adjusted result
			//	Adjust_Result[Round] = result_idx++;
			//Num_module += isSplit;		// Update the total number modules , old version
		}
		if (!isSplit)
			Round++;
		/**********************************************************************************/
		/************************** Find the next sub_module ***************************/
		NumG = 0;
		fill(index, index+N, -1);
		for (i = 0; i < N; i++)
		{
			if (Result[i] == Round)
			{	
				ind[NumG] = i;
				index[i] = NumG++;   // index[i] >= 0 if node i is involved in this round
				//NumG ++;
			}				
		}
		if (!NumG)
		{
			cout<<"No voxel in the submodule next round\n";
			cout<<"Round: "<<Round<<", and Num_module: "<<Num_module<<" should be equal.\n";
			continue;
		}
		int ii = 0;
		int jj = 0;
		R_new[0] = 0;
		for(i = 0;i < N; i++)
		{
			if(index[i] < 0)
				continue;
			K_new[ii] = K[i];                 
			for (j = R[i];j < R[i+1];j++)
			{
				if(index[C[j]] < 0)
					continue;
				C_new[jj] = index[C[j]];				
				if(C_new[jj] > NumG)                      //check flag
					cout<<C_new[jj]<<'\t'<<C[j]<<"C_new exceed NumG!\n";				
				jj++;
			}
			R_new[++ii] = jj;			
		}
		if (ii!=NumG)    //check flag
			cout<<"sub module voxel# not match!";

		//update R_new, C_new, V_new, and k (i.e., bsub);
		/**********************************************************************************/
		/******************************** diag(sum(bsub)) *********************************/
		double temp1 = 0, temp2 = 0;
		for (j = 0; j < NumG; j++)
				temp2 += (K_new[j]);
		
		for (i = 0; i < NumG; i++)
		{
			//if (!G[i])
				//continue;
			//sumBG[i] = 0;
			temp1 = R_new[i+1]-R_new[i];
			sumBG[i] = temp1 - K_new[i] * temp2 / m;
		}
		
		ITER++;

		//ofstream debugfile;
		//ostringstream s1;
		//
		//s1<<"Round"<<Round<<"_Iter"<<ITER;
		////string debugfilename = "round";
		//debugfile.open(s1.str().append("_ind"),ios::binary|ios::out);
		//debugfile.write((char *)ind, NumG*sizeof(int));
		//debugfile.close();
		//int Rlength = NumG+1;
		//int Clength = R_new[NumG];
		//debugfile.open(s1.str().append("_csr"),ios::binary|ios::out);
		//debugfile.write((char *)&Rlength, sizeof(int));
		//debugfile.write((char *)R_new, Rlength*sizeof(int));
		//debugfile.write((char *)&Clength, sizeof(int));
		//debugfile.write((char *)C_new, R_new[NumG]*sizeof(int));
		//debugfile.write((char *)&Clength, sizeof(int));
		//debugfile.write((char *)V_new, R_new[NumG]*sizeof(float));
		//debugfile.close();
		//debugfile.open(s1.str().append("_k"),ios::binary|ios::out);
		//debugfile.write((char *)K_new, NumG*sizeof(double));
		//debugfile.close();
		//debugfile.open(s1.str().append("_diag"),ios::binary|ios::out);
		//debugfile.write((char *)sumBG, NumG*sizeof(double));
		//debugfile.close();

		
	}
	
	double Q = 0;
	for (i = 0; i < N; i++)
		for (j = R[i]; j < R[i+1]; j++)
			Q += (Result[i] == Result[C[j]]);
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			Q -= (Result[i] == Result[j]) * (K[i]) * (K[j]) /m;
	Q = Q / m;

	cout<<"\nNumber of Modules: "<<Num_module<<",\tQ="<<Q<<endl;
	fout<<"\nNumber of Modules: "<<Num_module<<",\tQ="<<Q<<endl;

	//for (i = 0; i < N; i++)	
	//	Result[i] = Adjust_Result[Result[i]];
	
	//delete []Adjust_Result;
	delete []R_new;
	delete []C_new;
	delete []K_new;
	delete []K;
	delete []sumBG;
	delete []ind;
	delete []index;
}




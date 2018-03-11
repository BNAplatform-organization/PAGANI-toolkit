#include<stdio.h>
#include<math.h>
#define N 4
#define eps 1e-6
#define KM 30
float MaxValue(float x[],int n)
{     float Max=x[0];     
      int i;     
      for (i=1;i<n;i++)         
     if(fabs(x[i])>fabs(Max))Max=x[i];         
	  return Max; 
}
void PowerMethod(float *A)
{     float U[N],V[N],EC[N],r1,r2,temp;     
      int i,j,k=0;     
	  for(i=0;i<N;i++)U[i]=1;     
	  while(k<KM)     
	  {         
		  k++;         
		  for(i=0;i<N;i++)         
		  {             
			  temp=0;             
			  for(j=0;j<N;j++)
				  temp+=*(A+i*N+j)*U[j];             
			  V[i]=temp;         //V=A*U
		  }         
		  for(i=0;i<N;i++)
		  {
			  U[i]=V[i]/MaxValue(V,N);       //U=V/Lamda  
			  EC[i]=abs(U[i]);                //取绝对值
		  }
			  if(k==1)
			  r1=MaxValue(V,N);         
		  else 
			  r2=MaxValue(V,N);         
		      if
			  (fabs(r2-r1)<eps)
			  break;         
		   r1=r2;                        //收敛条件
	  }          
	  printf("r=%f\n",r2);          
	  for(i=0;i<N;i++)
		  printf("x[%d]=%f\n",i+1,EC[i]);
} 
void main()
{     
	//float A[N][N]={{2,-1,0},{-1,2,-1},{0,-1,2}} ;     
	//float A[N][N]={{1,0,0,0},{1,2,0,0},{-3,-3,-1,0},{-2,-2,-2,-3}} ;     
	float A[N][N]={{0,0,1,1},{0,0,1,0},{1,3,0,1},{1,0,1,0}} ;     
	PowerMethod(A[0]);
} 
#include <queue>
#include <stack>
#include <iostream>
using namespace std;

void Betweenness_CPU(int* row, int* col, int numVertices, int numEdges, float *BC)
{	
	int 	*dist 		= new int[numVertices];
	float 	*sigma 		= new float[numVertices];
	float 	*delta 		= new float[numVertices];
	stack<int> S;
	queue<int> Q;	

	int *P = new int[numEdges];
	int *P_cnt = new int[numVertices];

	for (int i = 0; i < numVertices; i += 1) 
	{
		BC[i] = 0;
	}

	for (int s = 0; s < numVertices; s += 1) 
	{
		for(int j = 0; j <numEdges; j++)
			P[j] = 0;

		for (int j = 0; j < numVertices; j++) 
		{
			sigma[j] = 0;
			dist[j] = -1; 
			delta[j] = 0;
			P_cnt[j] = 0;
		}
		sigma[s] = 1;
		dist[s] = 0;

		Q.push(s);
		while (!Q.empty())
		{
			int v = Q.front();
			Q.pop();			
			S.push(v);
			for(int i = row[v]; i < row[v+1]; i++)
			{
				int w = col[i];
				if (dist[w] < 0)
				{
					dist[w] = dist[v] + 1;
					Q.push(w);
				}
				if (dist[w] == dist[v] + 1)
				{
					sigma[w] = sigma[w] + sigma[v];
					P[row[w] + P_cnt[w]] = v;
					P_cnt[w] ++;
				}
			}
		}

		while (!S.empty())
		{
			int w = S.top();
			S.pop();
			int k;
			for ( int i = 0; i < P_cnt[w]; i ++)
			{
				delta[P[row[w] + i]] = delta[P[row[w] + i]] + 1.0*sigma[P[row[w] + i]]/sigma[w]*(1 + delta[w]);
			}			
			if (w != s)
			{		
				BC[w] = BC[w] + delta[w];
			}
		}
		if (s%10 == 0) cout<<s<<" done. "<<endl;
	}	

	delete []P;
	delete []sigma;
	delete []dist;
	delete []delta;
	delete []P_cnt;
	return;
}
#include <stack>
#include <queue>
#include <list>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <Timer.h>
#include <windows.h>
#include <process.h>


using namespace std;

typedef struct arg
{
	int id;
	int numVertices;
	int numThread;
	int *row;
	int *col;
	int *deg;
	int *sort_indices;
	float *APSP;
	
} ARG;

void APSP(void *i);

double Lp_CPU(int *row, int *col, int numVertices, int numEdges);
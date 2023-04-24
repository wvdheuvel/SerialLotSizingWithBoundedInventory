#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ilcplex/ilocplex.h>
#include "randomc.h"
#include "arrays.h"

ILOSTLBEGIN

void DataRnd(int m, int n, double *d, double **K, double **p, double **h, double &U, TRandomMersenne &ran);
void CplexMLS(int m, int n, int L, double *d, double **K, double **p, double **h, double U, IloEnv env, IloCplex cplex);

using namespace std;

int main(int argc, char **argv) 
{
	// Definitions
	int i, j, L, m, n;	
	double *d, **K, **p, **h, U=0.0; 
	//TRandomMersenne ran = TRandomMersenne(142857);		
	TRandomMersenne ran = TRandomMersenne(1);		

	// Initialize data
	m=5;	// number of levels
	n=10;	// number of periods
	L=2;	// bottleneck level

	d=new double[n];
	SetDim2Dbl(K,m,n);
	SetDim2Dbl(p,m,n);
	SetDim2Dbl(h,m,n);

	IloEnv env;
	try
	{	IloCplex cplex(env);	
	
		do
		{	//Generate data
			DataRnd(m,n,d,K,p,h,U,ran);

			//Solve problem with Cplex		
			
			CplexMLS(m,n,L,d,K,p,h,U,env,cplex);
		} 
		while(m>0);
	
		cplex.end();	
	}
	catch (IloException& ex) 
	{	cerr << "Error: " << ex << endl;
	}
	catch (...) 
	{	cerr << "Error" << endl;
	}
	env.end();

	// delete arrays
	delete[] d;
	DeleteDim2Dbl(K,m,n);
	DeleteDim2Dbl(p,m,n);
	DeleteDim2Dbl(h,m,n);
					
	cin >> i;
	return 0;
}
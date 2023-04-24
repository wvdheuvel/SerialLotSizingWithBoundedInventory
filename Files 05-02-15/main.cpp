#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "randomc.h"
#include "arrays.h"
#include <ctime>
#include <cmath>
#include "solchar.h"

void DataRnd(int m, int n, int L, double *d, double **K, double **p, double **h, double &U, TRandomMersenne &ran);
void CplexMLS(int m, int n, int L, double *d, double **K, double **p, double **h, double U, sol_char *solMIP);
double solveZangwill(int m, int n, double *d, double **consCost, double **varCost, double **invCost);
double solveZangwill2(int m, int n, double *d, double **K, double **p, double **h);
void WriteData(int m, int n, int L, double *d, double **K, double **p, double **h, double U);
void ReadData(int &m, int &n, int &L, double *d, double **K, double **p, double **h, double &U, char *file);
double ConcatenateRN(int m, int n, int L, double *d, double **K, double **p, double **h, double U);

using namespace std;

int main(int argc, char **argv) 
{
	// Definitions
	int i, j, L, m, n, tmpInt, nIt, nM, nT, nL, startM, startT, startL;	
	double *d, **K, **p, **h, U, optDP, optDP2, eps=0.00001, runtime; 
	double runtimeMIP=0.0, runtimeDP=0.0, avgMIP, avgDP, ssMIP, ssDP, sMIP, sDP, delta;
	long cnt=0;
	bool readFile=false;
	//char file[40]="instance_error.txt";
	//char file[40]="instance with DS path m=2 n=3.txt";
	char file[40]="instance error 01.txt";
	//char file[40]="instance error m=3 n=3 4.txt";
	TRandomMersenne ran = TRandomMersenne(142857);		
	clock_t t_start, t_end;
	//TRandomMersenne ran = TRandomMersenne(1);		
	sol_char *solMIP;
	solMIP = new sol_char;
	
	ofstream output;		
	output.open("runtimes.txt");	
	if (!output) 
		cout << "Error met openen output-file" << endl;	

	// Initialize data
	startM=3;
	startT=4;
	startL=0;
	nM=3; 
	nT=3;
	m=startM*pow(2.0,nM-1);	// number of levels
	n=startT*pow(2.0,nT-1);	// number of periods
	L=1;	// bottleneck level
	nIt=100;

	d=new double[n];
	SetDim2Dbl(K,m,n);
	SetDim2Dbl(p,m,n);
	SetDim2Dbl(h,m,n);

	// create headers output file	
	output << "T" << ' ' << "T" << ' '; 
	for (int iT=0; iT<nT; iT++)
	{	n=startT*pow(2.0,iT);
		output << n << ' ' << n << ' '  << n << ' ' << n << ' ';		
	}
	output << endl;

	output << "alg" << ' ' << "alg" << ' '; 
	for (int iT=0; iT<nT; iT++)
		output << "MIP" << ' ' << "MIP" << ' '  << "DP" << ' ' << "DP" << ' ';
	output << endl;	

	output << "meas" << ' ' << "meas" << ' '; 
	for (int iT=0; iT<nT; iT++)
		output << "avg" << ' ' << "std" << ' '  << "avg" << ' ' << "std" << ' ';			
	output << endl;	

	// start loops for multiple experiments
	for (int iM=0; iM<nM; iM++)
	{	m=startM*pow(2.0,iM);
		//m=startM+iM;			

	for (int iL=0; iL<=iM+1; iL++)
	{	L=floor(1.0*iL/(iM+1)*(m-2));
		if (m-L<2)
		{	cout << "Error: BN level too close to end level" << endl;
			cin >> tmpInt;
		}		
		//L=11;
		output << m << ' ' << L+1 << ' ';

	for (int iT=0; iT<nT; iT++)
	{	n=startT*pow(2.0,iT);
		//n=startT+iT;	

	avgMIP=0.0; avgDP=0.0;		
	ssMIP=0.0; ssDP=0.0;
	sMIP=0.0, sDP=0.0;
	cnt=0;
	/*solMIP.double_sourcing=0;
	solMIP.reg_network=0;
	solMIP.obj=-1.0;*/
	
	do
	{	cnt++;
		
		//Generate data
		if (readFile)
			ReadData(m,n,L,d,K,p,h,U,file);
		else
		{	DataRnd(m,n,L,d,K,p,h,U,ran);
			WriteData(m,n,L,d,K,p,h,U);
		}

		//Solve problem with Cplex	
		t_start = clock();
		CplexMLS(m,n,L,d,K,p,h,U,solMIP);
		t_end = clock();
		runtime = 1.0*(t_end - t_start)/CLOCKS_PER_SEC;
		delta=runtime-avgMIP;
		avgMIP+=delta/cnt;
		sMIP+=delta*(runtime-avgMIP);
		cout << "Running time MIP: " << runtime << " sec." << endl;
		
		//Solve problem with Zangwill's DP by Eric
		/*t_start = clock();				
		optDP=solveZangwill(m,n,d,K,p,h);
		t_end = clock();		
		runtime = 1.0*(t_end - t_start)/CLOCKS_PER_SEC;
		cout << "Running time ZW DP1: " << runtime << " sec." << endl;*/

		//Solve problem with Zangwill's DP by Wilco
		/*t_start = clock();				
		optDP=solveZangwill2(m,n,d,K,p,h);
		t_end = clock();		
		runtime = 1.0*(t_end - t_start)/CLOCKS_PER_SEC;
		delta=runtime-avgDP;
		avgDP+=delta/cnt;
		sDP+=delta*(runtime-avgDP);
		cout << "Running time ZW DP2: " << runtime << " sec." << endl;	
		cout << "Optimal obj. DP: " << optDP << endl;*/

		//Apply general DP
		t_start = clock();		
		optDP=ConcatenateRN(m,n,L,d,K,p,h,U);
		t_end = clock();		
		runtime = 1.0*(t_end - t_start)/CLOCKS_PER_SEC;		
		delta=runtime-avgDP;
		avgDP+=delta/cnt;
		sDP+=delta*(runtime-avgDP);
		cout << "Optimal obj. DP: " << optDP << endl;
		cout << "Running time DP: " << runtime << " sec." << endl;

		if (abs(optDP-solMIP->obj)>eps)
		{	cout << "Error: MIP vs DP" << endl;
			cin >> tmpInt;
			WriteData(m,n,L,d,K,p,h,U);
		}

		if (solMIP->double_sourcing_far)
			WriteData(m,n,L,d,K,p,h,U);
	} 
	//while(cnt);
	while(cnt<nIt);
	cout << "Number of instances checked: " << cnt << endl;
	//cnt = 16 geeft interessante instantie

	//output << m << ' ' << n << ' ' << L+1 << ' ' << nIt << endl;
	if (cnt>1)
	{	sMIP=sqrt(sMIP/(cnt-1));
		sDP=sqrt(sDP/(cnt-1));
	}
	cout << "Average (st dev) runtime MIP: " << avgMIP << " (" << sMIP << ")" << endl;
	//output << "Average (st dev) runtime MIP: " << avgMIP << ' ' << sMIP << endl;		
	output << avgMIP << ' ' << sMIP << ' ';		
	cout << "Average (st dev) runtime DP: " << avgDP << " (" << sDP << ")" << endl;
	//output << "Average (st dev) runtime DP: " << avgDP << ' ' << sDP << endl;
	output << avgDP << ' ' << sDP << ' ';

	}
	output << endl;
	}
	}

	output.close();

	// delete arrays
	delete[] d;
	DeleteDim2Dbl(K,m,n);
	DeleteDim2Dbl(p,m,n);
	DeleteDim2Dbl(h,m,n);
	delete solMIP;
					
	cin >> i;
	return 0;
}
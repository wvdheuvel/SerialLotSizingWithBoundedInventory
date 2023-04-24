#include <iostream>
#include <fstream>
#include "randomc.h"
int cnt=0;

using namespace std;

void DataRnd(int m, int n, double *d, double **K, double **p, double **h, double &U, TRandomMersenne &ran)
{	int i,j;

	// demand
	for (j=0; j<n; j++)
	{	d[j]=ran.IRandom(1,20);
	}	

	// cost parameters
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	K[i][j]=2.0;
			p[i][j]=0.0;
			//h[i][j]=1.0;
			//h[i][j]=1.0*(i+1)/m;
			h[i][j]=1.0*pow(i+1,2.0)/pow(m,2.0);
		}
	}

	// bottleneck
	U=10.0;
}

void PrintData(int T, double* d, double *Ks, double *ps, double *hs, double *Kr, double *pr, double *tr, double *hr)
{	cout << "d: ";
	for (int i=0; i<T; i++)
		cout << d[i] << ' ';
	cout << endl;

	cout << "Ks: ";
	for (int i=0; i<T; i++)
		cout << Ks[i] << ' ';
	cout << endl;
	
	cout << "ps: ";
	for (int i=0; i<T; i++)
		cout << ps[i] << ' ';
	cout << endl;

	cout << "hs: ";
	for (int i=0; i<T; i++)
		cout << hs[i] << ' ';
	cout << endl;

	cout << "Kr: ";
	for (int i=0; i<T; i++)
		cout << Kr[i] << ' ';
	cout << endl;

	cout << "pr: ";
	for (int i=0; i<T; i++)
		cout << pr[i] << ' ';
	cout << endl;

	cout << "tr: ";
	for (int i=0; i<T; i++)
		cout << tr[i] << ' ';
	cout << endl;		

	cout << "hr: ";
	for (int i=0; i<T; i++)
		cout << hr[i] << ' ';
	cout << endl;		
}

void WriteData(int T, double* d, double *Ks, double *ps, double *hs, double *Kr, double *pr, double *tr, double *hr)
{	ofstream output;		
	output.open("instance.txt");	
	if (!output) 
		cout << "Error met openen output-file" << endl;	
	
	for (int i=0; i<T; i++)
		output << d[i] << ' ';
	output << endl;
	
	for (int i=0; i<T; i++)
		output << Ks[i] << ' ';
	output << endl;
	
	for (int i=0; i<T; i++)
		output << ps[i] << ' ';
	output << endl;

	for (int i=0; i<T; i++)
		output << hs[i] << ' ';
	output << endl;

	for (int i=0; i<T; i++)
		output << Kr[i] << ' ';
	output << endl;

	for (int i=0; i<T; i++)
		output << pr[i] << ' ';
	output << endl;

	for (int i=0; i<T; i++)
		output << tr[i] << ' ';
	output << endl;		

	for (int i=0; i<T; i++)
		output << hr[i] << ' ';
	output << endl;	

	output.close();
}

void ReadData(int T, double* d, double *Ks, double *ps, double *hs, double *Kr, double *pr, double *tr, double *hr, char *file)
{	ifstream input;		
	input.open(file);		
	if (!input) 
		cout << "Error met openen input-file" << endl;
	
	for (int i=0; i<T; i++)
		input >> d[i];
		
	for (int i=0; i<T; i++)
		input >> Ks[i];	
	
	for (int i=0; i<T; i++)
		input >> ps[i];

	for (int i=0; i<T; i++)
		input >> hs[i];

	for (int i=0; i<T; i++)
		input >> Kr[i];

	for (int i=0; i<T; i++)
		input >> pr[i];

	for (int i=0; i<T; i++)
		input >> tr[i];

	for (int i=0; i<T; i++)
		input >> hr[i];

	input.close();
}
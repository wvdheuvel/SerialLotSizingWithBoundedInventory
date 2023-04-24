#include <iostream>
#include <fstream>
#include <cmath>
#include "randomc.h"
#include "data.h"
int cnt=0;

using namespace std;

double round(double x, int n);

void DataRnd(int m, int n, int L, double *d, double **K, double **p, double **h, double &U, TRandomMersenne &ran)
{	int i,j;

	// demand
	for (j=0; j<n; j++)
	{	d[j]=ran.IRandom(1,20);
	}	

	// cost parameters
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	
			// setup cost
			//K[i][j]=2.0;
			if (i==L)
				K[i][j]=2.0;			
			//else if (i<L)
			//	K[i][j]=0.0;			
			else
			{	//K[i][j]=1+2*ran.Random();
				K[i][j]=round(4*ran.Random(),2);
			}
			//K[1][1]=100.0;
			//K[2][1]=100.0;

			// production cost
			if (i<=L)
			{	p[i][j]=1.0;
				//p[i][j]=0.0;
			}
			else
				p[i][j]=round(ran.Random(),2);
			
			// holding cost
			//h[i][j]=1.0;
			//h[i][j]=1.0*(i+1)/m;
			//h[i][j]=1.0*pow(i+1,2.0)/pow(m,2.0);
			//h[i][j]=0.01+1.0*pow(i,2.0)/pow(m-1,2.0);
			if (i<=L)
				h[i][j]=0.01*(i+1);
			else
				h[i][j]=0.01+1.0*(i-L)/(m-L-1);
		}
	}

	// bottleneck cap
	U=10.0;
	//U=n*20.0;
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

void WriteData(int m, int n, int L, double *d, double **K, double **p, double **h, double U)
{	int i, j;
	ofstream output;		
	output.open("instance.txt");	
	if (!output) 
		cout << "Error met openen output-file" << endl;	
	
	// levels and periods
	output << m << ' ' << n << ' ' << L << ' ' << U << endl;
	
	// demand
	for (i=0; i<n; i++)
		output << d[i] << ' ';
	output << endl;
	
	// cost parameters
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	output << K[i][j] << ' ';
		}
		output << endl;
	}

	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	output << p[i][j] << ' ';
		}
		output << endl;
	}

	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	output << h[i][j] << ' ';
		}
		output << endl;
	}

	output.close();

}

void ReadData(int &m, int &n, int &L, double *d, double **K, double **p, double **h, double &U, char *file)
{	int i,j;
	double tmpDbl;
	ifstream input;		
	input.open(file);		
	if (!input) 
	{	cout << "Error met openen input-file" << endl;
		cin >> tmpDbl;
	}

	// levels and periods
	input >> m >> n >> L >> U;
	
	// demand
	for (i=0; i<n; i++)
		input >> d[i];
	
	// cost parameters
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
			input >> K[i][j];		
	}

	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
			input >> p[i][j];
	}

	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
			input >> h[i][j];
	}

	input.close();
}

void CopyData(int m, int n, int L, double *d, double **K, double **p, double **h, double U, Data *(&inst))
{	int i,j;
	double tmpDbl;

	// bottleneck and capacity
	inst->L=L;
	inst->U=U;
	
	// demand
	for (i=0; i<n; i++)
		inst->d[i]=d[i];
	
	// cost parameters
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	inst->K[i][j]=K[i][j];
			inst->p[i][j]=p[i][j];
			inst->h[i][j]=h[i][j];
		}
	}

}

/*double round(double x, int n)
{	double y=pow(10.0,n);
	double z=floor(x*y+0.5);
	return z/y;
}*/

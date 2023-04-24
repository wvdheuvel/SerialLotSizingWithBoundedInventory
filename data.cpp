#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include "data.h"
#include "randomc.h"
#include <string>

double round(double x, int n);

using namespace std;

Data::Data(int mTmp, int nTmp)
{	// create vectors/matrices of right length
	m=mTmp;
	n=nTmp;
	d = new double[n];
	SetDim2Dbl(K,m,n);
	SetDim2Dbl(p,m,n);
	SetDim2Dbl(h,m,n);
}

Data::Data(char *file)
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

	d = new double[n];
	SetDim2Dbl(K,m,n);
	SetDim2Dbl(p,m,n);
	SetDim2Dbl(h,m,n);
	
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


Data::~Data()
{	// delete vectors/matrices
	delete[] d;
	DeleteDim2Dbl(K,m,n);
	DeleteDim2Dbl(p,m,n);
	DeleteDim2Dbl(h,m,n);
}

void Data::DataRnd(TRandomMersenne &ran)
{
	int i, j;
	double baseK = 8.0;

	// demand
	for (j = 0; j < n; j++)
	{
		d[j] = ran.IRandom(1, 20);
	}

	// cost parameters
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			// setup cost
			//K[i][j]=2.0;
			if (i == L)
				K[i][j] = baseK;
			//else if (i<L)
			//	K[i][j]=0.0;			
			else
			{	//K[i][j]=1+2*ran.Random();
				K[i][j] = round(2 * baseK*ran.Random(), 2);
			}
			//K[1][1]=100.0;
			//K[2][1]=100.0;

			// production cost
			if (i <= L)
			{
				p[i][j] = 1.0;
				//p[i][j]=0.0;
			}
			else
			{	//p[i][j]=round(ran.Random(),2);
				p[i][j] = 1.0;
			}

			// holding cost
			//h[i][j]=1.0;
			//h[i][j]=1.0*(i+1)/m;
			//h[i][j]=1.0*pow(i+1,2.0)/pow(m,2.0);
			//h[i][j]=0.01+1.0*pow(i,2.0)/pow(m-1,2.0);
			if (i <= L)
				h[i][j] = 0.01*(i + 1);
			else
				h[i][j] = 1.0*(i - L) / (m - L - 1);
		}
	}

	// bottleneck cap
	U = 10.0;
	//U=n*20.0;
}

void Data::WriteData(int k)
{	int i, j;
	ofstream output;		
	string file = "Instances/instance";
	file.append("T");
	file.append(to_string(n));
	file.append("L");
	file.append(to_string(m));
	file.append("B");
	file.append(to_string(L));
	file.append("I");
	file.append(to_string(k));
	file.append(".txt");
	//output.open("Instances/instance.txt");	
	output.open(file);
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

void Data::PrintData()
{	int i, j;
	
	cout << "No of periods: " << n << endl;
	cout << "No of levels: " << m << endl;
	cout << "Bottleneck level: " << L << endl;
	cout << "Bottleneck capacity: " << U << endl;

	// demand
	cout << "demand:" << endl;
	for (i=0; i<n; i++)
		cout << d[i] << ' ';
	cout << endl;
	
	// cost parameters	
	cout << "setup cost:" << endl;
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	cout << K[i][j] << ' ';
		}
		cout << endl;
	}

	cout << "production cost:" << endl;
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	cout << p[i][j] << ' ';
		}
		cout << endl;
	}

	cout << "holding cost:" << endl;
	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	cout << h[i][j] << ' ';
		}
		cout << endl;
	}

}

void Data::CopyData(int & _L, double *_d, double **_K, double **_p, double **_h, double & _U)
{	int i,j;
	
	_L=L;

	for (j=0; j<n; j++)
	{	_d[j]=d[j];
	}

	for (i=0; i<m; i++)
	{	for (j=0; j<n; j++)
		{	_K[i][j]=K[i][j];
			_p[i][j]=p[i][j];
			_h[i][j]=h[i][j];
		}
	}

	_U=U;

}

double round(double x, int n)
{	double y=pow(10.0,n);
	double z=floor(x*y+0.5);
	return z/y;
}
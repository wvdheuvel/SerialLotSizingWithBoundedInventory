#include <ilcplex/ilocplex.h>
#include "definitions.h"

ILOSTLBEGIN


void CplexMLS(int m, int n, int L, double *d, double **K, double **p, double **h, double U, IloEnv env, IloCplex cplex)
{	//IloEnv env;
	try 
	{	IloInt i, j;
		double opt;
		
		// Variables		
		NumVarArray2 x, Inv;
		x = CreateNumVarArray2(env, m, n, "x", 0);
		Inv = CreateNumVarArray2(env, m, n, "I", 0);
		BoolVarArray2 y;
		y = CreateBoolVarArray2(env, m, n, "y");
		
		IloModel model(env);				
		
		// Constraints
		
		// 1b + 1e: inventory balance at all but retailer's level
		for (i=0; i<m-1; i++) // starting inv is zero
		{	model.add(Inv[i][0]==x[i][0]-x[i+1][0]);
		}
		for (i=0; i<m-1; i++)
		{	for (j=1; j<n; j++)
			{	model.add(Inv[i][j]==Inv[i][j-1]+x[i][j]-x[i+1][j]);
			}
		}				
		
		// 1c: inventory balance at retailer level
		model.add(Inv[m-1][0]==x[m-1][0]-d[0]);
		for (j=1; j<n; j++)
		{	model.add(Inv[m-1][j]==Inv[m-1][j-1]+x[m-1][j]-d[j]);
		}

		// 1d: bottleneck constraint
		for (j=0; j<n; j++)
		{	model.add(Inv[L-1][j]<=U);
		}				
				
		// setup forcing
		// cumulative demand
		double *M;
		M=new double[n];
		M[n-1]=d[n-1];		
		for (j=n-2; j>=0; j--)
			M[j]=M[j+1]+d[j];
		
		//constraints
		for (i=0; i<m; i++)
		{	for (j=0; j<n; j++)
			{	model.add(x[i][j]<=M[j]*y[i][j]);
			}
		}		
		delete[] M;		

		// Objective function
		IloExpr obj(env);
		for (i=0; i<m; i++)
		{	for (j=0; j<n; j++)
			{	obj+=K[i][j]*y[i][j]+p[i][j]*x[i][j]+h[i][j]*Inv[i][j];
			}
		}	
		model.add(IloMinimize(env, obj));
		obj.end();
	  
		// Solve model
		//IloCplex cplex(env); already created in main
		cplex.extract(model);
		cplex.solve();			
				
		// Output
		// Optimal solution
		opt=cplex.getObjValue();
		cplex.out() << endl << "Optimal solution: " << opt << endl;
		
		// show figure
		IloNum tolerance = cplex.getParam(IloCplex::EpInt);	
		cplex.out() << endl;
		for (i=0; i<m; i++)
		{	for (j=0; j<n; j++)
			{	//cplex.out() << (cplex.getValue(y[i][j])>1-tolerance) << ' ';
				if (cplex.getValue(y[i][j])>1-tolerance)
					cplex.out() << '|' << ' ';
				else
					cplex.out() << ' ' << ' ';
			}
			cplex.out() << endl;
			for (j=0; j<n; j++)
			{	//cplex.out() << ' ' << (cplex.getValue(Inv[i][j])>tolerance);
				cplex.out() << char(250);
				if (cplex.getValue(Inv[i][j])>tolerance)
					cplex.out() << '-';
				else
					cplex.out() << ' ';
			}
			cplex.out() << endl;
		}					
		cplex.out() << endl;

		// find special instances
		// count number of 'flow splits'
		int cnt=0;
		for (i=0; i<m-1; i++)
		{	for (j=0; j<n; j++)
			{	if (cplex.getValue(x[i+1][j])*cplex.getValue(Inv[i][j])>tolerance)
					cnt++;
			}
		}
		if (cnt>5)
		{	cout << "Interesting: flow split" << endl;
		}

		// count number of double sourcing periods
		cnt=0;
		for (i=0; i<m; i++)
		{	for (j=1; j<n; j++)
			{	if (cplex.getValue(x[i][j])*cplex.getValue(Inv[i][j-1])>tolerance)
				{	cnt++;
					if (i>L)
						cout << "Interesting: double sourcing far from bottleneck" << endl;
				}
			}
		}
		if (cnt>0)
		{	cout << "Interesting: double sourcing" << endl;			
		}

		model.end();
	} 
	catch(IloException& e) 
	{	cerr  << " ERROR: " << e << endl;   
	}
	catch(...) 
	{	cerr  << " ERROR" << endl;   
	}
	//env.end(); already ended in main
}
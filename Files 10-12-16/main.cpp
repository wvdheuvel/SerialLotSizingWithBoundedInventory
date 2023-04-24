#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <cmath>
#include "randomc.h"
#include "arrays.h"
#include "solchar.h"
#include "data.h"
#include "cplexMLS.h"
#include "dps.h"
#include "stats.h"

using namespace std;

int main(int argc, char **argv) 
{
	// Definitions and initialization of variables
	int i, L, m, n, tmpInt, nIt, nM, nT, startM, startT, startL; // see variable explanation below	
	double optDP, eps=0.00001, runtime; 
	long cnt=0; // counts number of instances solved
	bool readFile=false;	//true if data read from file, otherwise random data
	bool DPonly=false;		//true if only DP is solved (not MIP)
	//char file[40]="instance_error.txt";
	//char file[40]="instance error m=3 n=3 4.txt";
	//char file[40]="instance counter example LT6.txt"; // name input file when data read from file
	char file[50]="instance counter example LT6 T=5 m=5 NS.txt";
	//TRandomMersenne ran = TRandomMersenne(142857); // initialize random generator	
	TRandomMersenne ran = TRandomMersenne(1);	
	clock_t t_start, t_end;	// needed for measuring running times
	Data *inst;			// problem instance -> see class "data"
	sol_char *solMIP;	// solution characteristics -> see class "solchar"
	Stats *statsMIP, *statsDP;	// running time statistics -> see class "stats"
	
	// output file name for recording running times		
	ofstream output;	
	output.open("runtimes.txt");	
	if (!output) 
		cout << "Error: opening output file" << endl;	

	// Initialize parameter settings
	nM=3;		// no of different stages/levels
	nT=3;		// no of different time horizons
	nIt=10;		// no of instances per setting
	startM=3;	// starting number of levels
	startT=8;	// starting number of periods
	startL=1;	// starting bottleneck level
	L=1;		// bottleneck level (zero is upper level)
	m=startM*pow(2.0,nM-1);	// max number of levels
	n=startT*pow(2.0,nT-1);	// max number of periods
	
	if (readFile) // overwrite settings -> only one iteration is needed
	{	nM=1;
		nT=1;
		nIt=1;
		startL=1;
	}

	// create headers for output file	
	output << "T" << ' ' << "T" << ' '; 
	for (int iT=0; iT<nT; iT++)
	{	n=startT*pow(2.0,iT);
		output << n << ' ' << n << ' '  << n << ' ' << n << ' ';		
	}
	output << endl;

	output << "alg" << ' ' << "alg" << ' '; 
	for (int iT=0; iT<nT; iT++)
		output << "MIP" << ' ' << "MIP" << ' '  << "DP" << ' ' << "DP" << ' ' << "DS" << ' ';
	output << endl;	

	output << "m" << ' ' << "L" << ' '; 
	for (int iT=0; iT<nT; iT++)
		output << "avg" << ' ' << "std" << ' '  << "avg" << ' ' << "std" << ' ' << "DS" << ' ';			
	output << endl;	

	// start loops for multiple experiments
	for (int iM=0; iM<nM; iM++) 
	{	m=startM*pow(2.0,iM); // set no of stages
		//m=startM+iM;			

	for (int iL=startL; iL<=iM+1; iL++)
	{	// set bottleneck stage
		//L=floor(1.0*iL/(iM+1)*(m-2)); 
		L=floor(1.0*iL/(iM+2)*m);
		if (m-L<2)
		{	L=m-2;
			//cout << "Error: BN level too close to end level" << endl;
			//cin >> tmpInt;
		}		
		//L=2;
		output << m << ' ' << L+1 << ' ';

	for (int iT=0; iT<nT; iT++)
	{	n=startT*pow(2.0,iT); // set no of periods
		//n=startT+iT;	

	// iniatilize statistics
	cnt=0;
	statsMIP = new Stats();
	statsDP = new Stats();
	
	do
	{	cnt++;		

		// Read or generate data		
		if (readFile)
		{	inst = new Data(file);
			inst->PrintData();			
		}
		else
		{	inst = new Data(m,n);
			inst->L=L;
			inst->DataRnd(ran);			
			//inst->PrintData();			
			//inst->WriteData();
		}			
		solMIP = new sol_char;	

		//Solve problem with Cplex	
		if (!DPonly)
		{	t_start = clock();
			solveMLSbyCplex(inst,solMIP);
			t_end = clock();
			runtime = 1.0*(t_end - t_start)/CLOCKS_PER_SEC;
			statsMIP->Add(runtime,solMIP->double_sourcing);
			cout << "Running time MIP: " << runtime << " sec." << endl;
		}
		
		//Solve uncapacitated problem with Zangwill's DP written by Wilco
		/*t_start = clock();				
		optDP=solveZangwill2(m,n,d,K,p,h);
		t_end = clock();		
		runtime = 1.0*(t_end - t_start)/CLOCKS_PER_SEC;
		cout << "Running time ZW DP2: " << runtime << " sec." << endl;	
		cout << "Optimal obj. DP: " << optDP << endl;*/

		//Apply general DP
		t_start = clock();
		optDP=solveMLSbyDP(inst);
		t_end = clock();		
		runtime = 1.0*(t_end - t_start)/CLOCKS_PER_SEC;		
		statsDP->Add(runtime,0);
		cout << "Optimal obj. DP: " << optDP << endl;
		cout << "Running time DP: " << runtime << " sec." << endl;

		// give error when DP deviates from MIP
		if (!DPonly)
		{	if (solMIP->obj - optDP < -eps)
			{	cout << "Error: MIP vs DP" << endl;		
				inst->WriteData();
				cin >> tmpInt;
			}
		}

		// write instance to file when interesting
		if (solMIP->double_sourcing_far)
		{	inst->WriteData();		
		}

		delete inst;
		delete solMIP;	

		cout << "Number of instances checked: " << cnt << endl << endl;
	} 
	while(!readFile && cnt<nIt);
	//while(cnt<nIt);
	
	//cnt = 16 gives interesting instance

	//write output to screen and files
	//output << m << ' ' << n << ' ' << L+1 << ' ' << nIt << endl;
	cout << "Average (st dev) runtime MIP: " << statsMIP->getAvg() << " (" << statsMIP->getStDev() << ")" << endl;
	output << statsMIP->getAvg() << ' ' << statsMIP->getStDev() << ' ';	
	cout << "Average (st dev) runtime DP: " << statsDP->getAvg() << " (" << statsDP->getStDev() << ")" << endl;
	output << statsDP->getAvg() << ' ' << statsDP->getStDev() << ' ' << statsMIP->getDblSrc() << ' ';
	cout << "Number of instances with double sourcing: " << statsMIP->getDblSrc() << endl;

	delete statsMIP;
	delete statsDP;

	}
	output << endl;
	}
	}

	output.close();
					
	cin >> i;
	return 0;
}
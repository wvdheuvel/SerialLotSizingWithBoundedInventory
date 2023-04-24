#include <iostream>
#include <cmath>
#include "arrays.h"
#include "data.h"

using namespace std;

void cumDemands(int n, double *demand, double **matrix);
void cumHolding(int m, int n, double **cumD, double **h, double **cumH);
void costZWtrees(int m, int n, double *d, double **K, double **p, double **h, double ****cost);
void BottleneckTrees(int m, int n, double *d, double **K, double **p, double **h, double U, double ****costZW, double ****costBN);
void DownstreamTrees(int m, int n, double *d, double **K, double **p, double **h, double U, double ****costZW, double ****costBN, double ****costDS);
void OptDownstreamPath(int m, int n, double *d, double **K, double **p, double **h, double U, int b1, int t1, int t2, double ****costZW, double ****costBN, double ****costDS);
void UpstreamTrees(int L, int n, double *d, double **K, double **p, double **h, double ****costZW, double ****costBN, double ****costDS, double ****costUS);
void SetConstDim4Dbl(double ****x, int n1, int n2, int n3, int n4, double c);
void SetConstDim6Dbl(double ******x, int n1, int n2, int n3, int n4, int n5, int n6, double c);
void cumHoldU(int n, double U, double **h, double **cumH);
int maxI(int x, int y);
int minI(int x, int y);
void PrintCost(int n1, int n2, int n3, int n4, double cst, char* str);
void CopyPart(int m1, int m2, int n, double **x, double **y);

//double DPMLS(int m, int n, int L, double *d, double **K, double **p, double **h, double U)
double solveMLSbyDP(Data *inst)
{	// initialize variables
	double opt=-1.0, tmpCst, minCst, tmpDbl;
	double inf=pow(10.0,10.0);
	bool checkUsed=false;
	
	// copy data
	int m, n, L;
	double *d, **K, **p, **h, U;
	
	m=inst->m;
	n=inst->n;	
	d=new double[n];
	SetDim2Dbl(K,m,n);
	SetDim2Dbl(p,m,n);
	SetDim2Dbl(h,m,n);

	inst->CopyData(L,d,K,p,h,U);

	//create DP variables
	double **f, **cumD;									
	double **KDS, **pDS, **hDS;
	double ****costZW, ****costBN, ****costDS, ****costUS;
	int **OptIndb, **OptIndt;
	
	SetDim2Dbl(cumD,n,n);		// cumulative demands
	SetDim2Dbl(f,n+1,n+1);		//f[b][t] min cost to satisfy demands up till t 
								//using AT MOST bottleneck periods up till b
	SetDim2Dbl(KDS,m-L,n);		//copy for DS part
	SetDim2Dbl(pDS,m-L,n);		//copy for DS part
	SetDim2Dbl(hDS,m-L,n);		//copy for DS part
	SetDim2Int(OptIndb,n,n);	//optimal indices for backtraking purposes 
	SetDim2Int(OptIndt,n,n);	//optimal indices for backtraking purposes 
	SetDim4Dbl(costZW,m+1,n,n,n);//costs[i][j][t1][t2] min cost to satisfy demands [t1,t2] from node (i,j)	(dummy level 0 for convenience)
	SetDim4Dbl(costBN,n,n,n,n);  //costs[j1][j2][t1][t2] min cost to satisfy demands [t1,t2] through bottleneck periods [j1,j2]	
	SetDim4Dbl(costDS,n,n,n,n);  //costs[j1][j2][t1][t2] min cost to satisfy demands [t1,t2] through bottleneck periods [j1,j2]		
	SetDim4Dbl(costUS,n,n,n,n);
	
	if (checkUsed)
	{	SetConstDim4Dbl(costZW,m+1,n,n,n,-1.0);
		SetConstDim4Dbl(costBN,n,n,n,n,-1.0);
		SetConstDim4Dbl(costDS,n,n,n,n,-1.0);
	}
	
	CopyPart(m-L,m,n,K,KDS);
	CopyPart(m-L,m,n,p,pDS);
	CopyPart(m-L,m,n,h,hDS);
	cumDemands(n,d,cumD);
	
	costZWtrees(m-L,n,d,KDS,pDS,hDS,costZW);
	BottleneckTrees(m-L,n,d,KDS,pDS,hDS,U,costZW,costBN); 
	DownstreamTrees(m-L,n,d,KDS,pDS,hDS,U,costZW,costBN,costDS);
	UpstreamTrees(L,n,d,K,p,h,costZW,costBN,costDS,costUS);

	// arrays are not needed anymore
	DeleteDim4Dbl(costZW,m+1,n,n,n);		
	DeleteDim4Dbl(costBN,n,n,n,n);
	DeleteDim4Dbl(costDS,n,n,n,n);

	// Initialize
	f[0][0]=0;
	for (int t=1; t<=n; t++)
		f[0][t]=inf;

	// General DP
	for (int t2=0; t2<n; t2++)
	{	for (int b2=0; b2<=t2; b2++)
		{	// use exactly bottleneck periods up till b
			minCst=inf;
			for (int t1=0; t1<=t2; t1++)
			{	for (int b1=0; (b1<=t1)&&(b1<=b2); b1++)
				{	//BN trees only
					//tmpCst=f[b1][t1]+K[0][b1]+p[0][b1]*cumD[t1][t2]+costBN[b1][b2][t1][t2];
					
					//BN and DS paths  //replaced on 05-12-13
					/*tmpCst=K[L][b1]+p[L][b1]*cumD[t1][t2]+costBN[b1][b2][t1][t2];
					if ((b2>b1)&&(costDS[b1][b2][t1][t2]<tmpCst))
						tmpCst=costDS[b1][b2][t1][t2];

					if (tmpCst<costUS[b1][b2][t1][t2]-0.00001)
					{	cout << "Error: cost difference US/DS trees" << endl;
						cin >> tmpDbl;
						tmpDbl=tmpDbl;
					}*/

					tmpCst=costUS[b1][b2][t1][t2];
					tmpCst+=f[b1][t1];					
					
					//tmpCst=f[b1][t1]+costUS[b1][b2][t1][t2];

					/*if (tmpCst<1.0)
					{	cout << "Error: cost too small?" << endl;
						cin >> tmpDbl;
					}*/
					if (tmpCst<minCst)
					{	minCst=tmpCst;
						OptIndb[b2][t2]=b1;
						OptIndt[b2][t2]=t1;
					}
				}
			}
			
			// use at most bottleneck periods up till b
			if (f[b2][t2+1]<minCst)
			{	minCst=f[b2][t2+1];
				OptIndb[b2][t2]=-1;
				OptIndt[b2][t2]=t2+1;
			}
				

			// set min
			f[b2+1][t2+1]=minCst;
		}
	}

	// Find minimum
	int indOpt;
	opt=inf;
	for (int b=1; b<=n; b++)
	{	if (f[b][n]<opt)
		{	opt=f[b][n];
			indOpt=b-1;
		}
	}

	// Retrieve optimal RNs
	int indb=indOpt, indt=n-1, indbOld, indtOld;	
	while ((indb!=-1)&&(indt!=-1))
	{	indbOld=indb;
		indtOld=indt;
		indb=OptIndb[indbOld][indtOld]-1;
		indt=OptIndt[indbOld][indtOld]-1;
		if (indb==-2)
			indb=indbOld-1;		
		else
		{	//cout << indb+1 << ' ' << indbOld << ' ' << indt+1 << ' ' << indtOld << endl;		
		}
	}

	// delete arrays
	delete[] d;
	DeleteDim2Dbl(K,m,n);
	DeleteDim2Dbl(p,m,n);
	DeleteDim2Dbl(h,m,n);

	DeleteDim2Dbl(cumD,n,n);
	DeleteDim2Dbl(f,n+1,n+1);
	DeleteDim2Dbl(KDS,m-L,n);
	DeleteDim2Dbl(pDS,m-L,n);
	DeleteDim2Dbl(hDS,m-L,n);
	DeleteDim2Int(OptIndb,n,n);
	DeleteDim2Int(OptIndt,n,n);	
	DeleteDim4Dbl(costUS,n,n,n,n);
	return opt;
}

void DownstreamTrees(int m, int n, double *d, double **K, double **p, double **h, double U, double ****costZW, double ****costBN, double ****costDS)
{	double inf=pow(10.0,10.0), eps=0.00001, **cumD;
	SetDim2Dbl(cumD,n,n);
	cumDemands(n,d,cumD);
	
	/*for (int t2=0; t2<n; t2++)
	{	for (int b2=0; b2<=t2; b2++)
		{	for (int b1=0; b1<b2; b1++)
			{	for (int t1=b1; t1<=t2; t1++)
				{	//if ((b1==0) && (b2==2) && (t1==0) && (t2==2))
					//	cout << "Stop for debugging" << endl;
					
					if ((cumD[t1][t2]>U+eps)&&((m>3)||(t1==b1))) // case m=2 is specials case
						OptDownstreamPath(m,n,d,K,p,h,U,b1,t1,t2,costZW,costBN,costDS);
					else
						costDS[b1][b2][t1][t2]=inf;
				}
			}
		}
	}*/

	for (int t2=0; t2<n; t2++)
	{	for (int b1=0; b1<t2; b1++)
		{	for (int t1=b1; t1<=t2; t1++)
			{	/*if ((b1==0) && (b2==2) && (t1==0) && (t2==2))
					cout << "Stop for debugging" << endl;*/
				
				if ((cumD[t1][t2]>U+eps)&&((m>3)||(t1==b1))) // case m=2 is specials case
					OptDownstreamPath(m,n,d,K,p,h,U,b1,t1,t2,costZW,costBN,costDS);
				else
				{	for (int b2=b1+1; b2<=t2; b2++)
						costDS[b1][b2][t1][t2]=inf;
				}
			}
		}
		
	}

	DeleteDim2Dbl(cumD,n,n);
}

void OptDownstreamPath(int m, int n, double *d, double **K, double **p, double **h, double U, int b1, int t1, int t2, double ****costZW, double ****costBN, double ****costDS)
{	double inf=pow(10.0,10.0), eps=0.00001, tmpCst, minCst;
	int nDir=4, tUB, tmpI;
	double ****costDSP, **cumD, **cumHU;
	bool checkUsed=false;
	SetDim2Dbl(cumD,n,n);
	SetDim2Dbl(cumHU,n,n);
	cumDemands(n,d,cumD);
	cumHoldU(n,U,h,cumHU);
	SetDim4Dbl(costDSP,m,n,nDir,n+1);	//NB: extra period needed for cum demands
	if (checkUsed)
		SetConstDim4Dbl(costDSP,m,n,nDir,n+1,-1.0); // costDSPopt(i,j,dir,t) cost up til (i,j) with next node specified 
													// by dir and demand satisfied up till t-1	
	
	// compute UB on cum demand for presink nodes
	// NB is this correct?
	tUB=t1-1;
	while (cumD[t1][tUB+1]<=cumD[t1][t2]-U+eps)
		tUB++;
	if (abs(cumD[t1][tUB]-cumD[t1][t2]+U)<eps)
		tUB--;

	if (tUB>=t1)
	{

	/*if ((b1==0) && (b2==1) && (t1==0) && (t2==3))
		cout << "Stop for debugging" << endl;*/

	//Initialize pre-sink arcs
	
	// type 0 arc (b1,1)-(b1,2)
	costDSP[1][b1][0][t1]=K[0][b1]+p[0][b1]*cumD[t1][t2];	

	// type 0 arcs at period b1 and cum demands up to t1-1
	// NB arcs before t1 are allowed at final level
	for (int i=2; i<m; i++)
	{	costDSP[i][b1][0][t1]=costDSP[i-1][b1][0][t1]+K[i-1][b1]+p[i-1][b1]*(cumD[t1][t2]-U);
		if ((checkUsed) && (costDSP[i-1][b1][0][t1]<-eps))
		{	cout << "Error: var not set 00" << endl;
			cin >> tmpI;
			tmpI=tmpI;
		}
		if (checkUsed && (cumD[t1][t2]-U<-eps))
		{	cout << "Error: negative flow 01" << endl;
			cin >> tmpI;	
			tmpI=tmpI;
		}
	}		

	// type 1 arcs at period b1 before final level
	for (int i=1; i<m-1; i++)
	{	if (b1<t1)
		{	costDSP[i][b1][1][t1]=costDSP[i][b1][0][t1]+K[i][b1]+p[i][b1]*(cumD[t1][t2]-U);
			if ((checkUsed) && (costDSP[i][b1][0][t1]<-eps))
			{	cout << "Error: var not set 01" << endl;
				cin >> tmpI;
				tmpI=tmpI;
			}
			if (checkUsed && (cumD[t1][t2]-U<-eps))
			{	cout << "Error: negative flow 02" << endl;
				cin >> tmpI;	
				tmpI=tmpI;
			}
		}
		for (int t=t1+1; t<=t2+1; t++) // NB changed
		{	// case 8 immediate cost paper
			costDSP[i][b1][1][t]=costDSP[i][b1][0][t1]+K[i][b1]+p[i][b1]*(cumD[t1][t2]-U);		
			costDSP[i][b1][1][t]+=K[i+1][b1]+p[i+1][b1]*cumD[t1][t-1]+costZW[i+2][b1][t1][t-1];
			if ((checkUsed) && ((costDSP[i][b1][0][t1]<-eps)||(costZW[i+2][b1][t1][t-1]<-eps)))
			{	cout << "Error: var not set 02" << endl;
				cin >> tmpI;
				tmpI=tmpI;
			}
			if (checkUsed && (cumD[t1][t2]-U<-eps))
			{	cout << "Error: negative flow 03" << endl;
				cin >> tmpI;	
				tmpI=tmpI;
			}
		}
	}

	// horizontal arcs at period b1 at final level		
	if (b1==t1)
	{	costDSP[m-1][b1][1][t1+1]=costDSP[m-1][b1][0][t1]+K[m-1][b1]+p[m-1][b1]*(cumD[t1][t2]-U);
		if ((checkUsed) && (costDSP[m-1][b1][0][t1]<-eps))
		{	cout << "Error: var not set 03" << endl;
			cin >> tmpI;
			tmpI=tmpI;
		}
		if (checkUsed && (cumD[t1][t2]-U<-eps))
		{	cout << "Error: negative flow 04" << endl;
			cin >> tmpI;	
			tmpI=tmpI;
		}
	}

	//General DP pre-sink arcs
	for (int i=1; i<m-1; i++)
	{	for (int j=b1+1; j<=tUB; j++)
		{	//type 1 arcs
			for (int t=maxI(t1,j+1); t<=tUB+1; t++) // t=t1 changed 27-11-13
			{	minCst=inf;
				
				if (i>1) // else type 0 arcs do not exist 27-11-13
				{	// case 8 without ZW tree
					tmpCst=costDSP[i][j][0][t]+K[i][j]+p[i][j]*(cumD[t][t2]-U);
					if (tmpCst<minCst)
						minCst=tmpCst;
					if ((checkUsed) && (costDSP[i][j][0][t]<-eps))
					{	cout << "Error: var not set 04" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (checkUsed && (cumD[t][t2]-U<-eps))
					{	cout << "Error: negative flow 05" << endl;
						cin >> tmpI;	
						tmpI=tmpI;
					}

					for (int tau=maxI(j,t1); tau<=t-1; tau++) // changed 23-11
					{	// case 8 with ZW tree
						tmpCst=costDSP[i][j][0][tau]+K[i][j]+p[i][j]*(cumD[tau][t2]-U);
						tmpCst+=K[i+1][j]+p[i+1][j]*cumD[tau][t-1]+costZW[i+2][j][tau][t-1];
						if (tmpCst<minCst)
							minCst=tmpCst; 
						if ((checkUsed) && ((costDSP[i][j][0][tau]<-eps)||(costZW[i+2][j][tau][t-1]<-eps)))
						{	cout << "Error: var not set 06" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}
						if (checkUsed && (cumD[tau][t2]-U<-eps))
						{	cout << "Error: negative flow 06" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
					}
				}
				
				// case 9 without ZW tree
				tmpCst=costDSP[i][j-1][1][t]+h[i][j-1]*(cumD[t][t2]-U);
				if ((checkUsed) && (costDSP[i][j-1][1][t]<-eps))
				{	cout << "Error: var not set 05" << endl;
					cin >> tmpI;
					tmpI=tmpI;
				}
				if (tmpCst<minCst)
					minCst=tmpCst;
				
				for (int tau=maxI(j,t1); tau<=t-1; tau++) // changed 23-11
				{	// case 9 with ZW tree
					tmpCst=costDSP[i][j-1][1][tau]+h[i][j-1]*(cumD[tau][t2]-U);
					tmpCst+=K[i+1][j]+p[i+1][j]*cumD[tau][t-1]+costZW[i+2][j][tau][t-1];
					if (tmpCst<minCst)
						minCst=tmpCst;
					if ((checkUsed) && ((costDSP[i][j-1][1][tau]<-eps)||(costZW[i+2][j][tau][t-1]<-eps)))
					{	cout << "Error: var not set 07" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
				}
				
				costDSP[i][j][1][t]=minCst;
			}

			//type 0 arcs												
			for (int t=maxI(t1,j); t<=tUB+1; t++) // t=t1 changed 27-11
			{	minCst=inf;

				if (i>1) // else type 0 arcs that are needed are not defined 27-11-13
				{	// case 7
					tmpCst=costDSP[i][j][0][t]+K[i][j]+p[i][j]*(cumD[t][t2]-U);
					if (tmpCst<minCst)
						minCst=tmpCst;
					if ((checkUsed) && (costDSP[i][j][0][t]<-eps))
					{	cout << "Error: var not set 08" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (checkUsed && (cumD[t][t2]-U<-eps))
					{	cout << "Error: negative flow 07" << endl;
						cin >> tmpI;	
						tmpI=tmpI;
					}
				}
					
				// case 6
				tmpCst=costDSP[i][j-1][1][t]+h[i][j-1]*(cumD[t][t2]-U);
				if (tmpCst<minCst)
					minCst=tmpCst;
				if ((checkUsed) && (costDSP[i][j-1][1][t]<-eps))
				{	cout << "Error: var not set 09" << endl;
					cin >> tmpI;
					tmpI=tmpI;
				}
					
				costDSP[i+1][j][0][t]=minCst;	
			}			
			
		}
	}

	//Final level pre-sink arcs
	// case 8
	/*costDSP[m-1][b1][1][b1+1]=costDSP[m-1][b1][0][b1]+K[m-1][b1]+p[m-1][b1]*(cumD[t1][t2]-U);
	if (costDSP[m-1][b1][0][b1]<-eps)
	{	//cout << "Error: var not set" << endl;
		//cin >> tmpI;
	}*/
		
	for (int j=max(b1+1,t1); j<=tUB; j++) 
	{	minCst=inf;
		
		if (m>2)
		{	// case 8 without ZW tree		
			tmpCst=costDSP[m-1][j][0][j]+K[m-1][j]+p[m-1][j]*(cumD[j][t2]-U);
			if (tmpCst<minCst)
					minCst=tmpCst;
			if ((checkUsed) && (costDSP[m-1][j][0][j]<-eps))
			{	cout << "Error: var not set 10" << endl;
				cin >> tmpI;
				tmpI=tmpI;
			}
			if (checkUsed && (cumD[j][t2]-U<-eps))
			{	cout << "Error: negative flow 08" << endl;
				cin >> tmpI;	
				tmpI=tmpI;
			}
		}
		
		// case 9 without ZW tree
		if (j>t1) //changed 24-11
		{	tmpCst=costDSP[m-1][j-1][1][j]+h[m-1][j-1]*(cumD[j][t2]-U);
			if (tmpCst<minCst)
				minCst=tmpCst;
			if ((checkUsed) && (costDSP[m-1][j-1][1][j]<-eps))
			{	cout << "Error: var not set 11" << endl;
				cin >> tmpI;
				tmpI=tmpI;
			}
		}

		costDSP[m-1][j][1][j+1]=minCst;		
	}

	/*cout << costDSP[1][0][1][1] << endl;

	if ((b1==1) && (b2==2) && (t1==1) && (t2==2))
		cout << "Stop for debugging" << endl;*/

	//Initialize post-sink arcs

	// type 2 arcs at last level
	for (int j=t1+1; j<=tUB+1; j++) // j=b1+1 changed
	{	// case 5 without double ZW tree
		for (int t=tUB+2; t<=t2+1; t++) // t=j+1 changed 26-11-13
		{	costDSP[m-1][j][2][t]=costDSP[m-1][j-1][1][j]+h[m-1][j-1]*(cumD[j][t2]-U)+costZW[m][j][j][t-1];
			if ((checkUsed) && ((costDSP[m-1][j-1][1][j]<-eps)||(costZW[m][j][j][t-1]<-eps)))
			{	cout << "Error: var not set 12" << endl;
				cin >> tmpI;
 				tmpI=tmpI;
			}
			if (checkUsed && (cumD[j][t2]-U<-eps))
			{	cout << "Error: negative flow 09" << endl;
				cin >> tmpI;	
				tmpI=tmpI;
			}
		}
	}

	// type 2 arcs at period tUB+1	
	for (int i=m-2; i>=1; i--)
	{	for (int t=tUB+2; t<=t2+1; t++) // tUB+1 changed 26-11-13
		{	// case 1
			// no ZW tree
			if (t==t2+1)
				minCst=costDSP[i+1][tUB+1][2][t]+K[i+1][tUB+1]+p[i+1][tUB+1]*U;
			else
			{	minCst=costDSP[i+1][tUB+1][2][t]+K[i+1][tUB+1]+p[i+1][tUB+1]*(U-cumD[t][t2]);
				if (checkUsed && (U-cumD[t][t2]<-eps))
				{	cout << "Error: negative inventory 01" << endl;
					cin >> tmpI;	
					tmpI=tmpI;
				}
			}
			if ((checkUsed) && (costDSP[i+1][tUB+1][2][t]<-eps))
			{	cout << "Error: var not set 13" << endl;
				cin >> tmpI;
				tmpI=tmpI;
			}

			//with ZW tree
			for (int tau=tUB+2; tau<=t-1; tau++)
			{	tmpCst=costDSP[i+1][tUB+1][2][tau]+K[i+1][tUB+1]+p[i+1][tUB+1]*(U-cumD[tau][t2]);
				tmpCst+=h[i][tUB+1]*cumD[tau][t-1]+costZW[i+1][tUB+2][tau][t-1];
				if (tmpCst<minCst)
					minCst=tmpCst;
				if ((checkUsed) && ((costDSP[i+1][tUB+1][2][tau]<-eps)||(costZW[i+1][tUB+2][tau][t-1]<-eps)))
				{	cout << "Error: var not set 14" << endl;
					cin >> tmpI;
					tmpI=tmpI;
				}
				if (checkUsed &&(U-cumD[tau][t2]<-eps))
				{	cout << "Error: negative inventory 02" << endl;
					cin >> tmpI;	
					tmpI=tmpI;
				}
			}

			// case 5				
			for (int tau=tUB+1; tau<=tUB+1; tau++)
			{	if (t-1>=tau)
					tmpCst=costDSP[i][tUB][1][tau]+h[i][tUB]*(cumD[tau][t2]-U)+costZW[i+1][tUB+1][tau][t-1];
				else // no ZW tree
					tmpCst=costDSP[i][tUB][1][tau]+h[i][tUB]*(cumD[tau][t2]-U);
				//if ((checkUsed) && ((costDSP[i][tUB][1][tau]<-eps)||(costZW[i+1][tUB+1][tau][t-1])))
				if ((checkUsed) && (costDSP[i][tUB][1][tau]<-eps))
				{	cout << "Error: var not set 15" << endl;
					cin >> tmpI;
					tmpI=tmpI;
				}
				if (checkUsed && (cumD[tau][t2]-U<-eps))
				{	cout << "Error: negative flow 10" << endl;
					cin >> tmpI;	
					tmpI=tmpI;
				}
				
				if (tmpCst<minCst)
					minCst=tmpCst;				
			}

			costDSP[i][tUB+1][2][t]=minCst;

		}
	}

	// type 3 arcs at period tUB+1
	for (int i=2; i<m; i++)
	{	// case 4
		for (int t=tUB+2; t<=t2+1; t++) // tUB+1 changed 26-11-13
		{	// no ZW tree
			if (t==t2+1)
				minCst=costDSP[i][tUB+1][2][t]+K[i][tUB+1]+p[i][tUB+1]*U;
			else
			{	minCst=costDSP[i][tUB+1][2][t]+K[i][tUB+1]+p[i][tUB+1]*(U-cumD[t][t2]);
				if (checkUsed &&(U-cumD[t][t2]<-eps))
				{	cout << "Error: negative inventory 03" << endl;
					cin >> tmpI;	
					tmpI=tmpI;
				}
			}
			if ((checkUsed) && (costDSP[i][tUB+1][2][t]<-eps))
			{	cout << "Error: var not set 16" << endl;
				cin >> tmpI;
				tmpI=tmpI;
			}
			
			//with ZW tree
			for (int tau=tUB+2; tau<=t-1; tau++)
			{	tmpCst=costDSP[i][tUB+1][2][tau]+K[i][tUB+1]+p[i][tUB+1]*(U-cumD[tau][t2]);
				tmpCst+=h[i][tUB+1]*cumD[tau][t-1]+costZW[i][tUB+2][tau][t-1];
				if (tmpCst<minCst)
					minCst=tmpCst;
				if ((checkUsed) && ((costDSP[i][tUB+1][2][tau]<-eps)||(costZW[i][tUB+2][tau][t-1]<-eps)))
				{	cout << "Error: var not set 17" << endl;
					cin >> tmpI;
					tmpI=tmpI;
				}
				if (checkUsed &&(U-cumD[tau][t2]<-eps))
				{	cout << "Error: negative inventory 04" << endl;
					cin >> tmpI;	
					tmpI=tmpI;
				}
			}
			costDSP[i][tUB+1][3][t]=minCst;
		}
	}

	//cout << costDSP[1][1][2][3] << endl;	

	/*if ((b1==0) && (b2==1) && (t1==0) && (t2==2))
		cout << "Stop for debugging" << endl;*/

	//General DP post-sink arcs
	for (int i=m-1; i>=2; i--)
	{	
		//type 3 arcs
		for (int j=tUB; j>=b1+2; j--) // b1+1 changed 26-11-13
		{	for (int t=maxI(j+1,tUB+2); t<=t2+1; t++) // j+1 changed 26-11-13
			{	
				minCst=inf;
				
				if ((i<m-1)||(j>t1)) //else we need a type 2 arc at the final level that not exists 27-11-13 
				{	// case 4 no tree
					if (t==t2+1)
						tmpCst=costDSP[i][j][2][t]+K[i][j]+p[i][j]*U;
					else
					{	tmpCst=costDSP[i][j][2][t]+K[i][j]+p[i][j]*(U-cumD[t][t2]);
						if (checkUsed &&(U-cumD[t][t2]<-eps))
						{	cout << "Error: negative inventory 05" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
					}
					if (tmpCst<minCst)
							minCst=tmpCst;
					if ((checkUsed) && (costDSP[i][j][2][t]<-eps))
					{	cout << "Error: var not set 18" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}

					// case 4 with ZW tree				
					for (int tau=maxI(j+1,tUB+2); tau<=t-1; tau++) // j+1 changed 26-11-13
					{	tmpCst=costDSP[i][j][2][tau]+K[i][j]+p[i][j]*(U-cumD[tau][t2]);
						tmpCst+=h[i][j]*cumD[tau][t-1]+costZW[i][j+1][tau][t-1];
						if (tmpCst<minCst)
							minCst=tmpCst;
						if ((checkUsed) && ((costDSP[i][j][2][tau]<-eps)||(costZW[i][j+1][tau][t-1]<-eps)))
						{	cout << "Error: var not set 19" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}
						if (checkUsed &&(U-cumD[tau][t2]<-eps))
						{	cout << "Error: negative inventory 06" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
					}
				}

				// case 3
				if (t>j+1) 
				{	if (t==t2+1)
						tmpCst=costDSP[i][j+1][3][t]+h[i][j]*U;
					else
					{	tmpCst=costDSP[i][j+1][3][t]+h[i][j]*(U-cumD[t][t2]);
						if (checkUsed &&(U-cumD[t][t2]<-eps))
						{	cout << "Error: negative inventory 07" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
					}
					if ((checkUsed) && (costDSP[i][j+1][3][t]<-eps))
					{	cout << "Error: var not set 20" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (tmpCst<minCst)
						minCst=tmpCst;
				}

				costDSP[i][j][3][t]=minCst;
			}
		}

		// type 2 arcs
		for (int j=tUB; j>=b1+1; j--)
		{	for (int t=maxI(j+1,tUB+2); t<=t2+1; t++) // maxI(j+1,t1) changed 26-11-13
			{	//case 2
				minCst=inf;
				if (t>j+1) // type 3 arc satisfies at least demands up till period j+1
				{	if (t==t2+1)
						tmpCst=costDSP[i][j+1][3][t]+h[i-1][j]*U;
					else
					{	tmpCst=costDSP[i][j+1][3][t]+h[i-1][j]*(U-cumD[t][t2]);
						if (checkUsed &&(U-cumD[t][t2]<-eps))
						{	cout << "Error: negative inventory 08" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
					}
					if ((checkUsed) && (costDSP[i][j+1][3][t]<-eps))
					{	cout << "Error: var not set 21" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (tmpCst<minCst)
						minCst=tmpCst;
				}

				// case 1
				// no ZW tree
				if ((i<m-1)||(j>t1)) //else we need a type 2 arc at the final level that not exists 27-11-13 
				{	if (t==t2+1)
						tmpCst=costDSP[i][j][2][t]+K[i][j]+p[i][j]*U;
					else
					{	tmpCst=costDSP[i][j][2][t]+K[i][j]+p[i][j]*(U-cumD[t][t2]);
						if (checkUsed &&(U-cumD[t][t2]<-eps))
						{	cout << "Error: negative inventory 09" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
					}
					if ((checkUsed) && (costDSP[i][j][2][t]<-eps))
					{	cout << "Error: var not set 22" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (tmpCst<minCst)
						minCst=tmpCst;				

					//with ZW tree
					for (int tau=maxI(j+1,tUB+2); tau<=t-1; tau++) // changed 26-11-13
					{	tmpCst=costDSP[i][j][2][tau]+K[i][j]+p[i][j]*(U-cumD[tau][t2]);
						tmpCst+=h[i][j]*cumD[tau][t-1]+costZW[i][j+1][tau][t-1];
						if ((checkUsed) && ((costDSP[i][j][2][tau]<-eps)||(costZW[i][j+1][tau][t-1]<-eps)))
						{	cout << "Error: var not set 23" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}
						if (checkUsed &&(U-cumD[tau][t2]<-eps))
						{	cout << "Error: negative inventory 10" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
						if (tmpCst<minCst)
							minCst=tmpCst;
					}
				}

				// case 5				
				for (int tau=maxI(j,t1); tau<=minI(t-1,tUB+1); tau++)
				{	tmpCst=costDSP[i-1][j-1][1][tau]+h[i-1][j-1]*(cumD[tau][t2]-U)+costZW[i][j][tau][t-1];
					if ((checkUsed) && ((costDSP[i-1][j-1][1][tau]<-eps)||(costZW[i][j][tau][t-1]<-eps)))
					{	cout << "Error: var not set 24" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (checkUsed && (cumD[tau][t2]-U<-eps))
					{	cout << "Error: negative flow 11" << endl;
						cin >> tmpI;	
						tmpI=tmpI;
					}
					if (tmpCst<minCst)
						minCst=tmpCst;
				}

				costDSP[i-1][j][2][t]=minCst;
			}
		}
	}
	
	/*cout << costDSP[1][1][2][2] << endl;	

	if ((b1==0) && (b2==1) && (t1==0) && (t2==1))
		cout << "Stop for debugging" << endl;*/


	//First level post-sink arcs
	for (int b2=b1+1; b2<=t2; b2++)
	{	minCst=inf;
		for (int j=b1+1; (j<=tUB+1)&&(j<=b2); j++)	
		{	// no BN tree
			tmpCst=costDSP[1][j][2][t2+1]+cumHU[b1][j]+K[1][j]+p[1][j]*U;
			if ((checkUsed) && (costDSP[1][j][2][t2+1]<-eps))
			{	cout << "Error: var not set 25" << endl;
				cin >> tmpI;
				tmpI=tmpI;
			}
			if (tmpCst<minCst)
				minCst=tmpCst;

			// with BN tree
			if (j<b2) //otherwise no tree possible
			{	for (int t=j+1; t<=t2; t++)
				{	if (cumD[t][t2]<=U+eps)
					{	tmpCst=costDSP[1][j][2][t]+cumHU[b1][j]+K[1][j]+p[1][j]*(U-cumD[t][t2]);
						tmpCst+=h[0][j]*cumD[t][t2]+costBN[j+1][b2][t][t2];
						if ((checkUsed) && ((costDSP[1][j][2][t]<-eps)||(costBN[j+1][b2][t][t2]<-eps)))
						{	cout << "Error: var not set 26" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}
						if (checkUsed &&(U-cumD[t][t2]<-eps))
						{	cout << "Error: negative inventory 11" << endl;
							cin >> tmpI;	
							tmpI=tmpI;
						}
						if (tmpCst<minCst)
							minCst=tmpCst;
					}
				}
			}
		}
		costDS[b1][b2][t1][t2]=minCst;
		//cout << b1 << ' ' << b2 << ' ' << t1 << ' ' << t2 << ": " << minCst << endl;
	}

	}
	else
	{	for (int b2=b1+1; b2<=t2; b2++)
			costDS[b1][b2][t1][t2]=inf;
	}

	DeleteDim2Dbl(cumD,n,n);
	DeleteDim2Dbl(cumHU,n,n);
	DeleteDim4Dbl(costDSP,m,n,nDir,n+1);		
}

void BottleneckTrees(int m, int n, double *d, double **K, double **p, double **h, double U, double ****costZW, double ****costBN)
{	// define 4-dim cost array
	double inf, eps, **cumD, tmpCst, minCst, tmpDbl;	
	bool checkUsed=false;
	SetDim2Dbl(cumD,n,n);	// cumulative demands	
	
	// precomputation
	inf=pow(10.0,10.0);
	eps=pow(10.0,-5.0);
	cumDemands(n,d,cumD);
	
	// cout << costZW[0][0][0][n-1] << endl;

	// Initialization
	for (int b=0;b<n;b++)
	{	for (int t1=b;t1<n;t1++)	
		{	for (int t2=t1;t2<n;t2++)
			{	costBN[b][b][t1][t2]=K[1][b]+p[1][b]*cumD[t1][t2]+costZW[2][b][t1][t2];				
				if ((checkUsed) && (costZW[2][b][t1][t2]<-eps))
				{	cout << "Error: variable not set" << endl;
					cin >> tmpDbl;
				}
				/*cout << K[1][b] << ' ' << p[1][b] << ' ' << cumD[t1][t2] << endl;
				PrintCost(2,b,t1,t2,costZW[2][b][t1][t2],"ZW");
				PrintCost(b,b,t1,t2,costBN[b][b][t1][t2],"BN");*/
			}
		}
	}

	// General DP
	for (int t2=0; t2<n; t2++)
	{	for (int b2=0; b2<=t2; b2++)	
		{	for (int b1=b2-1; b1>=0; b1--)
			{	for(int t1=t2; t1>=b1; t1--)
				{	minCst=inf;

					// for debugging purposes
					/*if ((b1==0)&&(b2==2)&&(t1==0)&&(t2==2))
						b1=b1;
					*/
					
					// demand is split in operational and inventory flow at (a,b1)
					//for (int t=t1; (t<=t2-1)&&(cumD[t+1][t2]<=U+eps); t++)
					for (int t=t2-1; (t>=t1)&&(cumD[t+1][t2]<=U+eps); t--)
					{	tmpCst=K[1][b1]+p[1][b1]*cumD[t1][t]+costZW[2][b1][t1][t];
						tmpCst+=h[0][b1]*cumD[t+1][t2]+costBN[b1+1][b2][t+1][t2];
						if (tmpCst<minCst)
							minCst=tmpCst;
						if ((checkUsed) && (costZW[2][b1][t1][t]<-eps))
						{	cout << "Error: variable not set" << endl;
							cin >> tmpDbl;
						}
						if ((checkUsed) && (costBN[b1+1][b2][t+1][t2]<-eps))
						{	cout << "Error: variable not set" << endl;						
							cin >> tmpDbl;
						}
					}

					// inventory flow only at (a,b1) if possible
					if ((t1>b1)&&(cumD[t1][t2]<=U+eps))
					{	tmpCst=h[0][b1]*cumD[t1][t2]+costBN[b1+1][b2][t1][t2];
						if (tmpCst<minCst)
							minCst=tmpCst;
						if ((checkUsed) && (costBN[b1+1][b2][t1][t2]<-eps))
						{	cout << "Error: variable not set" << endl;
							cin >> tmpDbl;
						}
					}

					// set min cost
					costBN[b1][b2][t1][t2]=minCst;
					//PrintCost(b1,b2,t1,t2,costBN[b1][b2][t1][t2],"BN");
				}
			}
		}
	}	

	// compare ZW and BN trees for errors
	/*bool check;
	for (int b1=0; b1<n; b1++)
	{	for (int t1=b1; t1<n; t1++)	
		{	for (int t2=t1; t2<n; t2++)
			{	check=false;	
				for(int b2=b1; b2<=t2; b2++)
				{	if (costBN[b1][b2][t1][t2]<costZW[1][b1][t1][t2]-eps)
					{	cout << "Error: ZW vs BN" << endl;
						cin >> tmpDbl;
					}
					if (abs(costBN[b1][b2][t1][t2]-costZW[1][b1][t1][t2])<eps)
						check=true;
				}
				if (check==false)
				{	cout << "Error: ZW vs BN" << endl;
					cin >> tmpDbl;
				}
			}
		}
	}*/	
	
	DeleteDim2Dbl(cumD,n,n);
}

void costZWtrees(int m, int n, double *d, double **K, double **p, double **h, double ****cost)
{	// define variables	
	double tmpCost, minCost, opt, eps=0.00001;	
	double **cumD, **cumH;
	bool checkUsed=false;
	int tmpI;
	SetDim2Dbl(cumD,n,n);	// cumulative demands
	SetDim2Dbl(cumH,n,n);	// cumulative holding cost
	cumDemands(n,d,cumD);
	cumHolding(m,n,cumD,h,cumH);	

	// initialize cost at retailer level
	for (int t1=0;t1<n;t1++)	
	{	for (int t2=t1;t2<n;t2++)
		{	cost[m][t1][t1][t2]=cumH[t1][t2];
			//PrintCost(m,t1,t1,t2,cost[m][t1][t1][t2],"ZW");
			for (int j=t1-1;j>=0;j--)
			{	cost[m][j][t1][t2]=h[m-1][j]*cumD[t1][t2]+cost[m][j+1][t1][t2];
				//PrintCost(m,j,t1,t2,cost[m][j][t1][t2],"ZW");				
				if ((checkUsed) && (cost[m][j+1][t1][t2]<-eps))
				{	cout << "Error: var not set 27" << endl;
					cin >> tmpI;
					tmpI=tmpI;
				}
			}
		}
	}

	// initialize cost at last period
	for (int i=m-1;i>=0;i--)	
	{	cost[i][n-1][n-1][n-1]=cost[i+1][n-1][n-1][n-1]+K[i][n-1]+p[i][n-1]*d[n-1];
		//PrintCost(i,n-1,n-1,n-1,cost[i][n-1][n-1][n-1],"ZW");
		if ((checkUsed) && (cost[i+1][n-1][n-1][n-1]<-eps))
		{	cout << "Error: var not set 28" << endl;
			cin >> tmpI;
			tmpI=tmpI;
		}
	}

	// apply general DP
	for (int i=m-1;i>=0;i--)
	{	for (int j=n-2;j>=0;j--)
		{	for (int t1=j;t1<n;t1++)	
			{	for (int t2=t1;t2<n;t2++)
				{	// special case: production flow only 
					minCost=K[i][j]+p[i][j]*cumD[t1][t2]+cost[i+1][j][t1][t2];
					if ((checkUsed) && (cost[i+1][j][t1][t2]<-eps))
					{	cout << "Error: var not set 29" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
										
					// special case: inventory flow only
					if (t1>j) // otherwise not possible
					{	if (i==0) //dummy level
							tmpCost=cost[i][j+1][t1][t2];
						else
						{	tmpCost=h[i-1][j]*cumD[t1][t2]+cost[i][j+1][t1][t2];							
						}
						if ((checkUsed) && (cost[i][j+1][t1][t2]<-eps))
						{	cout << "Error: var not set 30" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}

						if (tmpCost<minCost)
							minCost=tmpCost;
					}

					// production and inventory flow
					for (int t=t1; t<t2; t++)
					{	tmpCost=K[i][j]+p[i][j]*cumD[t1][t]+cost[i+1][j][t1][t];
						if ((checkUsed) && (cost[i+1][j][t1][t]<-eps))
						{	cout << "Error: var not set 31" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}

						if (i==0) //dummy level
							tmpCost+=cost[i][j+1][t+1][t2];
						else
							tmpCost+=h[i-1][j]*cumD[t+1][t2]+cost[i][j+1][t+1][t2];
						if ((checkUsed) && (cost[i][j+1][t+1][t2]<-eps))
						{	cout << "Error: var not set 32" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}

						if (tmpCost<minCost)
							minCost=tmpCost;
					}

					if (minCost<-eps)
					{	cout << "Error: in DP" << endl;
						cin >> tmpCost;
					}
					cost[i][j][t1][t2]=minCost;
					//PrintCost(i,j,t1,t2,cost[i][j][t1][t2],"ZW");

				}
			}				
		}
	}

	opt=cost[0][0][0][n-1];

	DeleteDim2Dbl(cumD,n,n);
	DeleteDim2Dbl(cumH,n,n);	
}

void UpstreamTrees(int L, int n, double *d, double **K, double **p, double **h, double ****costZW, double ****costBN, double ****costDS, double ****costUS)
{	double **cumD, ******cost, tmpCst, minCst, eps=0.00001, inf=pow(10.0,10.0);	
	bool checkUsed=true;
	bool k=0; //only use two SC levels in DP
	int tmpI;
	SetDim2Dbl(cumD,n,n);
	SetDim6Dbl(cost,2,n,n,n,n,n);
	if (checkUsed)
		SetConstDim6Dbl(cost,2,n,n,n,n,n,-1.0);
	cumDemands(n,d,cumD);

	//Initialization
	//only a ZW tree
	for (int t2=0; t2<n; t2++)
	{	for (int t1=0; t1<=t2; t1++)
		{	for (int j=0; j<=t1; j++)
			{	//cost[L+1][j][j][j][t1][t2]=K[L+1][j]+p[L+1][j]*cumD[t1][t2]+costZW[2][j][t1][t2];
				cost[0][j][j][j][t1][t2]=K[L+1][j]+p[L+1][j]*cumD[t1][t2]+costZW[2][j][t1][t2];
				if ((checkUsed) && (costZW[2][j][t1][t2]<-eps))
				{	cout << "Error: var not set US 01" << endl;
					cin >> tmpI;
					tmpI=tmpI;
				}
			}
		}
	}

	//intervals at BN level through 1 period
	for (int t2=0; t2<n; t2++)
	{	for (int t1=0; t1<=t2; t1++)
		{	for (int b1=0; b1<=t1; b1++)
			{	for (int b2=b1+1; b2<=t2; b2++)
				{	minCst=inf;
					
					// add this code if at most up till b2 is used
					//tmpCst=cost[L+1][b1][b1][b2-1][t1][t2];
					tmpCst=cost[0][b1][b1][b2-1][t1][t2];
					if ((checkUsed) && (cost[0][b1][b1][b2-1][t1][t2]<-eps))
					{	cout << "Error: var not set US 02" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (tmpCst<minCst)
						minCst=tmpCst;

					//NB: cost is already counted in DS path (double counting is avoided)
					tmpCst=costDS[b1][b2][t1][t2]-K[L][b1]-p[L][b1]*cumD[t1][t2];
					if ((checkUsed) && (costDS[b1][b2][t1][t2]<-eps))
					{	cout << "Error: var not set US 03" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (tmpCst<minCst)
						minCst=tmpCst;

					// BN tree may be better
					tmpCst=costBN[b1][b2][t1][t2];
					if ((checkUsed) && (costBN[b1][b2][t1][t2]<-eps))
					{	cout << "Error: var not set US 03" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					if (tmpCst<minCst)
						minCst=tmpCst;

					cost[0][b1][b1][b2][t1][t2]=minCst;
				}	
			}
		}
	}

	// initialize cost at last period
	/*for (int i=L;i>=0;i--)	
	{	cost[i][n-1][n-1][n-1][n-1][n-1]=cost[i+1][n-1][n-1][n-1][n-1][n-1]+K[i][n-1]+p[i][n-1]*d[n-1];
		if ((checkUsed) && (cost[i+1][n-1][n-1][n-1][n-1][n-1]<-eps))
		{	cout << "Error: var not set US 04" << endl;
			cin >> tmpI;
			tmpI=tmpI;
		}
	}*/

	//General DP
	for (int i=L;i>=0;i--)
	{	k=1-k; //change row index
		
		cost[k][n-1][n-1][n-1][n-1][n-1]=cost[1-k][n-1][n-1][n-1][n-1][n-1]+K[i][n-1]+p[i][n-1]*d[n-1];
		if ((checkUsed) && (cost[1-k][n-1][n-1][n-1][n-1][n-1]<-eps))
		{	cout << "Error: var not set US 04" << endl;
			cin >> tmpI;
			tmpI=tmpI;
		}
		
		for (int j=n-2;j>=0;j--)
		{	for (int b1=j; b1<n; b1++)
			{	for (int b2=b1; b2<n; b2++)
				{	for (int t1=b1;t1<n;t1++)	
					{	for (int t2=maxI(t1,b2);t2<n;t2++)
						{	
							/*if ((i==1)&&(j==0)&&(b1==0)&&(b2==2)&&(t1==0)&&(t2==2))
							{	cout << "Stop for debugging" << endl;
							}*/

							minCst=inf;

							// special case: production flow only 
							if ((i<L)||(j==b1))
							{	tmpCst=K[i][j]+p[i][j]*cumD[t1][t2]+cost[1-k][j][b1][b2][t1][t2];
								if (tmpCst<minCst)
									minCst=tmpCst;
								if ((checkUsed) && (cost[1-k][j][b1][b2][t1][t2]<-eps))
								{	cout << "Error: var not set US 05" << endl;
									cin >> tmpI;
									tmpI=tmpI;
								}
							}
							
							// special case: inventory flow only
							if (b1>j) // otherwise not possible
							{	if (i==0) //dummy level
									tmpCst=cost[k][j+1][b1][b2][t1][t2];
								else
									tmpCst=h[i-1][j]*cumD[t1][t2]+cost[k][j+1][b1][b2][t1][t2];

								if ((checkUsed) && (cost[k][j+1][b1][b2][t1][t2]<-eps))
								{	cout << "Error: var not set US 06" << endl;
									cin >> tmpI;
									tmpI=tmpI;
								}

								if (tmpCst<minCst)
									minCst=tmpCst;
							}

							
							// production and inventory flow
							if ((i<L)||(j==b1))
							{	for (int b=b1; b<=b2-1; b++)								
								{	for (int t=maxI(b,t1); t<=t2-1; t++)
									{	tmpCst=K[i][j]+p[i][j]*cumD[t1][t]+cost[1-k][j][b1][b][t1][t];
										if ((checkUsed) && (cost[1-k][j][b1][b][t1][t]<-eps))
										{	cout << "Error: var not set US 07" << endl;
											cin >> tmpI;
											tmpI=tmpI;
										}

										if (i==0) //dummy level
											tmpCst+=cost[k][j+1][b+1][b2][t+1][t2];
										else
											tmpCst+=h[i-1][j]*cumD[t+1][t2]+cost[k][j+1][b+1][b2][t+1][t2];

										if ((checkUsed) && (cost[k][j+1][b+1][b2][t+1][t2]<-eps))
										{	cout << "Error: var not set US 08" << endl;
											cin >> tmpI;
											tmpI=tmpI;
										}

										if (tmpCst<minCst)
											minCst=tmpCst;
									}
								}
							}

							if (minCst<-0.00001)
							{	cout << "Error: in US tree DP" << endl;
								cin >> tmpCst;
							}

							cost[k][j][b1][b2][t1][t2]=minCst;

						}
					}
				}
			}
		}
	}

	// find min over all possible sources	
	for (int b1=0; b1<n; b1++)
	{	for (int t1=b1; t1<n; t1++)
		{	for (int t2=t1; t2<n; t2++)
			{	for (int b2=b1; b2<=t2; b2++)
				{	minCst=cost[k][0][b1][b2][t1][t2];
					if ((checkUsed) && (cost[k][0][b1][b2][t1][t2]<-eps))
					{	cout << "Error: var not set US 09" << endl;
						cin >> tmpI;
						tmpI=tmpI;
					}
					for (int j=1; j<=b1; j++)
					{	tmpCst=cost[k][j][b1][b2][t1][t2];
						if ((checkUsed) && (cost[k][j][b1][b2][t1][t2]<-eps))
						{	cout << "Error: var not set US 10" << endl;
							cin >> tmpI;
							tmpI=tmpI;
						}

						if (tmpCst<minCst)
							minCst=tmpCst;
					}
					costUS[b1][b2][t1][t2]=minCst;
				}
			}
		}
	}

	DeleteDim6Dbl(cost,2,n,n,n,n,n);
	DeleteDim2Dbl(cumD,n,n);
}

void SetConstDim4Dbl(double ****x, int n1, int n2, int n3, int n4, double c)
{	for (int i1=0; i1<n1; i1++)
	{	for (int i2=0; i2<n2; i2++)
		{	for (int i3=0; i3<n3; i3++)
			{	for (int i4=0; i4<n4; i4++)
					x[i1][i2][i3][i4]=c;
			}
		}
	}

}

void SetConstDim6Dbl(double ******x, int n1, int n2, int n3, int n4, int n5, int n6, double c)
{	for (int i1=0; i1<n1; i1++)
	{	for (int i2=0; i2<n2; i2++)
		{	for (int i3=0; i3<n3; i3++)
			{	for (int i4=0; i4<n4; i4++)
				{	for (int i5=0; i5<n5; i5++)
					{	for (int i6=0; i6<n6; i6++)
						{	x[i1][i2][i3][i4][i5][i6]=c;
						}
					}
				}
			}
		}
	}
}

void cumHoldU(int n, double U, double **h, double **cumH)
{	for (int t1=0; t1<n-1; t1++)	
	{	cumH[t1][t1+1]=h[0][t1]*U;
		for (int t2=t1+2; t2<n; t2++)
		{	cumH[t1][t2]=cumH[t1][t2-1]+h[0][t2-1]*U;
		}
	}
}

int maxI(int x, int y)
{	if (x>y)
		return x;
	else
		return y;
}

int minI(int x, int y)
{	if (x<y)
		return x;
	else
		return y;
}

void PrintCost(int n1, int n2, int n3, int n4, double cst, char* str)
{	cout << str << ":" << ' ' << n1 << ' ' << n2 << ' ' << n3 << ' ' << n4 << ": " << cst << endl;
}

void CopyPart(int m1, int m2, int n, double **x, double **y)
{	for (int i=0; i<m1; i++)
	{	for (int j=0; j<n; j++)
		{	y[i][j]=x[m2-m1+i][j];
		}
	}
}

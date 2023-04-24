#include <iostream>
#include <cmath>
#include "arrays.h"

using namespace std;

void cumDemands(int n, double *demand, double **matrix);
void cumHolding(int m, int n, double **cumD, double **h, double **cumH);
int max(int m, int n);

double solveZangwill(int m, int n, double *d, double **consCost, double **varCost, double **invCost)
{	double inf=pow(10.0,10.0);

	// compute cumuative demands
	double **demands;
	SetDim2Dbl(demands,n,n);
	cumDemands(n,d,demands);

	// define 4-dim cost
	double ****costs; //costs[a][b][c][d] is the minimal cost of making production for period c up til d from point a,b 
	SetDim4Dbl(costs,m,n,n,n);

	for (int level=m-1;level>=0;level--)	
	{	for (int time=n-1;time>=0;time--)
		{	for (int beginperiod=n-1;beginperiod>=0;beginperiod--)
			{	for (int endperiod=n-1;endperiod>=0;endperiod--)
				{	if (endperiod<beginperiod || time>beginperiod)// we have a negative amount to give so costs=0.
					{	costs[level][time][beginperiod][endperiod]=0;
					}
					else if (level==m-1) //j=m
					{	if (time==n-1) // j=m and i=n
						{	costs[level][time][beginperiod][endperiod]=0;
						}
						else //j=m and i<n
						{ 	if (beginperiod+1>endperiod)
							{	costs[level][time][beginperiod][endperiod]=0;
							}
							else
							{	costs[level][time][time][endperiod]=invCost[level][time]*demands[time+1][endperiod]+costs[level][time+1][time+1][endperiod];
							}
						}						
					}
					else if (level==m-2) //j=m-1
					{	if (time==n-1) //j=m-1, i=n
						{	if (demands[beginperiod][endperiod]>0)
							{	costs[level][time][beginperiod][endperiod]=consCost[level][time]+varCost[level][time]*demands[beginperiod][endperiod];
							}
							else
							{	costs[level][time][beginperiod][endperiod]=0;
							}
						}
						else // j=m-1, i<n
						{	costs[level][time][beginperiod][endperiod]=inf;
							for (int i=time;i<=endperiod;i++)
							{	double tempcost=0;
								if (beginperiod>endperiod || time>=n-2)
								{	tempcost=0;
								}
								else
								{	if (demands[beginperiod][endperiod]>0)
									{	if(i+1>endperiod)
										{	tempcost=consCost[level][time]+varCost[level][time]*demands[beginperiod][i]+costs[level+1][time][beginperiod][i];
										}
										else
										{	tempcost=consCost[level][time]+varCost[level][time]*demands[beginperiod][i]+costs[level+1][time][beginperiod][i]+varCost[level][time]*demands[i+1][endperiod]+costs[level][time+1][i+1][endperiod];
										}
									}
									else
									{	tempcost=costs[level+1][time][beginperiod][i]+varCost[level][time]*demands[i+1][endperiod]+costs[level][time+1][i+1][endperiod];
									}
								}
								if (tempcost<costs[level][time][beginperiod][endperiod])
								{	costs[level][time][beginperiod][endperiod]=tempcost;
								}
							}
						}
					}
					else
					{	if(time==n-1)
						{	if (demands[beginperiod][endperiod]>0)
							{	costs[level][time][beginperiod][endperiod]=consCost[level][time]+varCost[level][time]*demands[beginperiod][endperiod];
							}
							else
							{	costs[level][time][beginperiod][endperiod]=costs[level+1][time][beginperiod][endperiod];
							}
						}
						else
						{	costs[level][time][beginperiod][endperiod]=inf;
							if (max(time,beginperiod-1)>endperiod-1)
							{	if (demands[beginperiod][endperiod]>0)
								{	costs[level][time][beginperiod][endperiod]=consCost[level][time]+varCost[level][time]*demands[beginperiod][endperiod]+costs[level+1][time][beginperiod][endperiod]+varCost[level][time]*demands[endperiod][endperiod]+costs[level][time+1][endperiod][endperiod];
								}
								else
								{	costs[level][time][beginperiod][endperiod]=costs[level+1][time][beginperiod][endperiod]+varCost[level][time]*demands[endperiod][endperiod]+costs[level][time+1][endperiod][endperiod];
								}
							}
							else
							{	for (int i=max(time,beginperiod-1);i<=endperiod-1;i++)
								{	double tempcost=0;
									if (beginperiod>endperiod || time>=n)
									{	tempcost=0;
									}
									else
									{	if (demands[beginperiod][endperiod]>0)
										{	tempcost=consCost[level][time]+varCost[level][time]*demands[beginperiod][i]+costs[level+1][time][beginperiod][i]+varCost[level][time]*demands[i+1][endperiod]+costs[level][time+1][i+1][endperiod];
										}
										else
										{	tempcost=costs[level+1][time][beginperiod][i]+varCost[level][time]*demands[i+1][endperiod]+costs[level][time+1][i+1][endperiod];
										}
									}
									if (tempcost<costs[level][time][beginperiod][endperiod])
									{	costs[level][time][beginperiod][endperiod]=tempcost;
									}
								}
							}
						}
					}
				}
			}
		}
	}	
	
	return costs[0][0][0][n-1];

	DeleteDim2Dbl(demands,n,n);
	DeleteDim4Dbl(costs,m,n,n,n);	
}

double solveZangwill2(int m, int n, double *d, double **K, double **p, double **h)
{	// define variables	
	double tmpCost, minCost;	
	double **cumD, **cumH;
	SetDim2Dbl(cumD,n,n);	// cumulative demands
	SetDim2Dbl(cumH,n,n);	// cumulative holding cost
	cumDemands(n,d,cumD);
	cumHolding(m,n,cumD,h,cumH);

	// define 4-dim cost array
	double ****cost; //costs[i][j][t1][t2] min cost to satisfy demands in [t1,t2] from node (i,j)
	SetDim4Dbl(cost,m+1,n,n,n); // add dummy level 0 for convenience

	// initialize cost at retailer level
	for (int t1=0;t1<n;t1++)	
	{	for (int t2=t1;t2<n;t2++)
		{	cost[m][t1][t1][t2]=cumH[t1][t2];
			for (int j=t1-1;j>=0;j--)
			{	cost[m][j][t1][t2]=h[m-1][j]*cumD[t1][t2]+cost[m][j+1][t1][t2];
			}
		}
	}

	// initialize cost at last period
	for (int i=m-1;i>=0;i--)	
	{	cost[i][n-1][n-1][n-1]=cost[i+1][n-1][n-1][n-1]+K[i][n-1]+p[i][n-1]*d[n-1];
	}

	// apply general DP
	for (int i=m-1;i>=0;i--)
	{	for (int j=n-2;j>=0;j--)
		{	for (int t1=j;t1<n;t1++)	
			{	for (int t2=t1;t2<n;t2++)
				{	// special case: production flow only 
					minCost=K[i][j]+p[i][j]*cumD[t1][t2]+cost[i+1][j][t1][t2];
					
					// special case: inventory flow only
					if (t1>j) // otherwise not possible
					{	if (i==0) //dummy level
							tmpCost=cost[i][j+1][t1][t2];
						else
							tmpCost=h[i-1][j]*cumD[t1][t2]+cost[i][j+1][t1][t2];

						if (tmpCost<minCost)
							minCost=tmpCost;
					}

					// production and inventory flow
					for (int t=t1; t<t2; t++)
					{	tmpCost=K[i][j]+p[i][j]*cumD[t1][t]+cost[i+1][j][t1][t];
						if (i==0) //dummy level
							tmpCost+=cost[i][j+1][t+1][t2];
						else
							tmpCost+=h[i-1][j]*cumD[t+1][t2]+cost[i][j+1][t+1][t2];

						if (tmpCost<minCost)
							minCost=tmpCost;
					}

					if (minCost<0)
					{	cout << "Error: in DP" << endl;
						cin >> tmpCost;
					}
					cost[i][j][t1][t2]=minCost;

				}
			}				
		}
	}

	tmpCost=cost[0][0][0][n-1];

	DeleteDim2Dbl(cumD,n,n);
	DeleteDim2Dbl(cumH,n,n);
	DeleteDim4Dbl(cost,m+1,n,n,n);	

	return tmpCost;
}

void cumDemands(int n, double *demand, double **matrix)
{	// we create a matrix with the demands where matrix[a][b] represents the sum of demands from a to b. Note: if a>b matrix[a][b]=0
	for (int i=0;i<n;i++)	
	{	for (int j=0;j<n;j++)
		{	matrix[i][j]=0.0;
			if (i<=j)
			{	for (int h=i;h<=j;h++)
				{	matrix[i][j]=matrix[i][j]+demand[h];
				}		
			}
		}
	}
}

int max(int m, int n)
{	if (m>n)
		return m;
	else
		return n;
}

void cumHolding(int m, int n, double **cumD, double **h, double **cumH)
{	for (int t2=0;t2<n;t2++)	
	{	cumH[t2][t2]=0;
		for (int t1=t2-1;t1>=0;t1--)
		{	cumH[t1][t2]=h[m-1][t1]*cumD[t1+1][t2]+cumH[t1+1][t2];
		}
	}
}

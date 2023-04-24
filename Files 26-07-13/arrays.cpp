void SetDim2Dbl(double **(&x), int n1, int n2)
{	x=new double*[n1];
	for (int i1=0; i1<n1; i1++)
		x[i1]=new double[n2];

}

void SetDim3Dbl(double ***(&x), int n1, int n2, int n3)
{	x=new double**[n1];
	for (int i1=0; i1<n1; i1++)
	{	x[i1]=new double*[n2];
		for (int i2=0; i2<n2; i2++)
			x[i1][i2]=new double[n3];		
	}

}


void SetDim4Dbl(double ****(&x), int n1, int n2, int n3, int n4)
{	x=new double***[n1];
	for (int i1=0; i1<n1; i1++)
	{	x[i1]=new double**[n2];
		for (int i2=0; i2<n2; i2++)
		{	x[i1][i2]=new double*[n3];
			for (int i3=0; i3<n3; i3++)
				x[i1][i2][i3]=new double[n4];			
		}
	}

}

void DeleteDim2Dbl(double **(&x), int n1, int n2)
{	for (int i1=0; i1<n1; i1++)
	{	delete[] x[i1];			
	}	
	delete[] x;
}

void DeleteDim3Dbl(double ***(&x), int n1, int n2, int n3)
{	for (int i1=0; i1<n1; i1++)
	{	for (int i2=0; i2<n2; i2++)
		{	delete[] x[i1][i2];			
		}
		delete[] x[i1];			
	}	
	delete[] x;

}

void DeleteDim4Dbl(double ****(&x), int n1, int n2, int n3, int n4)
{	for (int i1=0; i1<n1; i1++)
	{	for (int i2=0; i2<n2; i2++)
		{	for (int i3=0; i3<n3; i3++)
			{	delete[] x[i1][i2][i3];			
			}
			delete[] x[i1][i2];			
		}
		delete[] x[i1];			
	}	
	delete[] x;
}
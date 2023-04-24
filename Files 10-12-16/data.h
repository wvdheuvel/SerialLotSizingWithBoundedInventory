#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "arrays.h"
#include "randomc.h"

#ifndef Data_H
#define Data_H

using namespace std;

// Data: objects where information from the data file are stored
class Data
{
public:
	// Constructors
	Data(int m, int n); // used when data is created randomly
	Data(char *file); // used when data is read from file
	// Destructor
	~Data();
	
	// Data manipulation
	void DataRnd(TRandomMersenne &ran); // create random data
	void WriteData(); // write instance to file
	void PrintData(); // print instance on screen

	// Copy the instance for easy of use
	void CopyData(int & _L, double *_d, double **_K, double **_p, double **_h, double & _U);

	// Instance parameters
	int m; // number of stages/levels
	int n; // number of periods
	int L; // bottleneck level
	double U; // inventory capacity
	double* d; // vector of demands
	double** K; // matrix of setup costs
	double** p; // matrix of unitary production costs
	double** h; // matrix of unitary inventory costs
};

#endif
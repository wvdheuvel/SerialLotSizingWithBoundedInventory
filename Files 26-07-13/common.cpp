#include "definitions.h"

ILOSTLBEGIN

void SetName(IloExtractable obj, char* prefix, int i)
{
	char buffer[50];
	sprintf(buffer, "_%u", i);

	char name[100];
	strcpy(name, prefix);
	strcat(name, buffer);
	obj.setName(name);
}

void SetName2(IloExtractable obj, char* prefix, int i, int j)
{
	char buffer[50];
	sprintf(buffer, "_%u_%u", i, j);

	char name[100];
	strcpy(name, prefix);
	strcat(name, buffer);
	obj.setName(name);
}

void SetName3(IloExtractable obj, char* prefix, int i, int j, int k)
{
	char buffer[50];
	sprintf(buffer, "_%u_%u_%u", i, j, k);

	char name[100];
	strcpy(name, prefix);
	strcat(name, buffer);
	obj.setName(name);
}

void SetName4(IloExtractable obj, char* prefix, int i, int j, int k, int l)
{
	char buffer[50];
	sprintf(buffer, "_%u_%u_%u_%u", i, j, k, l);

	char name[100];
	strcpy(name, prefix);
	strcat(name, buffer);
	obj.setName(name);
}

IloBoolVarArray CreateBoolVarArray(IloEnv env, IloInt m, char* prefix)
{
	IloBoolVarArray array = IloBoolVarArray(env, m);
	for (int i = 0; i < m; i++)
		SetName(array[i], prefix, i);
	return array;
}

IloNumVarArray CreateNumVarArray(IloEnv env, IloInt m, char* prefix, IloNum LB, IloNum UB)
{
	IloNumVarArray array = IloNumVarArray(env, m, LB, UB);
	for (int i = 0; i < m; i++)
		SetName(array[i], prefix, i);
	return array;
}

IloIntVarArray CreateIntVarArray(IloEnv env, IloInt m, char* prefix, IloNum LB, IloNum UB)
{
	IloIntVarArray array = IloIntVarArray(env, m, LB, UB);
	for (int i = 0; i < m; i++)
		SetName(array[i], prefix, i);
	return array;
}

NumArray2 CreateNumArray2(IloEnv env, IloInt m, IloInt n)
{
	NumArray2 array =  NumArray2(env, m);
	for (int i = 0; i < m; i++)
		array[i] = IloNumArray(env, n);
	return array;
}

IntArray2 CreateIntArray2(IloEnv env, IloInt m, IloInt n)
{
	IntArray2 array =  IntArray2(env, m);
	for (int i = 0; i < m; i++)
		array[i] = IloIntArray(env, n);
	return array;
}

NumVarArray2 CreateNumVarArray2(IloEnv env, IloInt m, IloInt n, char* prefix, IloNum LB, IloNum UB)
{
	NumVarArray2 array =  NumVarArray2(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = IloNumVarArray(env, n, LB, UB);
		for (int j = 0; j < n; j++)
			SetName2(array[i][j], prefix, i, j);
	}
	return array;
}

IntVarArray2 CreateIntVarArray2(IloEnv env, IloInt m, IloInt n, char* prefix, IloInt LB, IloInt UB)
{
	IntVarArray2 array =  IntVarArray2(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = IloIntVarArray(env, n, LB, UB);
		for (int j = 0; j < n; j++)
			SetName2(array[i][j], prefix, i, j);
	}
	return array;
}
BoolVarArray2 CreateBoolVarArray2(IloEnv env, IloInt m, IloInt n, char* prefix)
{
	BoolVarArray2 array =  BoolVarArray2(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = IloBoolVarArray(env, n);
		for (int j = 0; j < n; j++)
			SetName2(array[i][j], prefix, i, j);
	}
	return array;
}

NumVarArray3 CreateNumVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, char* prefix, IloNum LB, IloNum UB)
{
	NumVarArray3 array = NumVarArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateNumVarArray2(env, n, k, "", LB, UB);
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				SetName3(array[i][j][c], prefix, i, j, c);
	}
	return array;
}

NumArray3 CreateNumArray3(IloEnv env, IloInt m, IloInt n, IloInt k)
{
	NumArray3 array = NumArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateNumArray2(env, n, k);
	}
	return array;
}

IntArray3 CreateIntArray3(IloEnv env, IloInt m, IloInt n, IloInt k)
{
	IntArray3 array = IntArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateIntArray2(env, n, k);
	}
	return array;
}

BoolVarArray3 CreateBoolVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, char* prefix)
{
	BoolVarArray3 array = BoolVarArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateBoolVarArray2(env, n, k, "");
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				SetName3(array[i][j][c], prefix, i, j, c);
	}
	return array;
}

IntVarArray3 CreateIntVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, char* prefix, IloNum LB, IloNum UB)
{
	IntVarArray3 array = IntVarArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateIntVarArray2(env, n, k, "", LB, UB);
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				SetName3(array[i][j][c], prefix, i, j, c);
	}
	return array;
}

IntVarArray4 CreateIntVarArray4(IloEnv env, IloInt m, IloInt n, IloInt k, IloInt l, char* prefix, IloNum LB, IloNum UB)
{
	IntVarArray4 array = IntVarArray4(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateIntVarArray3(env, n, k, l, "", LB, UB);
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				for (int d = 0; d < l; d++)
					SetName4(array[i][j][c][d], prefix, i, j, c, d);
	}
	return array;
}

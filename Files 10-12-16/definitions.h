#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <ilcplex/ilocplex.h>

#include <vector>
#include <algorithm>
using namespace std;

typedef std::pair<int, int> Coordinate;
typedef vector<Coordinate> CoordinateArray;


typedef IloArray<IloNumArray> NumArray2;
typedef IloArray<IloBoolArray> BoolArray2;
typedef IloArray<IloIntArray> IntArray2;
typedef IloArray<IntArray2> IntArray3;
typedef IloArray<NumArray2> NumArray3;

typedef IloArray<IloNumVarArray> NumVarArray2;
typedef IloArray<IloBoolVarArray> BoolVarArray2;
typedef IloArray<IloIntVarArray> IntVarArray2;

typedef IloArray<NumVarArray2> NumVarArray3;
typedef IloArray<BoolVarArray2> BoolVarArray3;
typedef IloArray<IntVarArray2> IntVarArray3;

typedef IloArray<IntVarArray3> IntVarArray4;

NumArray2 CreateNumArray2(IloEnv env, IloInt m, IloInt n);
IntArray2 CreateIntArray2(IloEnv env, IloInt m, IloInt n);

IloNumVarArray CreateNumVarArray(IloEnv env, IloInt m, char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
IloIntVarArray CreateIntVarArray(IloEnv env, IloInt m, char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
IloBoolVarArray CreateBoolVarArray(IloEnv env, IloInt m, char* prefix);

NumVarArray2 CreateNumVarArray2(IloEnv env, IloInt m, IloInt n, char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
IntVarArray2 CreateIntVarArray2(IloEnv env, IloInt m, IloInt n, char* prefix, IloInt LB = 0, IloInt UB = IloInfinity);
BoolVarArray2 CreateBoolVarArray2(IloEnv env, IloInt m, IloInt n, char* prefix);

NumVarArray3 CreateNumVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
IntVarArray3 CreateIntVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
BoolVarArray3 CreateBoolVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, char* prefix);

IntArray3 CreateIntArray3(IloEnv env, IloInt m, IloInt n, IloInt k);
NumArray3 CreateNumArray3(IloEnv env, IloInt m, IloInt n, IloInt k);

IntVarArray4 CreateIntVarArray4(IloEnv env, IloInt m, IloInt n, IloInt k, IloInt l, char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);

void SetName(IloExtractable obj, char* prefix, int i);
void SetName2(IloExtractable obj, char* prefix, int i, int j);
void SetName3(IloExtractable obj, char* prefix, int i, int j, int k);
void SetName4(IloExtractable obj, char* prefix, int i, int j, int k, int l);

struct BoundingRectangle
{
	int top;
	int bottom;
	int left;
	int right;
};

struct BoundingRectangleLess
{
	bool operator()(const BoundingRectangle& br1, const BoundingRectangle& br2)
	{
		if(br1.top < br2.top)
			return true;
		else if (br1.top == br2.top 
				&& br1.bottom < br2.bottom)
			return true;
		else if (br1.top == br2.top && br1.bottom == br2.bottom
				&& br1.left < br2.left)
			return true;
		else if (br1.top == br2.top && br1.bottom == br2.bottom && br1.left == br2.left
			&& br1.right < br2.right)
			return true;
		else
			return false;
	}
};

struct BoundingRectangleEqual
{
	bool operator()(const BoundingRectangle& br1, const BoundingRectangle& br2)
	{
		if (br1.top == br2.top && br1.bottom == br2.bottom && br1.left == br2.left && br1.right == br2.right)
			return true;
		else
			return false;
	}
};

#endif


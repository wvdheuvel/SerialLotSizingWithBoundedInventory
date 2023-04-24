#ifndef DPS_H
#define DPS_H

// solve MLS by DP
double solveMLSbyDP(Data *inst);

// solve uncpacitated problem: programmed by Eric (may contain errors!)
double solveZangwill(int m, int n, double *d, double **consCost, double **varCost, double **invCost);

// solve uncpacitated problem: programmed by Wilco
double solveZangwill2(int m, int n, double *d, double **K, double **p, double **h);

#endif
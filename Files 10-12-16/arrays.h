#ifndef ARRAYS_H
#define ARRAYS_H

// creates arrays of different dimensions
void SetDim2Dbl(double **(&x), int n1, int n2);
void SetDim2Int(int **(&x), int n1, int n2);
void SetDim3Dbl(double ***(&x), int n1, int n2, int n3);
void SetDim4Dbl(double ****(&x), int n1, int n2, int n3, int n4);
void DeleteDim6Dbl(double ******(&x), int n1, int n2, int n3, int n4, int n5, int n6);
void DeleteDim6Int(int ******(&x), int n1, int n2, int n3, int n4, int n5, int n6);
void DeleteDim2Dbl(double **(&x), int n1, int n2);
void DeleteDim2Int(int **(&x), int n1, int n2);
void DeleteDim3Dbl(double ***(&x), int n1, int n2, int n3);
void DeleteDim4Dbl(double ****(&x), int n1, int n2, int n3, int n4);
void SetDim6Dbl(double ******(&x), int n1, int n2, int n3, int n4, int n5, int n6);
void SetDim6Int(int ******(&x), int n1, int n2, int n3, int n4, int n5, int n6);

#endif
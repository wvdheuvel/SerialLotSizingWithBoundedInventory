#ifndef ARRAYS_H
#define ARRAYS_H

void SetDim2Dbl(double **(&x), int n1, int n2);
void SetDim3Dbl(double ***(&x), int n1, int n2, int n3);
void SetDim4Dbl(double ****(&x), int n1, int n2, int n3, int n4);
void DeleteDim2Dbl(double **(&x), int n1, int n2);
void DeleteDim3Dbl(double ***(&x), int n1, int n2, int n3);
void DeleteDim4Dbl(double ****(&x), int n1, int n2, int n3, int n4);

#endif
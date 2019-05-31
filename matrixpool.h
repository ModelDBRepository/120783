#ifndef _MATRIXPOOL_H
#define _MATRIXPOOL_H
void NormTopMat(double **Mat, double A, double wd, const int row, const int col);
void RandTopMat(double **Mat, double A, double wd, const int row, const int col);
void nRandTopMat(double **Mat, double A, double wd, const int row, const int col);
void AntiTopMat(double **Mat, double wd,double depth, const int row, const int col);
void CosTopMat(double **Mat, double wd, const int row, const int col);
void NulMat(double **Mat, const int row, const int col);
void iMat(double **Mat, const int row, int col);
void ShowMat(double **Mat, const int row, const int col);
double dotprod(int n, double *vec1, double *vec2);

#endif
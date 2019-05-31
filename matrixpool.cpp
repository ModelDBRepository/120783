#include "matrixpool.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "rng.h"
#include <ctime>
using namespace std;

void NormTopMat(double **Mat, double A, double wd, const int row, const int col)
{
	int i=0,j=0;
	for (i=1;i<row;i++)
		Mat[i]=Mat[i-1]+col;
	for (i=0;i<row;i++){
		for(j=0;j<col;j++){
			 Mat[i][j]=A*exp(-pow(double(i-j),2)/(2*pow(wd,2)));
		}
	}
}

//void nRandTopMat(double **Mat,  double A, double wd, const int row, const int col)
//{
//	int i=0,j=0;
//	extern int seed;
//	for (i=1;i<row;i++)
//		Mat[i]=Mat[i-1]+col;
//	for (i=0;i<row;i++){
//		for(j=0;j<col;j++){
//			if(abs(i-j)<wd)
//			 Mat[i][j]=A*fabs(gsdev(seed));
//			else
//			Mat[i][j]=0;
//		}
//	}
//}
//void RandTopMat(double **Mat,  double A, double wd, const int row, const int col)
//{
//	int i=0,j=0;
//	 //int seed=int(time(0));
//	extern int seed;
//	for (i=1;i<row;i++)
//		Mat[i]=Mat[i-1]+col;
//	for (i=0;i<row;i++){
//		for(j=0;j<col;j++){
//			if(abs(i-j)<wd)
//			 Mat[i][j]=A*ran0(seed);
//			else
//			Mat[i][j]=0;
//		}
//	}
//}

void AntiTopMat(double **Mat, double wd,double depth, const int row, const int col)
{
	int i=0,j=0;
	for (i=1;i<row;i++)
		Mat[i]=Mat[i-1]+col;
	for (i=0;i<row;i++){
		for(j=0;j<col;j++){
			Mat[i][j]=(1-depth*exp(-pow(double(i-j),2)/(2*pow(wd,2))));
		}
	}
}

void CosTopMat(double **Mat, double wd, const int row, const int col)
{
	int i=0,j=0;
	for (i=1;i<row;i++)
		Mat[i]=Mat[i-1]+col;
	for (i=0;i<row;i++){
		for(j=0;j<col;j++){
			Mat[i][j]=(1+cos(2*3.1415/wd*(i-j+wd/2)))/2;
		}
	}

}

void NulMat(double **Mat, const int row, const int col)
{
	int i=0,j=0;
	for (i=1;i<row;i++)
		Mat[i]=Mat[i-1]+col;
	for (i=0;i<row;i++){
		for(j=0;j<col;j++){
			Mat[i][j]=0;
		}
	}
}

void iMat(double **Mat, const int row, int col)
{
	int i=0,j=0;
	for (i=1;i<row;i++)
		Mat[i]=Mat[i-1]+col;
	for (i=0;i<row;i++){
		for(j=0;j<col;j++){
			Mat[i][j]=(i==j);
		}
	}
}


void ShowMat(double **Mat, const int row, const int col)
{
	int i,j;
	ofstream wMatrix("wMatrix.txt", ios::app);
	for (i=0;i<row;i++){
		for(j=0;j<col;j++)
		{
			wMatrix<<Mat[i][j]<<' ';
		}
		wMatrix<<'\n';
	}
	wMatrix.close();
}


double dotprod(int n, double *vec1, double *vec2){
       int k, m;
       double sum=0.0;
       k=n/6;
       m=n%6;
       while(k--){
       sum+=*vec1*(*vec2);
       sum+=*(vec1+1)**(vec2+1);
       sum+=*(vec1+2)**(vec2+2);
       sum+=*(vec1+3)**(vec2+3);
       sum+=*(vec1+4)**(vec2+4);
       sum+=*(vec1+5)**(vec2+5);
       vec1+=6;
       vec2+=6;}
       while(m--)
       sum+=*vec1++**vec2++;
       return sum;
       }


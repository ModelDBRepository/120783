#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <complex>
using namespace std;

template<class T>
class vect
{
private:
	int dim;
	T *v;
public:
	vect();
	vect(int size, T val=0);
	vect(int size, const T *a);
	~vect();
	vect & operator=(const vect & vrhs);
	T & operator()(int id);
	const T & operator()(int id) const;
	T operator*(const vect &vout);
	int ndim(){return dim;}
	void disp();
	void disp(ofstream & outfile);
	//vect operator * (const matrix &vmat);

};


template<class T>
class matrix{
	T **pmat;
	int row, col;
public:
	matrix();
	matrix(int rowsize, int colsize, const T *a);
	matrix(int rowsize, int colsize, T val=0);
	~matrix();	
	T& operator()(int ri, int ci);
	matrix& operator=(const matrix & mrhs);
	matrix& operator*(const matrix & mrhs);	
	void disp();
	void disp(ofstream & outfile);
	int nrows(){return row;}
	int ncols(){return col;}
};

template<class T>
vect<T>::vect()
:dim(0),v(0)
{}



template<class T>
vect<T>::vect(int size, T val){
	v=new T[dim=size];
	for(int i=0;i<dim;i++)
	{
		v[i]=val;
	}
}


template<class T>
vect<T>::vect(int size, const T *a){
	v=new T[dim=size];
	for(int i=0;i<dim;i++)
	{
		v[i]=*a++;
	}

}

template<class T>
vect<T>::~vect(){
	delete [] v;
}

template<class T>
void vect<T>::disp(ofstream & outfile)
{
	int i=0;
	for(i=0;i<dim; i++)
	outfile<<v[i]<<" "<<"\n";
	//outfile<<"\n";
}


template<class T>
void vect<T>::disp()
{
	int i=0;
	for(i=0;i<dim; i++)
	cout<<v[i]<<" ";
	cout<<"\n";
}


template<class T>
T & vect<T>::operator()(int id)
{
	return v[id];
}

template<class T>
vect<T> & vect<T>::operator =(const vect<T> &vrhs)
{
	if(this!=&vrhs){
		int i;
		if(dim!=vrhs.dim)
		{
			if(v!=0){
				delete[] v;				
			}
		}
		dim=vrhs.dim;
		v=new T [dim];

		for(i=0;i<dim;i++)
			v[i]=vrhs.v[i];
	}
	return *this;
}

//template<class T>
//T vect<T>::operator *(const vect<T> &vout)
//{
//	T sum=0;
//	if(dim!=vout.dim)
//		cerr<< "vectors of different dimensions";
//	else
//	{
//		for (int i=0;i<vout.dim;i++)
//			sum+=v[i]*vout.v[i];
//	}
//	return sum;
//}

template<class T>
matrix<T>::matrix()
:row(0), col(0), pmat(0){}



template<class T>
matrix<T>::matrix( int rowsize, int colsize, const T *a){

	int i,j;
	row=rowsize;
	col=colsize;
	pmat=new T *[rowsize];
	pmat[0]=new T [rowsize*colsize];
	for(i=1;i<rowsize;i++)
		pmat[i]=pmat[i-1]+colsize;
	for (i=0;i<rowsize;i++)
		for (j=0;j<colsize;j++)
			pmat[i][j]=*a++;


}

template<class T> //constant matrix
matrix<T>::matrix(int rowsize, int colsize, T val){
	int i,j;
	row=rowsize;
	col=colsize;
	pmat=new T *[colsize];
	pmat[0]=new T [rowsize*colsize];
	for(i=1;i<rowsize;i++)
		pmat[i]=pmat[i-1]+colsize;
	for (i=0;i<rowsize;i++)
		for (j=0;j<colsize;j++)
			pmat[i][j]=val;
	
}


template<class T> //destructor
matrix<T>::~matrix(){
	delete [] pmat[0];
	delete [] pmat;
}

template<class T> //refer to an element in row ri column ci
T & matrix<T>::operator ()(int ri, int ci)
{
	return pmat[ri][ci];
}


//template<class T>
//vect matrix<T>::rowslice(int rn)
//{
//	vect temp(col);
//		for(int i=0;i<col;i++)
//			temp.v[i]=pmat[rn][i];
//	return temp;
//
//}

//template<class T>
//vect matrix<T>::colslice(int cn)
//{
//	vect temp(row);
//		for(int i=0;i<row;i++)
//			temp.v[i]=pmat[i][cn];
//	return temp;
//
//}

template<class T>
matrix<T> & matrix<T>::operator =(const matrix<T> &mrhs)
{
	if(this!=&mrhs){
		int i,j;
		if(row!=mrhs.row|| col!=mrhs.col)
		{
			if(pmat!=0){
				delete[](pmat[0]);
				delete[](pmat);
			}
		}
		row=mrhs.row;
		col=mrhs.col;
		pmat=new T *[row];
		pmat[0]=new T [row*col];
	
	for(i=1;i<row;i++)
		pmat[i]=pmat[i-1]+col;
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
			pmat[i][j]=mrhs.pmat[i][j];
	}
	return *this;
}

//template<class T>
//matrix<T>& matrix<T>::operator *(const matrix<T> & mrhs)
//{
//	int i, j, k;
//	matrix temp(row,mrhs.col);
//	if (mrhs.row!=col)
//      cerr<<"size is not matched";
//	else{
//		for (i=0; i <row; i++){
//			for (j=0; j<mrhs.col; j++)
//			{
//				temp.pmat[i][j]=0;
//				for (k=0; k<col; k++)
//					temp.pmat[i][j]+=pmat[i][k] * mrhs.pmat[k][j];
//			}
//		}
//	}
//	*this=temp;
//	return *this;
//}
//

template<class T>
void matrix<T>::disp(ofstream & outfile)
{
	int i,j;
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			outfile<<pmat[i][j]<<" "<<"\n";}
		//outfile<<"\n";
	}
}

template<class T>
void matrix<T>::disp()
{
	int i,j;
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			cout<<pmat[i][j]<<" ";}
		cout<<"\n";
	}
}


#endif
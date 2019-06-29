#include "../include/DNewtonSolver.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using std::cout;
using std::endl;
using std::vector;



//Calculate the Jacobian matrix
void Jacobian(void(*pFunc)(double *, void *, double *), double * input, void * para, const int n, double ** jac){
	double precision = 2.2E-16;
	double delta;
	double *x = new double[n];
	double *fi = new double[n];
	double *ff = new double[n];
	(*pFunc)(input, para, fi);
	for(int i=0; i<n; ++i){
		x[i] = input[i];
	}

	for(int j=0; j<n; ++j){
		delta = sqrt(precision*fabs(x[j]));
		if (delta==0) delta=precision;
		x[j] += delta;
		(*pFunc)(x, para, ff);
		for(int i=0; i<n; ++i){
			jac[i][j] = (ff[i]-fi[i])/delta;
		}
		x[j] = input[j];
	}

	delete[] x;
	delete[] fi;
	delete[] ff;

}
//Check wehther the residual is small enough, less than epsabs
int checkResidual(double *f, const int n, double epsabs){
	double residual = 0.0;
	if(epsabs<0.0){
		cout<<"Absolute tolerence is negative!"<<endl;
		exit(EXIT_FAILURE);
	}

	for(int i=0; i<n; ++i){
		residual += fabs(f[i]);
	}

	if(residual<epsabs) return 1;
	return 0;
}

// LU decomposition, this algorithm can be found in Numerical Recipes in C, 2nd ed.
// a[n][n] is the matrix for decomposition.
// idx[n] records the row permutation effected by the partial pivoting
// d = 1 or -1, number of row interchanges is even or odd
void LUDcmp(double ** a, const int n, int *idx, int &d){
	int imax;
	double big, dum, sum, temp;
	double *vv = new double[n];		//the implicit scaling of each row.
	double tiny = 1.0e-20;

	d = 1; 		//No row interchanges yet.
	for (int i=0; i<n; ++i){
		big = 0.0;
		for (int j=0; j<n; ++j){
			if((temp=fabs(a[i][j]))>big) big = temp;
		}
		if (big==0) {			//The largest element is zero
			cout<<"Singular matrix in routine LUDcmp"<<endl;
			delete[] vv;
			exit(EXIT_FAILURE);
		}

		vv[i] = 1/big;		//Scaling
	}

	for (int j=0; j<n; ++j){		//Loop over columns of Crout's method
		for (int i=0; i<j; ++i){
			sum = a[i][j];
			for (int k=0; k<i; ++k) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}

		//Search for largest pivot element
		big = 0.0;
		for (int i=j; i<n; ++i){
			sum = a[i][j];
			for (int k=0; k<j; ++k) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ( (dum=vv[i]*fabs(sum))>=big){
				big = dum;
				imax = i;
			}
		}

		if (j != imax){			// Do we need to interchange rows?
			for (int k=0; k<n; ++k){	//Yes
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;				//change the parity of d;
			vv[imax] = vv[j];	//Interchange the scale factor;
		}

		idx[j]=imax;

		if (a[j][j]==0)	a[j][j] = tiny;	//If the pivot element is zero, submitted by a tiny value

		if	(j!=(n-1)){						//Divide by the pivot element
			dum = 1.0/(a[j][j]);
			for (int i=j+1; i<n; ++i) a[i][j] *= dum;
		}
	}										//Go back to the next column in the reduction

	delete[] vv;
}

//Solves the set of n linear equations AX = B. Here the input a[n][n] is not the matrix A, but
//its LU decomposition. And idx[n] is the permutation vector. b[n] is the right-hand size vector
//B, and returns the solution vector X.
//This algorithm can be found in Numerical Recipes in C, 2nd ed.
void lubksb(double **a, const int n, int *idx, double *b){

	int ii=-1, ip;
	double sum;
	for(int i=0; i<n; ++i){
		ip = idx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(ii+1)
			for(int j=ii; j<=i-1; ++j)	sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i] = sum;
	}

	for (int i=n-1; i>=0; --i){
		sum = b[i];
		for(int j=i+1;j<n; ++j)	{
			sum -= a[i][j]*b[j];
		}

		b[i] = sum/a[i][i];
	}
}

//Use Newton Iterator to find the roots of the linear equation f(x) = 0.
int NewtonIter(void(*pFunc)(double *, void *, double *), double * input, void * para, const int n, const int nIter, const double delta){
	double * f = new double[n];
	int d;
	int * idx = new int[n];
	double **jac = new double *[n];
	for(int i=0; i<n; ++i) jac[i] = new double[n];

	(*pFunc)(input, para, f);
	for(int i=0; i<nIter;++i){
		if (checkResidual(f,n,delta)){
			for(int i=0; i<n;++i) delete[] jac[i];
			delete[] jac;
			delete[] idx;
			delete[] f;
			return 0;
		}
		else{
			Jacobian(pFunc, input, para, n, jac);
			LUDcmp(jac, n, idx, d);
			lubksb(jac, n, idx, f);

			for(int i=0; i<n; ++i) input[i] -= f[i];
			(*pFunc)(input, para, f);
		}
	}

	if (checkResidual(f,n,delta)){
		for(int i=0; i<n;++i) delete[] jac[i];
		delete[] jac;
		delete[] idx;
		delete[] f;
		return 0;
	}

	for(int i=0; i<n;++i) delete[] jac[i];
	delete[] jac;
	delete[] idx;
	delete[] f;
	cout<<"Reached the maximum iteration number: " <<nIter<<endl;
	return 1;

}






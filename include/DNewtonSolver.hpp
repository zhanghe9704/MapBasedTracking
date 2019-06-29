#ifndef DNEWTONSOLVER_HPP
#define DNEWTONSOLVER_HPP

#include <vector>
#include <iostream>
#include "map.hpp"


using std::vector;

//Calculate the Jacobian matrix
void Jacobian(void(*pFunc)(double *, void *, double *), double * input, void * para, const int n, double ** jac);
//Check wehther the residual is small enough, less than epsabs
int checkResidual(double *f, const int n, double epsabs);
// LU decomposition, this algorithm can be found in Numerical Recipes in C, 2nd ed. 
// a[n][n] is the matrix for decomposition. 
// idx[n] records the row permutation effected by the partial pivoting
// d = 1 or -1, number of row interchanges is even or odd 
void LUDcmp(double ** a, const int n, int *idx, int &d);
//Solves the set of n linear equations AX = B. Here the input a[n][n] is not the matrix A, but
//its LU decomposition. And idx[n] is the permutation vector. b[n] is the right-hand size vector 
//B, and returns the solution vector X. 
//This algorithm can be found in Numerical Recipes in C, 2nd ed. 
void lubksb(double **a, const int n, int *idx, double *b);

//Use Newton Iterator to find the roots of the linear equation f(x) = 0.
int NewtonIter(void(*pFunc)(double *, void *, double *), double * input, void * para, const int n, const int nIter, const double delta);

#endif

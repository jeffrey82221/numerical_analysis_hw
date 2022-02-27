//
//  MAT.h
//  Numerical_Analysis
//
//  Created by jeffrey_MAC on 2015/3/6.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#ifndef __Numerical_Analysis__MAT__
#define __Numerical_Analysis__MAT__
#include <iostream>
using namespace std;
#include <stdio.h>
#include "VEC.h"
class MAT {
private:
    int n;      // define nxn matrix
    VEC **va;   // array of n pointers to vectors

public:
    MAT(int dim);   // uninit constructor
    MAT(const MAT &m1); // copy constructor
    MAT(int dim,double *v); // init constructor
    ~MAT(); // destructor
    int dim();  // return dimension of the matrix
    MAT tpose();    // transpose
    MAT &operator-();   // unary operator, negative value
    MAT &operator=(MAT m1); // assignment
    MAT &operator+=(MAT &m1);   // m += m1;
    MAT &operator-=(MAT &m1);   // m -= m1;
    MAT &operator*=(double a);  // m *= dbl;
    MAT &operator/=(double a);  // m /= dbl
    MAT operator+(MAT m1);      //m1+m2
    MAT operator-(MAT m1);      //m1-m2
    MAT operator*(MAT m1);      //m1*m2
    VEC & operator[](int m);    //m'th row
    VEC operator*(VEC v1);      //m*v1
    MAT operator*(double a);    //m*dbl
    MAT operator/(double a);    //m/dbl
    friend MAT operator*(double a,MAT &m1); // dbl x m
    friend VEC operator*(VEC &v1,MAT &m1);  // vT x m
    friend ostream& operator<<(ostream&, const MAT&);
    
    int p; //norm of p
};
MAT operator*(double a,const MAT &m1);      // dbl x m
VEC operator*(VEC &v1,MAT &m1);             // vT x m
MAT eye(int dim);
MAT &luFact(MAT &m1); // LU decomposition
MAT &luBanded(MAT &m1,int m); //banded LU decomposition
VEC fwdSubs(MAT &m1,VEC b); // forward substitution
VEC bckSubs(MAT &m1,VEC b); // backward substitution
MAT &cholesky(MAT &A); // Cholesky decomposition
VEC choSolve(MAT &L,VEC b); // forward and backward substitutions
int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol); //Jacobi Iteration Method
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol); //Gauss-Seidel Iteration Method
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol); //Symmetric Gauss-Seidel Iteration Method
int cg(MAT &A,VEC b,VEC &x,int maxIter, double tol); //Conjugate Gradient Method
//Function below implement the power methods of finding different eigen value and eigen vector
int powerMethod(MAT &A,double &v,VEC &q,double maxIter,double tol); //Power Method
int invPowerMethod(MAT &A,double &v,VEC &q,double maxIter,double tol);//Inverse Power Method
int sinvPowerMethod(MAT &A,double &v,VEC &q,double maxIter,double tol,double omega);//Inverse Powever Method with shifting
void QRDecomposition(MAT &A,MAT &Q,MAT &R); //QR decomposition
int EVqr(MAT &A,double tol,int maxiter); //QR iteration for finding all eigen value and eigen vector
int EVqrShifted(MAT &A,double mu,double tol,int maxiter); //QR iteration for finding all eigen value and eigen vector with shifting
double EVqrError(MAT &A); //calculate the error for qr iteration

void splineM(int N,VEC &X,VEC &Y,VEC &M);   // generate spline momentum M
double spline(double x,int N,VEC &X,VEC &Y,VEC &M); // spline interp at x
double integ(VEC &X,VEC &Y,int n);  // composite nth order Newton-Cotes integral
VEC newtoncotesWeights(int n); //compute the weight vector of Newton-cotes
VEC coefficientofNCWeights(VEC &in);
int chosensize(int m,int n);
void recursivechoice(VEC &result,VEC &in,int chosen_size,int begin,double product_value,int &result_index);
#endif /* defined(__Numerical_Analysis__MAT__) */

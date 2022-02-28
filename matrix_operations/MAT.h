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
};
MAT operator*(double a,const MAT &m1);      // dbl x m
VEC operator*(VEC &v1,MAT &m1);             // vT x m


#endif /* defined(__Numerical_Analysis__MAT__) */

//
//  VEC.h
//  Numerical_Analysis
//
//  Created by jeffrey_MAC on 2015/3/6.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#ifndef Numerical_Analysis_VEC_h
#define Numerical_Analysis_VEC_h
#include <iostream>
using namespace std;
class VEC {
private:
    int dim;    // vector length
    double *val;// array to store vector
public:
    VEC(int n);                     //uninit constructor , vec set to 0
    VEC(const VEC &v1);             //copy constructor
    VEC(int n,double *v);           //init constructor
    ~VEC();                         //destructor
    int len();                      //dimension of the vector
    VEC operator-();                //unary operator , negative value
    VEC &operator=(const VEC v1);   //assignment
    VEC &operator+=(const VEC v1);  // V += v1;
    VEC &operator-=(const VEC v1);  // V -= v1;
    VEC &operator*=(double a);      // V *= dbl;
    VEC &operator/=(double a);      // V /= dbl;
    VEC operator+(const VEC v1);    // V + v1
    VEC operator-(const VEC v1);    // V - v1
    double operator*(VEC v1);       //inner product
    VEC operator*(double a);        //scale product
    VEC operator/(double a);        //scale division
    VEC operator()(int begin,int end);
    double &operator[](int n);      //indexing
    double norm(int p);             //p-norm of the vector
    friend VEC operator*(double a,const VEC v1);    // dbl x V
    friend VEC *newVEC(int n);      // alloc memory for VEC
    friend ostream& operator<<(ostream&, const VEC&);
    void sort();
};
VEC operator*(double a,const VEC v1);
VEC *newVEC(int n); // alloc memory for VEC .
VEC abs(VEC v1);
double max(VEC v1);
double Lagrange(double x,VEC &XDATA,VEC &YDATA);
double DDif(VEC &XS, VEC &YS, VEC &A,int i0,int ik);
double PiecewiseLinear(double x,VEC &XS,VEC &YS);

#endif
//
//  VEC.cpp
//  Numerical_Analysis
//
//  Created by jeffrey_MAC on 2015/3/6.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "VEC.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
VEC::VEC(int n){
    dim = n;
    val=(double *) calloc(n,sizeof(double));
}//uninit constructor , vec set to 0
VEC::VEC(const VEC &v1){
    dim=v1.dim;
    val=(double *) calloc(dim,sizeof(double));
    for (int i=0; i<dim; i++) {
        val[i]=v1.val[i];
    }
}//copy constructor
VEC::VEC(int n,double *v){
    dim=n;
    val=(double *) calloc(n,sizeof(double));
    for (int i=0; i<n; i++) val[i]=v[i];
}//init constructor
VEC::~VEC(){
    free(val);
}//destructor
int VEC::len(){
    return dim;
}//dimension of the vector
VEC VEC::operator-(){
    for (int i=0; i<dim; i++) val[i]=-val[i];
    return *this;
}//unary operator , negative value
VEC &VEC::operator=(const VEC v1){
    dim=v1.dim;
    for (int i=0; i<dim; i++) {
        val[i]=v1.val[i];
    }
    return *this;
}//assignment
VEC &VEC::operator+=(const VEC v1){
    for (int i=0; i<dim; i++) {
        val[i]+=v1.val[i];
    }
    return *this;
}// V += v1;
VEC &VEC::operator-=(const VEC v1){
    for (int i=0; i<dim; i++) {
        val[i]-=v1.val[i];
    }
    return *this;
}// V -= v1;
VEC &VEC::operator*=(double a){
    for (int i=0; i<dim; i++) {
        val[i]*=a;
    }
    return *this;
}// V *= dbl;
VEC &VEC::operator/=(double a){
    for (int i=0; i<dim; i++) {
        val[i]/=a;
    }
    return *this;
}// V /= dbl;
VEC VEC::operator+(const VEC v1){
    VEC s(*this);
    for (int i=0; i<dim; i++) s.val[i]+=v1.val[i];
    return s;
}// V + v1
VEC VEC::operator-(const VEC v1){
    VEC s(*this);
    for (int i=0; i<dim; i++) s.val[i]-=v1.val[i];
    return s;
}// V - v1
double VEC::operator*(VEC v1){
    double ans=0;
    for (int i =0; i<dim; i++) {
        ans+=val[i]*v1.val[i];
    }
    return ans;
}//inner product
VEC VEC::operator*(double a){
    VEC s(*this);
    for (int i =0; i<dim; i++) {
        s.val[i] = val[i]*a;
    }
    return s;
}//scale product
VEC VEC::operator/(double a){
    VEC s(*this);
    for (int i =0; i<dim; i++) {
        s.val[i] = val[i]/a;
    }
    return s;
}//scale division
double &VEC::operator[](int n){
    if (n<0) n=0;
    else if (n>=dim) n=dim-1;
    return val[n];
}//indexing
VEC VEC::operator()(int beg,int end){
    VEC r(end-beg+1);
    for (int i = 0; i<end-beg+1; i++) {
        r[i] = val[i+beg];
    }
    return r;
}
double VEC::norm(int p){ //for p = 0 : p = Infinite
    double result = 0;
    if (p==0) { //for p = 0 : we execute the norm of Infinite (find the max of the vector)
        for (int i = 0; i<dim; i++){
            double max = fabs(val[i]);
            if (max>result) {
                result = max;
            }
    
        }
    }
    else{ //else we perform the normal p-norm
        for (int i = 0; i<dim; i++) {
            result+=pow(fabs(val[i]),p);
        }
        result = pow(result,1.0/p);
    
    }
    return result;
    
}//p-norm of the vector
void VEC::sort(){
    double tmp;
    for (int i = 1; i<len(); i++) {
        for (int j = 0; j<i; j++) {
            if (val[i]>val[j]) {
                tmp = val[j];
                val[j] = val[i];
                val[i] = tmp;
            }
        }
    }
    
    
}
VEC operator*(double a,const VEC v1){
    VEC s(v1.dim);
    for (int i =0; i<v1.dim; i++) {
        s.val[i] = v1.val[i]*a;
    }
    return s;
}
VEC *newVEC(int n){
    VEC *vptr;
    vptr=(VEC *)malloc(sizeof(VEC));
    vptr->dim=n;
    vptr->val=(double*)calloc(n,sizeof(double));
    return vptr;
}// alloc memory for VEC .
ostream& operator<<(ostream& os, const VEC& v1){

    for (int i=0; i<v1.dim; i++) {
        os<<fixed << setw( 11 ) << setprecision( 10 ) <<v1.val[i]<<"   ";
    }
    return os;
}

double Lagrange(double x,VEC &XS,VEC &YS){
    int n = XS.len();
    VEC NS(YS);
    int j,k;
    for (k=1; k<n; k++) {
        for (j=0; j<n-k; j++) {
            NS[j] = ((x-XS[j])*NS[j+1]-(x-XS[k+j])*NS[j])/(XS[j+k]-XS[j]);
        }
    }
    return NS[0];
}
double PiecewiseLinear(double x,VEC &XS,VEC &YS){
    int n = XS.len();
    VEC xs2(2),ys2(2);
    
    for (int i = 0; i<n-1; i++) {
        if (x>XS[i]&&x<XS[i+1]) {
            xs2 = XS(i,i+1);
            ys2 = YS(i,i+1);
            return Lagrange(x,xs2,ys2);
            break;
        }else if(x==XS[i]){
            return YS[i];
        }
    }
    return YS[n-1];
}
VEC abs(VEC v1){
    VEC r(v1.len());
    for (int i = 0; i<v1.len(); i++) {
        r[i] = fabs(v1[i]);
    }
    return r;
}
double max(VEC v1){
    double r=-INFINITY;
    for (int i = 0; i<v1.len(); i++) {
        if (v1[i]>r) {
            r = v1[i];
        }
    }
    return r;
}
double DDif(VEC &XS, VEC &YS, VEC &A,int i0,int ik){//Divided Difference : A is the coefficient of interpolated polynomial
    double result;
    if (i0==ik) {
        result = YS[i0];
    }else{
        result = (DDif(XS,YS, A, i0+1, ik)-DDif(XS,YS, A, i0, ik-1))/(XS[ik]-XS[i0]);
    }
    if (i0==0) {
        A[ik]=result;
    }
    return result;
    
}

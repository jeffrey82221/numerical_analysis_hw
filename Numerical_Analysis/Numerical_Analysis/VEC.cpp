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
    os<<"[";
    for (int i=0; i<v1.dim-1; i++) {
        os<<v1.val[i]<<",";
    }
    os<<v1.val[v1.dim-1];
    os<<"]";
    return os;
}

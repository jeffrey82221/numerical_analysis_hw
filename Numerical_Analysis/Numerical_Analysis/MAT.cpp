//
//  MAT.cpp
//  Numerical_Analysis
//
//  Created by jeffrey_MAC on 2015/3/6.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "MAT.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
MAT::MAT(int dim)
{
    n=dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++) {
        va[i]=newVEC(n);
    }
}
MAT::MAT(const MAT &m1){
    VEC **vsrc=m1.va;
    n = m1.n;
    va = (VEC **)malloc(n*sizeof(VEC*));
    for (int i = 0; i<n; i++) {
        va[i] = newVEC(n);
        (*va[i])= (*vsrc[i]);
    }
}
MAT::MAT(int dim,double *v){
    n = dim;
    va=(VEC **)malloc(n*sizeof(VEC*));
    for (int i=0; i<n; i++) {
        va[i]=newVEC(n);
        for (int j=0; j<n; j++) {
            (*va[i])[j]=*(v++);
        }
    }
}
MAT::~MAT(){
    for (int i=n-1; i>=0; i--) free(va[i]);
    free(va);
}
int MAT::dim() {
    return n;
}

MAT MAT::tpose(){
    MAT mnew(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            mnew[i][j]=(*va[j])[i];
        }
    }
    return mnew;
}
MAT & MAT::operator-() {
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            (*va[i])[j]=-(*va[i])[j];
    return *this;
}
MAT &MAT::operator=(MAT m1)
{
    for (int i=0; i<n; i++)
        (*va[i])=m1[i];
    return *this;
}

MAT &MAT::operator+=(MAT &m1){
    for (int i=0; i<n; i++)
        (*va[i])+=m1[i];
    return *this;
}
MAT &MAT::operator-=(MAT &m1){
    for (int i=0; i<n; i++)
        (*va[i])-=m1[i];
    return *this;
}
MAT &MAT::operator*=(double a){
    for (int i=0; i<n; i++)
        (*va[i])*=a;
    return *this;
}
MAT &MAT::operator/=(double a){
    for (int i=0; i<n; i++)
        (*va[i])/=a;
    return *this;
}

MAT MAT::operator+(MAT m1)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])+m1[i];
    return s;
}
MAT MAT::operator-(MAT m1)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])-m1[i];
    return s;
}
MAT MAT::operator*(MAT m1)  {
    MAT z(n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            z[i][j]=0;
            for (int k=0; k<n; k++)
                z[i][j]+=((*va[i])[k]*m1[k][j]);
        }
    }
    return z;
}


VEC &MAT::operator[](int m){
    if (m<0) m=0;
    else if (m>=n) m=n-1;
    return *va[m];
}//indexing

VEC MAT::operator*(VEC v1)
{
    VEC s(n);
    for (int i=0; i<n; i++) {
        s[i]=(*va[i])*v1;
    }
    return s;
}
MAT MAT::operator*(double a)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])*a;
    return s;
}
MAT MAT::operator/(double a)
{
    MAT s(n);
    for (int i=0; i<n; i++)
        s[i]=(*va[i])/a;
    return s;
}
MAT operator*(double a,MAT &m1){
    MAT s(m1.n);
    for (int i=0; i<m1.n; i++) {
        s[i] = m1[i]*a;
    }
    return s;
}


VEC operator*(VEC &v1,MAT &m1){
    VEC v2(m1.n);
    for (int i=0; i<m1.n; i++) {
        v2[i]=0;
        for (int j=0; j<m1.n; j++) {
            v2[i] += v1[j]*m1[j][i];
        }
    }
    return v2;
}

ostream& operator<<(ostream& os, const MAT& m1){
    for (int i =0; i<m1.n-1; i++) {
        VEC *v = m1.va[i];
        os<<*v<<","<<endl;
    }
    os<<*(m1.va[m1.n-1]);
    return os;
}

MAT eye(int dim){
    MAT eye(dim);
    for (int i = 0; i<dim; i++) {
        eye[i][i]=1;
    }
    return eye;
}

MAT &luFact(MAT &m1){

    for (int i =0; i<m1.dim(); i++) {
    //copy a[i][j] to u[i][j] need no action due to in-place LU
        for (int j = i+1; j<m1.dim(); j++) {
            m1[j][i] /= m1[i][i];
            
        }
        for (int j = i+1; j<m1.dim(); j++) {
            for (int k = i+1; k<m1.dim(); k++) {
                m1[j][k] -= m1[j][i]*m1[i][k];
            }
        }
    }
    return m1;
}

MAT &luBanded(MAT &m1,int m){
    for (int i =0; i<m1.dim(); i++) {
        //copy a[i][j] to u[i][j] need no action due to in-place LU
        for (int j = i+1; j<i+m && j<m1.dim(); j++) {
            m1[j][i] /= m1[i][i];
            
        }
        for (int j = i+1;j<i+m && j<m1.dim(); j++) {
            for (int k = i+1;k<i+m && k<m1.dim(); k++) {
                m1[j][k] -= m1[j][i]*m1[i][k];
            }
        }
    }
    return m1;
}

VEC fwdSubs(MAT &m1,VEC b){
    VEC y(b);
    for (int i = 0; i<m1.dim(); i++) {
        for (int j = i+1; j<m1.dim(); j++) {
            y[j] -= m1[j][i]*y[i];
        }
    }
    return y;
}// forward substitution

VEC bckSubs(MAT &m1,VEC y){
    VEC x(y);
    for (int i = m1.dim()-1; i>=0; i--) {
        x[i] /= m1[i][i];
        for (int j = i-1; j>=0; j--) {
            x[j] -= m1[j][i]*x[i];
        }
        
    }
    return x;
}

MAT &cholesky(MAT &A){
    int i,j,k;
    for (i=0; i<A.dim(); i++) {
        A[i][i]=sqrt(A[i][i]);
        for (j=i+1; j<A.dim(); j++) {
            A[j][i] /= A[i][i];
        }
        for (j=i+1; j<A.dim(); j++) {
            for (k=i+1; k<=j; k++)
                A[j][k] -= A[j][i]*A[k][i];
        }
    }
    return A;
} // Cholesky decomposition
VEC choSolve(MAT &L,VEC b){
    //forward sub.
    VEC y(b);
    for (int i = 0; i<L.dim(); i++) {
        y[i]/=L[i][i];
        for (int j = i+1; j<L.dim(); j++) {
            y[j] -= L[j][i]*y[i];
        }
    }
    MAT LT = L.tpose();
    return bckSubs(LT, y);
} // forward and backward substitutions

int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol){
    int iter = 0;
    double error = 1+tol;
    double s = 0;
    //VEC r = b - A*x;
    VEC x_new(x.len());
    while (error>tol&iter<maxIter) {
        iter++;
        for (int i = 0; i<A.dim(); i++) {
            s = 0;
            for (int j = 0; j<i; j++) {
                s+=A[i][j]*x[j];
            }
            for (int j = i+1; j<A.dim() ; j++) {
                s+=A[i][j]*x[j];
            }
            x_new[i] = (b[i]-s)/A[i][i];
        }
        error = (x_new-x).norm(A.p);
        x = x_new;
    }
    return iter;
}//Jacobi Iteration Method
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol){
    int iter = 0;
    double error = 1+tol;
    double s = 0;
    //VEC r = b - A*x;
    VEC x_old(x.len());
    
    while (error>tol&iter<maxIter) {
        iter++;
        for (int i = 0; i<A.dim(); i++) {
            s = 0;
            for (int j = 0; j<i; j++) {
                s+=A[i][j]*x[j];
            }
            for (int j = i+1; j<A.dim() ; j++) {
                s+=A[i][j]*x_old[j];
            }
            x[i] = (b[i]-s)/A[i][i];
        }
        error = (x-x_old).norm(A.p);
        x_old = x;
    }
    return iter;
} //Gauss-Seidel Iteration Method
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol){
    int iter = 0;
    double error = 1+tol;
    double s = 0;
    //VEC r = b - A*x;
    VEC x_old(x.len());
    
    while (error>tol&iter<maxIter) {
        iter++;
        VEC x_hold(x_old);
        for (int i = 0; i<A.dim(); i++) {
            s = 0;
            for (int j = 0; j<i; j++) {
                s+=A[i][j]*x[j];
            }
            for (int j = i+1; j<A.dim() ; j++) {
                s+=A[i][j]*x_old[j];
            }
            x[i] = (b[i]-s)/A[i][i];
        }
        x_old = x;
        for (int i = 0; i<A.dim(); i++) {
            s = 0;
            for (int j = 0; j<i; j++) {
                s+=A[i][j]*x_old[j];
            }
            for (int j = i+1; j<A.dim() ; j++) {
                s+=A[i][j]*x[j];
            }
            x[i] = (b[i]-s)/A[i][i];
        }
        error = (x-x_hold).norm(A.p);
        x_old = x;
    }
    return iter;
} //Symmetric Gauss-Seidel Iteration Method
int cg(MAT &A,VEC b,VEC &x,int maxIter, double tol){
    
    //algorithm initial parameter declare here
    int iter = 0;
    double error = 1+tol;
    VEC r = b - A*x;
    VEC p=r;
    //main parameter declare here
    double alpha;
    double beta;
    double n = x.len();
    
    //reuse temporal parameter declare here
    double rr_; //number for content r*r temporaly
    VEC Ap(p); //vector for content A*p temporaly

    //main algorithm:
    while (error>=tol&iter<maxIter) {
        rr_ = r*r;
        error = sqrt(rr_/n);
        Ap = A*p;
        alpha = (rr_)/(p*(Ap));
        x = x+(alpha*p);
        r = r-(alpha*Ap);
        beta = (r*r)/(rr_);
        p = r + beta*p;
        iter++;
        //cout<<error;
    }
    return iter;
}//Conjugate Gradient Method

int powerMethod(MAT &A,double &v,VEC &q,double maxIter,double tol){
    
    //initialize parameter
    for (int i = 0; i<q.len(); i++) {
        q[i]=1;
    }
    
    int iter = 0;
    double err = 1+tol;
    q = q/q.norm(2);
    VEC z(q.len());
    VEC r(q.len());
    VEC w(q.len());
    
    //Iteration
    while (err>=tol&&iter<=maxIter) {
        z = A*(q);
        iter++;
        q = z/z.norm(2);
        v = q*(A*q);
        r = A*q-v*q;
        //w=A.tpose()*q;
        w = q*A;
        w = w/w.norm(2);
        err = r.norm(2)/fabs(w*q);
    }
    
    return iter;
}//Power Method

int invPowerMethod(MAT &A,double &v,VEC &q,double maxIter,double tol){
    //initialize parameter
    for (int i = 0; i<q.len(); i++) {
        q[i]=1;
    }
    
    int iter = 0;
    double err = 1+tol;
    q = q/q.norm(2);
    VEC z(q.len());
    VEC r(q.len());
    VEC w(q.len());
    MAT LU(A);
    LU=luFact(LU);
    while (err>=tol&&iter<=maxIter) {
        z = bckSubs(LU,fwdSubs(LU, q));
        iter++;
        q = z/z.norm(2);
        v = q*(A*q);
        r = A*q-v*q;
        //w=A.tpose()*q;
        w = q*A;
        w = w/w.norm(2);
        err = r.norm(2)/fabs(w*q);
    }
    
    return iter;
}//Inverse Power Method

int sinvPowerMethod(MAT &A,double &v,VEC &q,double maxIter,double tol,double omega){
    //initialize parameter
    for (int i = 0; i<q.len(); i++) {
        q[i]=1;
    }
    
    int iter = 0;
    double err = 1+tol;
    q = q/q.norm(2);
    VEC z(q.len());
    VEC r(q.len());
    VEC w(q.len());
    MAT I = eye(q.len());
    MAT LU(A-omega*I);
    LU=luFact(LU);
    while (err>=tol&&iter<=maxIter) {
        z = bckSubs(LU,fwdSubs(LU, q));
        iter++;
        q = z/z.norm(2);
        v = q*(A*q);
        r = A*q-v*q;
        //w=A.tpose()*q;
        w = q*A;
        w = w/w.norm(2);
        err = r.norm(2)/fabs(w*q);
    }
    
    return iter;
}//Inverse Power Method with Shifting

void QRDecomposition(MAT &AA,MAT &Q,MAT &R){
    MAT A(AA.tpose());
    R[0][0] = sqrt(A[0]*A[0]);
    Q[0] = A[0]/R[0][0];
    VEC A_tmp(A.dim());
    for (int j = 1; j<A.dim(); j++) {
        for (int i = 0; i<=j; i++) {
            R[j][i] = Q[i]*A[j];
        }
        A_tmp = A[j];
        for (int i=0; i<=j-1; i++) {
            A_tmp-=R[j][i]*Q[i];
        }
        R[j][j]=sqrt(A_tmp*A_tmp);
        Q[j] = A_tmp/R[j][j];
    }
}//QR decomposition

int EVqr(MAT &A,double tol,int maxiter){
    int iter;
    
    
    
    return iter;
}//QR iteration for finding all eigen value and eigen vector

double EVqrError(MAT &A){
    double err = 0.0;
    for (int i = 1; i<A.dim(); i++) {
        if (fabs(A[i][i-1])>err) {
            err = fabs(A[i][i-1]);
        }
    }
    return err;
    
}//calculate the error for qr iteration

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
        os<<*v<<endl;
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
        //cout<<"x\n"<<x<<endl;
        //cout<<"r\n"<<r<<endl;
        rr_ = r*r;
        //std::cout<<"rr_"<<rr_<<std::endl;
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

void QRDecomposition(MAT &A_,MAT &Q,MAT &R){
    MAT A(A_);
    Q = MAT(A.dim());
    R = MAT(A.dim());
    R[0][0] = sqrt(A[0]*A[0]);
    Q[0] = A[0]/R[0][0];
    for (int j = 1; j<A.dim(); j++) {
        for (int i = 0; i<=j; i++) {
            R[j][i] = Q[i]*A[j];
            if(i<j)  A[j]-=R[j][i]*Q[i];
        }
//        A_tmp = *new VEC(A.dim());
//        for (int i=0; i<=j-1; i++) {
//            A_tmp+=R[j][i]*Q[i];
//        }
//        A_tmp = A[j]-A_tmp;
        
//        A_tmp = A[j];
//        for (int i=0; i<j; i++) {
//            A_tmp-=R[j][i]*Q[i];
//        }
        R[j][j]=sqrt(A[j]*A[j]);
        Q[j] = A[j]/R[j][j];
    }
    //Q = Q.tpose();
    //R = R.tpose();
}//QR decomposition

int EVqr(MAT &T,double tol,int maxiter){
    int iter=0;
    int dim = T.dim();
    MAT Q(dim);
    MAT R(dim);
    //float err = 1 + tol;
    QRDecomposition(T, Q, R);
    double err = EVqrError(T);
    while (err>tol&&iter<maxiter) {

        cout<<iter<<endl;
        cout<<err<<endl;
        iter++;
        T = (Q*R).tpose();
        QRDecomposition(T, Q, R);
        err = EVqrError(T);

    }

    return iter;
}//QR iteration for finding all eigen value and eigen vector
int EVqrShifted(MAT &T,double mu,double tol,int maxiter){
    int iter=0;
    int dim = T.dim();
    MAT I = eye(dim);
    MAT muI = I*mu;
    MAT Q(dim);
    MAT R(dim);
    //float err = 1 + tol;
    MAT T_muI = T-muI;
    QRDecomposition(T_muI, Q, R);
    double err = EVqrError(T);
    while (err>tol&&iter<maxiter) {
        iter++;
        cout<<iter<<endl<<err<<endl;
        T = (Q*R).tpose()+muI;
        T_muI = T-muI;
        QRDecomposition(T_muI, Q, R);
        err = EVqrError(T);
    }
    T = (Q*R).tpose()+muI;
    return iter;
    
}//QR iteration for finding all eigen value and eigen vector with shifting
double EVqrError(MAT &A){
    double err = 0.0;
    for (int i = 1; i<A.dim(); i++) {
        if (fabs(A[i][i-1])>=err) {
            err = fabs(A[i][i-1]);
        }
    }
    return err;
}//calculate the error for qr iteration
void splineM(int N,VEC &X,VEC &Y,VEC &M){
    MAT A(N);
    VEC d(N);
    double hi;
    double hio1;
    double hi_hio1;
    A[0][0]=2;
    for (int i = 1; i<N-1; i++) {
        hi = X[i]-X[i-1];
        hio1 = X[i+1]-X[i];
        hi_hio1 = hi+hio1;
        A[i][i-1]=hi/hi_hio1;
        A[i][i]=2;
        A[i][i+1]=hio1/hi_hio1;
        d[i] = 6/hi_hio1*((Y[i+1]-Y[i])/hio1-(Y[i]-Y[i-1])/hi);
    }
    A[N-1][N-1]=2;
    A = luBanded(A,2);
    M = bckSubs(A,fwdSubs(A, d));
}// generate spline momentum M
double spline(double x,int N,VEC &X,VEC &Y,VEC &M){
    double B;
    double A;
    double hi;
    double x_xi_1;
    double hi6;
    for (int i = 1; i<N; i++) {
        if (x>X[i-1]&&x<X[i]) {
            // interpolation here
            hi = X[i]-X[i-1];
            B = Y[i-1] - M[i-1]*hi*hi/6;
            A = (Y[i]-B)/hi - M[i]*hi/6;
            x_xi_1 = x - X[i-1];
            hi6 = 6*hi;
            return M[i-1]*pow(X[i]-x, 3)/hi6 + M[i]*pow(x_xi_1, 3)/hi6 + A*(x_xi_1) + B;
            break;
            
        }else if (x==X[i-1]){
            return Y[i-1];
            break;
        }
    }
    return Y[N-1];
    
    
}// spline interp at x
double integ(VEC &X,VEC &Y,int n){
    //set up the weights vector
    double result;
    VEC W = newtoncotesWeights(n);
    double h = (X[n]-X[0])/n;
    for (int i = 0; i<X.len()-1; i = i+n) {
        for (int j = 0; j<W.len(); j++) {
            result+=h*W[j]*Y[j+i];
        }
    }
    return result;
}// composite nth order Newton-Cotes integral


VEC newtoncotesWeights(int n){
    VEC result(n+1);
    for (int i =0; i<=n; i++) {
        //construct remain vector for calculating the divider and the coefficient of polynomial
        VEC remain_vec(n);
        int count = 0;
        for (int j = 0; j<=n; j++) {
            if (i!=j) {
                remain_vec[count] = -j;
                count++;
            }
        }
        //construct the divider of each weight
        int divider = 1;
        for (int j=0;j<n; j++) {
            divider*=(i+remain_vec[j]); //(0-1)(0-2)(0-3)
        }
        VEC coefficient = coefficientofNCWeights(remain_vec);
        double dividend = 0;
        for (int j=0; j<coefficient.len(); j++) {
            dividend+=pow(n, coefficient.len()-j)*coefficient[j]/(coefficient.len()-j);
        }
        result[i] = dividend/divider;
    }
    return result;
}//compute the weight vector of Newton-cotes
VEC coefficientofNCWeights(VEC &in){
    VEC result(in.len()+1);
    int n = in.len();
    for (int i = 0; i<=n; i++) {
        int sum = 0;
        //construct the vector contain the products of each choice of input vector
        int choice_number = chosensize(n, i);
        //cout<<choice_number<<endl;
        VEC products(choice_number);
        int tmp = 0;
        recursivechoice(products, in, i, 0, 1.0,tmp);
        for (int j = 0; j<choice_number; j++) {
            sum += int(products[j]);
        }
        result[i] = sum;
    }
    return result;
}//return all possible choice of specific size from a set of number
int chosensize(int m,int n){
    double result = 1;
    for (int i = m,j=n; i>m-n&&j>0; i--,j--) {
        result *= i;
        result /= j;
    }
    return int(result);
}

void recursivechoice(VEC &result,VEC &in,int chosen_size,int begin,double product_value,int &result_index){
    // initial:
    // result: the result product vector : all is 1
    // chosen size: the original chosen size
    // begin: 0
    // result_index = 0
    
    if (chosen_size == 0) {
        result[result_index] = product_value;
        result_index++;
        return;
    }else{
        int end = in.len() - chosen_size;
        for (int i = begin; i<=end; i++) {
            recursivechoice(result,in,chosen_size-1,i+1,product_value*in[i],result_index);
        }
    }
}


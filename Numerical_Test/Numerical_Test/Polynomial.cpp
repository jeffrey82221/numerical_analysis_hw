//
//  Polynomial.cpp
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/1.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "Polynomial.h"
#include <math.h>

ostream& operator<<(ostream& os, Polynomial& P){
    for (int i = P.n(); i>=2; i--) {
        os<<P.a(i)<<"*x^"<<i<<"+";
    }
    os<<P.a(1)<<"*x+";
    os<<P.a(0);
    return os;
}
Polynomial deflation(Polynomial &P,Complex &z,Complex &reminder){
    int n = P.n();
    Complex *b = new Complex[n]; //from n-1 to 0
    b[n-1] = P.a(n);
    for (int i = n-1; i>=1; i--) {
        b[i-1] = P.a(i)+b[i]*z;
    }
    reminder = P.a(0)+b[0]*z;
    return Polynomial(b,n-1);
}
Polynomial deflation(Polynomial &P,Complex p,Complex q,Complex &R,Complex &S){
    int n = P.n();
    Complex *b = new Complex[n-1]; //from n-2 to 0
    b[n-2] = P.a(n);
    b[n-3] = P.a(n-1)-p*b[n-2];
    for (int i = n-4; i>=0; i--) {
        b[i] = P.a(i+2)-p*b[i+1]-q*b[i+2];
    }
    R = P.a(1)-p*b[0]-q*b[1];
    S = P.a(0)-q*b[0];
    return Polynomial(b, n-2);
}
Polynomial Differentiate(Polynomial &P){
    int n = P.n();
    Complex *b = new Complex[n]; //from n-1 to 0
    for (int i = 0; i<n; i++) {
        b[i] = P.a(i+1)*(i+1);
    }
    return Polynomial(b,n-1);
}
Polynomial newtonsMethod(Polynomial &P,Complex &x,double e,int maxIter){
    double k = 0;
    double err = 1+e;
    Complex f;
    Complex f_;
    Polynomial F=deflation(P,x,f);
    deflation(F, x, f_);
    while (err>=e && k<maxIter) {
        x = x - f/f_;
        //err = fabs(P.eval(x));
        err = fabs(f);
        k++;
        F = deflation(P, x, f);
        deflation(F, x, f_);
    }
    cout<<"k = "<<k<<endl;
    return F;
}
Complex* roots(Polynomial &P,double e,double maxIter){
    Complex *roots = new Complex[P.n()];
    //find the first guess 1+mu
    double max = -INFINITY;
    for (int i = 0; i<=P.n(); i++) {
        if (fabs(P.a(i))>max) {
            max = fabs(P.a(i));
        }
    }
    Complex x(max+1,0);
    Polynomial F = P;
    for (int i = 0; i<P.n(); i++) {
        F = newtonsMethod(F, x, e, maxIter);
        roots[i] = x;
    }
    return roots;
}
Polynomial quadraticMethod(Polynomial &P,Complex &p,Complex &q,double e,int maxIter){
    Complex R,S;
    Polynomial F = deflation(P, p, q, R, S);
    int k = 0;
    double err = 1+e;
    while (err>=e && k<maxIter) {
        p = p + R/F.a(0);
        q = q + S/F.a(0);
        F = deflation(P, p,q,R,S);
        err = fabs(R)>fabs(S)?fabs(R):fabs(S);
        k++;
    }
    return F;
}
Polynomial bairstowsMethod(Polynomial &P,Complex &p,Complex &q,double e,int maxIter){
    double err = 1+e;
    int k = 0;
    int n = P.n();
    Complex *b = new Complex[n-1];
    Complex *c = new Complex[n-1];
    Complex *d = new Complex[n-1];
    Complex R,S;
    Complex dRdp,dRdq,dSdp,dSdq;
    Complex div;
    while (err>=e&&k<maxIter) {
        b[n-2] = P.a(n);
        b[n-3] = P.a(n-1) - p*b[n-2];
        for (int j = n-4; j>=0; j--) {
            b[j] = P.a(j+2)-p*b[j+1]-q*b[j+2];
        }
        R = P.a(1)-p*b[0]-q*b[1];
        S = P.a(0) - q*b[0];
        c[n-2] = 0;
        c[n-3] = -b[n-2]-p*c[n-2];
        for (int j = n-4; j>=0; j--) {
            c[j] = -b[j+1] - p*c[j+1] - q*c[j+2];
        }
        dRdp = -b[0] - p*c[0]-q*c[1];
        dSdp = -q*c[0];
        d[n-2] = 0;
        d[n-3] = -p*d[n-2];
        for (int j = n-4; j>=0; j--) {
            d[j] = -p*d[j+1] - b[j+2] - q*d[j+2];
        }
        dRdq = -p*d[0] - b[1]-q*d[1];
        dSdq = -b[0] - q*d[0];
        div = (dRdp)*(dSdq)-(dRdq)*(dSdp);
        p = p - (R*dSdq-S*dRdq)/div;
        q = q - (-R*dSdp+S*dRdp)/div;
        k++;
        err = fabs(R)>fabs(S)?fabs(R):fabs(S);
    }
    return deflation(P, p, q, R, S);
}
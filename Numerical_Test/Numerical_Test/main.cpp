//
//  main.cpp
//  Numerical_Test
//
//  Created by 100061212 林奕勳 on 2015/4/26.
//  Copyright (c) 2015年 100061212 林奕勳. All rights reserved.
//
#define LEN 301;
#include <iostream>
#include "VEC.h"
#include "MAT.h"
#include <iomanip>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double L1 = 1e-12;
double L2 = 1e-9;
double R1 = 1;
double R2 = 1;
VEC zero2(2);
double V(double t){
    if (t>=0) {
        return 1;
    }else return 0;
}
double f1(VEC &i, double V){
    double i3 = i[0]-i[1];
    //if ((i3)*R1>V) {
      //  return V/L1;
    //}else{
        return (V-R1*(i3))/L1;
    //}
        //return -R2/L2*i[0] + R1/L2*i[1];
}
double f2(VEC &i,double V){
    return (R1*i[0]-(R1+R2)*i[1])/L2;
    //return R2/L2*i[0] - R1*(1/L1+1/L2)*i[1] + V/L1;
}
VEC f(VEC &i ,double V){
    VEC f(2);
    f[0] = f1(i, V);
    f[1] = f2(i, V);
    return f;
}
VEC xnp12ndOrderForwardIntegration(VEC &xn,VEC fn,VEC fn_1,double h){
    VEC xnp1(2);
    //cout<<fn<<","<<fn_1<<endl;
    xnp1 = (xn + h*(3*fn-fn_1)/2);
    return xnp1;
}

////////////////// for trapezoidal method ////////////////////////////////////////////
// for calcualted the left term
double F1(VEC i,double h,double V){
    return i[0] - h/2*f1(i, V);
}
double F2(VEC i,double h,double V){
    return i[1] - h/2*f2(i, V);
}
MAT Jacobian(double h){
    MAT r(2);
    VEC i1(2);
    VEC i2(2);
    VEC i0(2);
    i1[0] = 1;
    i2[1] = 1;
    r[0][0] = F1(i1,h, V(1))-F1(i0,h,V(1));
    r[0][1] = F1(i2,h, V(1))-F1(i0,h,V(1));
    r[1][0] = F2(i1,h, V(1))-F2(i0,h,V(1));
    r[1][1] = F2(i2,h, V(1))-F2(i0,h,V(1));
    return r;
}
VEC Fconst(double h,double V){
    VEC r(2);
    r[0] = F1(zero2, h, V);
    r[1] = F2(zero2, h, V);
    return r;
}
VEC trapezoidalRule(VEC in,double h,double t){
    MAT J = Jacobian(h);
    VEC inp1(2);
    VEC b = in + h/2*f(in, V(t)) - Fconst(h, V(t+h));
    luFact(J);
    inp1 = bckSubs(J,  fwdSubs(J, b));
    return inp1;
}

// 3rd order forward //////////////////////////////////////////////////////////////////////
VEC xnp13rdOrderForwardIntegration(VEC &xn,VEC fn,VEC fn_1,VEC fn_2,double h){
    VEC xnp1(2);
    //cout<<fn<<","<<fn_1<<endl;
    xnp1 = (xn + h*(23*fn-16*fn_1+5*fn_2)/12);
    return xnp1;
}
//// general approach for backward integration and forward integration/////////////////////
// backward integration :
VEC b_bckward(int k){ // generate the coefficients b0 ... bk-1 for k'th order
    MAT M(k);
    VEC b(k);
    VEC y(k);
    // contructed the formula matrix
    for (int i = 0; i<k; i++) {
        for (int j = 0; j<k; j++) {
            if(i==0)    M[i][j] = 1;
            else    M[i][j] = M[i-1][j]*(1-j);
        }
    }
    // contructed the right hand side vec
    for (int i = 0; i<k; i++) {
        y[i] = 1.0/(i+1);
    }
    luFact(M);
    b = bckSubs(M,  fwdSubs(M, y));
    return b;
}
VEC F(VEC x1, double h,double b0,double t){ // the left hand side term as an vector function
    return x1 - h*b0*f(x1,V(t+h)); // here we call the system formula f , which contain the formulas of the circuit
}
MAT JF(double h,int rank,double b0){ // Jacobian of F , finding the slope of the system in all direction as system matrix
    MAT J(rank);
    VEC zeros(rank);
    MAT I = eye(rank);
    for (int i = 0; i<rank; i++) {
        J[i] = F(I[i], h, b0, 0) - F(zeros, h, b0, 0);
    }
    return J.tpose();
}
VEC F0(double h,int rank,double b0,double t){ // the constant of the system is add as F with input 0_
    VEC zeros(rank);
    return F(zeros, h, b0, t);
}
VEC backwardMethod(VEC **x_,double h,double t, int order,int rank){
    VEC x(rank);
    VEC b = b_bckward(order);
    MAT J = JF(h, rank, b[0]);
    VEC y(rank);
    for (int i = 0; i<order-1; i++) {
        y+=b[i+1]*f(*x_[i], V(t-i*h));
    }
    y *= h;
    y += *x_[0];
    y -= F0(h, rank, b[0], t+h);
    //cout<<y<<endl;
    luFact(J);
    x = bckSubs(J,  fwdSubs(J, y));
    //update the new x_
    for (int i = order-2; i>=1; i--) {
        *x_[i] = *x_[i-1];
    }
    *x_[0] = x;
    return x;
}

// forward integration :
VEC b_forward(int k){ // generate the coefficients b0 ... bk-1 for k'th order
    MAT M(k);
    VEC b(k);
    VEC y(k);
    // contructed the formula matrix
    for (int i = 0; i<k; i++) {
        for (int j = 0; j<k; j++) {
            if(i==0)    M[i][j] = 1;
            else    M[i][j] = M[i-1][j]*(-j);
        }
    }
    // contructed the right hand side vec
    for (int i = 0; i<k; i++) {
        y[i] = 1.0/(i+1);
    }
    luFact(M);
    b = bckSubs(M,  fwdSubs(M, y));
    return b;
}

VEC forwardMethod(VEC **x_,double h,double t, int order,int rank){
    VEC x(rank);
    VEC b = b_forward(order);
    for(int i = 0;i<order;i++){
        x+= b[i]*f(*x_[i], V(t-i*h));
    }
    x*=h;
    x+=*x_[0];
    //update the new x_
    for (int i = order-1; i>=1; i--) {
        *x_[i] = *x_[i-1];
    }
    *x_[0] = x;
    return x;
}
///// Gear's Method ////////////////////////////////////////////////////
VEC a_gears(int k){ // generate the coefficients b0 ... bk-1 for k'th order
    MAT M(k+1);
    VEC b(k+1);
    VEC y(k+1);
    // contructed the formula matrix
    for (int i = 0; i<k+1; i++) {
        for (int j = 0; j<k; j++) {
            if(i==0)    M[i][j] = 1;
            else    M[i][j] = M[i-1][j]*(-j);
        }
    }
    for (int i = 0; i<k+1; i++) {
        M[i][k] = i;
    }
    // contructed the right hand side vec
    for (int i = 0; i<k+1; i++) {
        y[i] = 1.0;
    }
    luFact(M);
    b = bckSubs(M,  fwdSubs(M, y));
    return b;
}
VEC gearsMethod(VEC **x_,double h,double t, int order,int rank){
    VEC x(rank);
    VEC a = a_gears(order);
    MAT J = JF(h, rank, a[order]);
    VEC y(rank);
    for (int i = 0; i<order; i++) {
        y+=a[i]**x_[i];
    }
    y -= F0(h, rank, a[order], t+h);
    luFact(J);
    x = bckSubs(J,  fwdSubs(J, y));
    //update the new x_
    for (int i = order-1; i>=1; i--) {
        *x_[i] = *x_[i-1];
    }
    *x_[0] = x;
    return x;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) {
    /// trapezoidal Rule
    cout << "trapezoidal Rule : " <<endl;
    auto h = 1e-13;
    double top = 5e-8;
    int count = 0;
    int mod = 10000;
    VEC in(2);  //i(0) = (0,0)
    for (int i = 0; i<3; i++) {
        in = zero2;
        count = 0;
        cout<<"with h = "<<h<<" : "<<endl;
        for (double t = 0.; t<=top+h; t+=h) {
            if (count%mod==0) {
                cout<<"t = "<<t<<", i1 = "<<scientific<<setprecision(11)<<in[0]<<", i2 = "<<in[1]<<endl;
                //cout<<setprecision( 15 )<<in[1]<<endl;
            }
            in = trapezoidalRule(in, h, t);
            count++;
        }
        h = h*100;
        mod = mod/100;
    }
    cout<<"2nd order forward integration : "<<endl;
    h = 1e-11;
    mod = 100;
    cout<<"with h = "<<h<<" : "<<endl;
    VEC **i;
    int order = 2;
    int rank = 2;
    count = 0;
    i=(VEC **)malloc((order)*sizeof(VEC*));
    for (int ii = 0; ii<order; ii++) {
        i[ii] = newVEC(rank);
    }
    for (double t = 0; t<=top+h; t+=h) {
        //cout<<setprecision( 11 )<<"t = "<<t<<" i = "<<*i[0]<<endl;
        if (count%mod==0) {
            cout<<"t = "<<t<<", i1 = "<<scientific<<(*i[0])[0]<<", i2 = "<<(*i[0])[1]<<endl;
        }
        forwardMethod(i, h, t, order, rank);
        count++;
    }
    h = 1e-9;
    mod = 1;
    cout<<"with h = "<<h<<" : "<<endl;
    count = 0;
    i=(VEC **)malloc((order)*sizeof(VEC*));
    for (int ii = 0; ii<order; ii++) {
        i[ii] = newVEC(rank);
    }
    for (double t = 0; t<=top+h; t+=h) {
        //cout<<setprecision( 11 )<<"t = "<<t<<" i = "<<*i[0]<<endl;
        if (count%mod==0) {
            cout<<"t = "<<t<<", i1 = "<<scientific<<(*i[0])[0]<<", i2 = "<<(*i[0])[1]<<endl;
        }
        forwardMethod(i, h, t, order, rank);
        count++;
    }
    cout<<"2nd order gears integration : "<<endl;
    h = 1e-11;
    mod = 100;
    cout<<"with h = "<<h<<" : "<<endl;
    count = 0;
    i=(VEC **)malloc((order)*sizeof(VEC*));
    for (int ii = 0; ii<order; ii++) {
        i[ii] = newVEC(rank);
    }
    for (double t = 0; t<=top+h; t+=h) {
        //cout<<setprecision( 11 )<<"t = "<<t<<" i = "<<*i[0]<<endl;
        if (count%mod==0) {
            cout<<"t = "<<t<<", i1 = "<<scientific<<(*i[0])[0]<<", i2 = "<<(*i[0])[1]<<endl;
        }
        gearsMethod(i, h, t, order, rank);
        count++;
    }
    h = 1e-9;
    mod = 1;
    cout<<"with h = "<<h<<" : "<<endl;
    count = 0;
    i=(VEC **)malloc((order)*sizeof(VEC*));
    for (int ii = 0; ii<order; ii++) {
        i[ii] = newVEC(rank);
    }
    for (double t = 0; t<=top+h; t+=h) {
        //cout<<setprecision( 11 )<<"t = "<<t<<" i = "<<*i[0]<<endl;
        if (count%mod==0) {
            cout<<"t = "<<t<<", i1 = "<<scientific<<(*i[0])[0]<<", i2 = "<<(*i[0])[1]<<endl;
        }
        gearsMethod(i, h, t, order, rank);
        count++;
    }
    return 0;
}
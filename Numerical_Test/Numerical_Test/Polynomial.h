//
//  Polynomial.h
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/1.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#ifndef __Numerical_Test__Polynomial__
#define __Numerical_Test__Polynomial__
#include "Complex.h"
#include <stdio.h>
#include <iostream>
class Polynomial{
private:
    Complex *coefficients;
    int order;
public:
    Polynomial(double co[],int order){
        this->order = order;
        this->coefficients = new Complex[order+1];
        for (int i = 0; i<=order; i++) {
            this->coefficients[i]+= co[i];
        }
    }
    Polynomial(Complex co[],int order){
        this->order = order;
        this->coefficients = co;
    }
    
    Complex a(int i){
        return coefficients[i];
    }
    int n(){
        return order;
    }
    Complex eval(Complex x){
        Complex result;
        for (int i = 0; i<=order; i++) {
            result=result+(coefficients[i]*(x^i));
        }
        return result;
    }
    
};

ostream& operator<<(ostream& os, Polynomial& P);
Polynomial deflation(Polynomial &P,Complex &root,Complex &reminder);
Polynomial deflation(Polynomial &P,Complex p,Complex q,Complex &R,Complex &S); //for quadratic method
Polynomial Differentiate(Polynomial &P);
Polynomial newtonsMethod(Polynomial &P,Complex &in,double e);
Complex* roots(Polynomial &P,double e,double maxIter);
Polynomial quadraticMethod(Polynomial &P,Complex &p,Complex &q,double e,int maxIter);
Polynomial bairstowsMethod(Polynomial &P,Complex &p,Complex &q,double e,int maxIter);

#endif /* defined(__Numerical_Test__Polynomial__) */


//TEST::
//    //double y = bisectionMethod(0.01,2.5,0.0000000001);
//    int order = 4;
//    double co[5] =  {1.,1.,2.,1.,1.};
//    Polynomial P1(co,order);
//    cout<<"F(x) = "<<P1<<"=0"<<endl;
////    Complex p(0,0),q(0,0);
////    Polynomial P2 = quadraticMethod(P1,p,q,0.00000000001,1000000);
////    cout<<"P2 = "<<P2<<endl;
////    cout<<"x^2+"<<p<<"x+"<<q<<endl;
////    Complex *rs = roots(P1,1e-17,10000000);
////    cout<<"roots = "<<endl;
////    cout<<fixed<<setprecision( 10 );
////    for (int i = 0; i<P1.n(); i++) {
////        cout<<"z"<<i<<" = "<<rs[i]<<endl;
////    }
////
////    cout<<"check:"<<endl;
////    for (int i = 0; i<P1.n(); i++) {
////        cout<<"F("<<"z"<<i<<") = "<<P1.eval(rs[i])<<endl;
////    }
//    Complex p(1,0),q(2,0);
//    Polynomial P2 = bairstowsMethod(P1, p, q,0.0000000000001, 1000000000);
//    cout<<"x^2+"<<p<<"x+"<<q<<endl;
//    cout<<P2<<endl;
////    cout<<"P1 = "<<P1<<endl;
////    Complex b_1;
////    Complex z10(12,0);
////    Complex reminder;
////    Polynomial P2 = newtonsMethod(P1,z10,0.000000001,10000000);
////    cout<<"root = "<<z10<<endl;
////    cout<<"P2 = "<<P2<<endl;
////    Complex x0(0,0.1);
////    Complex z1 = newtonsMethod(x0,0.0000000001);
////    cout<<fixed<<setprecision( 7 ) <<"z1 = "<<z1<<endl;
////    cout<<"F(z1) = "<<P1.eval(z1)<<endl;
////    Polynomial P2 = deflation(P1, z1);
////    cout<<"P2 = "<<P2<<endl;
////    Complex z2 = newtonsMethod(Complex(0.1,-1.2),0.0000000001);
////    cout<<"z2 = "<<z2<<endl;
////    cout<<"F(z2) = "<<Func(z2)<<endl;
////    Polynomial P3 = deflation(P2, z2);
////    cout<<"P3 = "<< P3<<endl;
//    //check


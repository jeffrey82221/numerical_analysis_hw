//
//  Complex.cpp
//  Numerical_Analysis
//
//  Created by jeffrey_MAC on 2015/3/4.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "Complex.h"
#include <math.h>
Complex::Complex(double r,double i){
    this->x = r;
    this->y = i;
}
Complex::Complex(const Complex &C){ //copy constructor
    x = C.r();
    y = C.i();
}
double Complex::r() const{
    return x;
}//get real part
double Complex::i() const{
    return y;
}//get imaginary part
Complex& Complex::operator+= (Complex &z){
    x+=z.r();
    y+=z.i();
    return *this;
}
Complex& Complex::operator+= (double dbl){
    x+=dbl;
    return *this;
}// C1 += dbl;
Complex& Complex::operator-= (Complex &z){
    x-=z.r();
    y-=z.i();
    return *this;
}// C1 -= C2;
Complex& Complex::operator-= (double dbl){
    x-=dbl;
    return *this;
}// C1 -= dbl;
Complex& Complex::operator*= (Complex &z){
    double xx = x;
    double yy = y;
    x = xx*z.x-yy*z.y;
    y = yy*z.x+xx*z.y;
    return *this;
}
Complex& Complex::operator*= (double dbl){
    x*=dbl;
    y*=dbl;
    return *this;
}// C1 += dbl;
Complex& Complex::operator/= (Complex &z){
    double xx = x;
    double yy = y;
    x = xx*z.x+yy*z.y;
    y = yy*z.x-xx*z.y;
    double div =(z.x*z.x+z.y*z.y);
    x/=div;
    y/=div;
    return *this;
}// C1 -= C2;
Complex& Complex::operator/= (double dbl){
    x/=dbl;
    y/=dbl;
    return *this;
}// C1 -= dbl;
//non-member functions

ostream& operator<<(ostream& os, const Complex& C){
    if (C.i()>0) {
        os<<C.r()<<"+"<<C.i()<<"i";
    }else if(C.i()<0){
        os<<C.r()<<C.i()<<"i";
    }else{
        os<<C.r();
    }
    return os;
}

Complex operator+(Complex z){
    Complex z1(z); //copy to another z
    return z1;
}
Complex operator+(Complex z1,Complex z2){
    Complex z(z1);
    return z+=z2;
}
Complex operator+(Complex z1,double d){
    Complex z(z1);
    return z+=d;
}
Complex operator+(double d,Complex z1)
{
    Complex z(z1);
    return z+=d;
}
//minus
Complex operator-(Complex z){
    Complex z1 = Complex(-z.r(),-z.i()); //copy to another z
    return z1;
}
Complex operator-(Complex z1,Complex z2){
    Complex z(z1);
    return z-=z2;
}
Complex operator-(Complex z1,double d){
    Complex z(z1);
    return z-=d;
}
Complex operator-(double d,Complex z1)
{
    Complex z(z1);
    return z-=d;
}
//multiply
Complex operator*(Complex z1,Complex z2){
    Complex z(z1);
    return z*=z2;
}
Complex operator*(Complex z1,double d){
    Complex z(z1);
    return z*=d;
}
Complex operator*(double d,Complex z1)
{
    Complex z(z1);
    return z*=d;
}
//divide
Complex operator/(Complex z1,Complex z2){
    Complex z(z1);
    return z/=z2;
}
Complex operator/(Complex z1,double d){
    Complex z(z1);
    return z/=d;
}
Complex operator/(double d,Complex z1)
{
    Complex z(z1);
    return z/=d;
}


double fabs(Complex z)
{
    return sqrt(z.r()*z.r()+z.i()*z.i());
}
double arg(Complex z)
{
    double theta = atan(z.i()/z.r());
    if (z.r()>0) {
        return theta;
    }else if(z.r()<0&&z.i()>=0){
        return theta+M_PI;
    }else if(z.r()<0&&z.i()<0){
        return theta-M_PI;
    }else if(z.r()==0&&z.i()>0){
        return M_PI/2;
    }else if(z.r()==0&&z.i()<0){
        return -M_PI/2;
    }else{// if(z.r()==0&&z.i()==0){
        return NAN;
    }
}
Complex &operator++(Complex &z1)
{
    return z1+=1.0;
}
Complex operator++(Complex &z1,int)
{
    Complex z(z1);
    ++z1;
    return z;
}
int operator==(Complex z1,Complex z2){
    if (z1.r()==z2.r()&&z1.i()==z2.i()) return 1;
    return 0;
}
int operator!=(Complex z1,Complex z2){
    if (z1.r()!=z2.r()||z1.i()!=z2.i()) return 1;
    return 0;
}
Complex operator^(Complex z1,int r){
    
    Complex z(1,0);
    if(r>=0){
        for (int i=0; i<r; i++) {
            z*=z1;
        }
    }else{
        for (int i=0; i<-r; i++) {
            z/=z1;
        }
    }
    return z;
}
Complex cong(Complex z1){
    Complex z(z1.r(),-z1.i());
    return z;
}
Complex sqrt(Complex z1){
    double r = fabs(z1);
    return sqrt(r)*((z1+r)/fabs(z1+r));
}
Complex log(Complex z1){
    double r = log(fabs(z1));
    double i = arg(z1);
    return Complex(r,i);
    
}





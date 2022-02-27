//
//  Complex.h
//  Numerical_Analysis
//
//  Created by jeffrey_MAC on 2015/3/3.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#ifndef Numerical_Analysis_Complex_h
#define Numerical_Analysis_Complex_h
#include <iostream>
using namespace std;

class Complex {
public:
    Complex(double r=0,double i=0); //initial constructor
    Complex(const Complex &C); //copy constructor
    double r() const; //get real part
    double i() const; //get imaginary part
    Complex& operator+= (Complex&);
    Complex& operator+= (double);       // C1 += dbl;
    Complex& operator-= (Complex&);  // C1 -= C2;
    Complex& operator-= (double);       // C1 -= dbl;
    Complex& operator*= (Complex&);  // C1 *= C2;
    Complex& operator*= (double);       // C1 *= dbl;
    Complex& operator/= (Complex&);  // C1 /= C2;
    Complex& operator/= (double);       // C1 /= dbl;
    friend ostream& operator<<(ostream&, const Complex&);
private:
    double x,y;
};

//non-member functions
Complex operator+(Complex);         //unary plus
Complex operator+(Complex,Complex); //plus
Complex operator+(Complex,double);  //plus
Complex operator+(double,Complex);  //plus
Complex operator-(Complex);         //unary minus
Complex operator-(Complex,Complex); //minus
Complex operator-(Complex,double);  //minus
Complex operator-(double,Complex);  //minus
//multiply
Complex operator*(Complex,Complex);
Complex operator*(Complex,double);
Complex operator*(double,Complex);
//divide
Complex operator/(Complex,Complex);
Complex operator/(Complex,double);
Complex operator/(double,Complex);
double fabs(Complex);
double arg(Complex);
Complex& operator++(Complex&);      //prefix increment
Complex operator++(Complex&,int);   //postfix increment
int operator==(Complex,Complex);    //equal
int operator!=(Complex,Complex);    //not equal
Complex operator^(Complex,int);
Complex cong(Complex); //congugate
Complex sqrt(Complex);
Complex log(Complex);


#endif

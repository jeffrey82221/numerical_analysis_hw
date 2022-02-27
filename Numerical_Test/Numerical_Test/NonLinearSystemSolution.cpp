//
//  NonLinearSystemSolution.cpp
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/1.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "NonLinearSystemSolution.h"
#include <math.h>
#include "Complex.h"

double Func(double x){
    //return pow(x,1.0/3.0)+pow(x,1.0/2.0)-1.0;
    //return pow(log(x),3.0)+log(x)-1.0;
    //return pow(tan(x),3.0)+tan(x)-1.0;
    //return 1.0/(pow(x,1.0/2.0)+1)+1.0/(pow(x,1.0/2.0)+2)-1.0;
    //return exp(x)-pow(log(x),2.0);
    //return (tan(x)+1)/(pow(tan(x),2.0)+1)-x;
    //return pow(x,3.0)+pow(x,1.0/2.0)-10.0;
    return -24.0 + 94.0*x -120.0*pow(x,2.0) + 62.5*pow(x,3.0)-13.5*pow(x,4.0)+pow(x,5.0);
    
}
double Func_(double x){
    return 94.0 - 2.*120.*pow(x, 1.0) + 3*62.5*pow(x, 2.0) - 4*13.5*pow(x, 3.0) + 5*pow(x, 4.0);
}
Complex Func(Complex x){
    return (x^2)+1;
}
Complex Func_(Complex x){
    return 2*x;
}


double bisectionMethod(double a,double b,double e){
    double x = (a+b)/2.0;
    int k = 0;
    while(fabs(a-x)>e){
        if(Func(a)*Func(x)<=0){
            a = a;
            b = x;
        }else{
            a = x;
            b = b;
        }
        k++;
        x = (a+b)/2.0;
        
    }
    return x;
    
}

double chordMethod(double a,double b,double e){
    double g = (Func(b)-Func(a))/(b-a);
    double x = b;
    double k = 0;
    double err = 1+e;
    while(err>e){
        x = x - Func(x)/g;
        k++;
        err = fabs(Func(x));
    }
    return x;
}
double regulaFalsiMethod(double a,double b,double e){
    double k = 0;
    double err = 1+e;
    double x=0;
    while(err>e){
        x = a - Func(a)*(b-a)/(Func(b)-Func(a));
        if(Func(x)*Func(a)<=0){
            b = x;
        }else{
            a = x;
        }
        k++;
        err = fabs(Func(x));
        
    }
    return x;
}
double secantMethod(double a,double b,double e){
    double k = 0;
    double err = 1 + e;
    double x_ = a; //xk-1
    double x = b;  //xk
    double x_new = 0;//xk+1
    while (err>e) {
        x_new = x - Func(x)*(x-x_)/(Func(x)-Func(x_));
        k++;
        x_ = x;
        x = x_new;
        err = fabs(Func(x));
    }
    return x;
}
double newtonsMethod(double x,double e){
    double k = 0;
    double err = 1+e;
    while (err>e) {
        x = x - Func(x)/Func_(x);
        k++;
        err = fabs(Func(x));
    }
    return x;
}
Complex newtonsMethod(Complex x,double e){
    double k = 0;
    double err = 1+e;
    while (err>e) {
        x = x - Func(x)/Func_(x);
        k++;
        err = fabs(Func(x));
    }
    return x;
}
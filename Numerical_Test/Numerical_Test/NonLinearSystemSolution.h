//
//  NonLinearSystemSolution.h
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/1.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#ifndef __Numerical_Test__NonLinearSystemSolution__
#define __Numerical_Test__NonLinearSystemSolution__

#include <stdio.h>
#include "Complex.h"
double Func(double x);
double Func_(double x);
Complex Func(Complex x);
Complex Func_(Complex x);
double bisectionMethod(double a,double b,double e);
double chordMethod(double a,double b,double e);
double regulaFalsiMethod(double a,double b,double e);
double secantMethod(double a,double b,double e);
double newtonsMethod(double x,double e);
Complex newtonsMethod(Complex x,double e);
#endif /* defined(__Numerical_Test__NonLinearSystemSolution__) */

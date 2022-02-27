//
//  nonlinearCircuit.h
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/9.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//  This Class implement the nonlinear resistor network depend only on voltage of HW12

#ifndef __Numerical_Test__nonlinearCircuit__
#define __Numerical_Test__nonlinearCircuit__
#include "Circuit.h"
#include <stdio.h>
class nonlinearCircuit : public Circuit{
private:
    double R0;
    double alpha;
    //non-linear resistor equation
public:
    //non-linear resistor equation
    virtual double I(double V);
    //differentiation of non-linear resistor equation
    virtual double dIdV(double V);
    void setParameter(double R0,double alpha);
};

#endif /* defined(__Numerical_Test__nonLinearCircuit__) */

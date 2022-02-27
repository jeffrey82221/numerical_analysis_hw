//
//  TempCircuit.h
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/9.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//  This Class implement the nonlinear resistor network depend on temperature of HW12

#ifndef __Numerical_Test__TempCircuit__
#define __Numerical_Test__TempCircuit__
#include "Circuit.h"
#include <stdio.h>
class TempCircuit : public Circuit{
private:
    double R0;
    double k;
    double beta;
    //non-linear resistor equation
protected:
public:
    //non-linear resistor equation
    virtual double I(double V);
    //differentiation of non-linear resistor equation
    virtual double dIdV(double V);
    void setParameter(double R0,double k,double beta);
    double temperature(int node1,int node2);
};
#endif /* defined(__Numerical_Test__TempCircuit__) */

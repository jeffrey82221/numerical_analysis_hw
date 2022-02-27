//
//  TempCircuit.cpp
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/9.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "TempCircuit.h"
double TempCircuit::I(double V){
    //return (-R0+sqrt(R0*R0+4*V*V*k*beta))/(2*k*beta*V);
    //return (-1+sqrt(1+4*V*V))/(2*V);
    return 2*V/(R0+sqrt(R0*R0+4*k*beta*V*V));
}
//differentiation of non-linear resistor equation
double TempCircuit::dIdV(double V){
    //return (R0/(2*k*beta*V*V))*(1-R0/sqrt(R0*R0+4*k*beta*V*V));
    //return (1/(2*V*V))*(1-1/sqrt(1+4*V*V));
    double tmp = 1+4*k*beta*pow(V/R0,2);
    return 2/R0/(tmp+sqrt(tmp));
}
void TempCircuit::setParameter(double R0,double k,double beta){
    this->R0 = R0;
    this->k = k;
    this->beta = beta;
}

double TempCircuit::temperature(int node1,int node2){
    if (node1==node2||(*this->adjecentM)[node1-1][node2-1]==0.) {
        cout<<"no divice between node "<<node1<<" and node "<<node2 <<endl;
        return NAN;
    }
    return (-R0+sqrt(R0*R0+4*k*beta*pow(((*voltages)[node1-1]-(*voltages)[node2-1]),2)))/(2*k);
}

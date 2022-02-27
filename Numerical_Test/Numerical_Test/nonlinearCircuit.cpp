//
//  nonlinearCircuit.cpp
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/9.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "nonlinearCircuit.h"
double nonlinearCircuit::I(double V){
    return V/(R0+alpha*fabs(V));
}
double nonlinearCircuit::dIdV(double V){
    return R0/pow(R0+alpha*fabs(V), 2);
}

void nonlinearCircuit::setParameter(double R0,double alpha){
    this->R0 = R0;
    this->alpha = alpha;
}

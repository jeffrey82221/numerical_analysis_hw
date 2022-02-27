//
//  Circuit.h
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/3.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#ifndef __Numerical_Test__Circuit__
#define __Numerical_Test__Circuit__

#include <stdio.h>
#include "MAT.h"
#include "VEC.h"
#include <math.h>
#include "CircuitConfiguration.h"
class Circuit{
protected:
    int node_num;
    int vdd_index;
    int gnd_index;
    MAT *adjecentM;
    double VDD;
    VEC *voltages;
public:
    virtual double I(double V);
    //differentiation of non-linear resistor equation
    virtual double dIdV(double V);
    //constructor
    Circuit();
    void setCircuit(int num,int vdd,int gnd,MAT &adjecentM);
    void resetVoltage();
    void setVoltage(VEC &v){ this->voltages = &v;}
    void setVDD(double vdd){ VDD = vdd;}
    //non-linear function for each node //node = 1~9
    double F(int node,VEC &Vs);
    //build the dFi(x)/dxj Jacobian function //fi = 1~9 //vj = 1~9
    double dFdV(int fi,int vj,VEC &vs);
    //Fs
    VEC Fs(VEC &vs);
    virtual MAT Jacobian(VEC &vs);
    //solve the non-linear system by newtons method
    VEC eval(double e,int maxIter);
    double totalCurrentSupply();
    double current(int node1,int node2);
    
};




#endif /* defined(__Numerical_Test__Circuit__) */


/*
 //read in the file for comparision
 MAT configuration = adjMatrix(3);
 srand((unsigned)time(NULL));
 nonlinearCircuit circuit1;
 circuit1.setParameter(1,0.1);
 circuit1.setCircuit(9,1,8,configuration);
 //circuit1.setVDD(2);
 cout<<"voltage dependent resistance : "<<endl;
 //circuit1.resetVoltage();
 VEC v(9);
 for (double vdd = 0; vdd<=5; vdd+=0.0001) {
 circuit1.setVDD(vdd);
 //circuit1.resetVoltage();
 v = circuit1.eval(0.00000001, 10000000);
 cout<<fixed << setw( 2 ) << setprecision( 4 ) <<" for vdd = "<<vdd<<" (Volt) , total current supplied = "
 <<fixed << setw( 7 ) << setprecision( 6 ) <<circuit1.totalCurrentSupply()<<"  , ";
 cout<<"I(r2) = "<<circuit1.current(4, 7)<<" , I(r7) = "<<circuit1.current(5, 8)<<"  , I(r12) = "<<circuit1.current(6, 9)<<"  ";
 //cout<<" Req =  "<<vdd/circuit1.totalCurrentSupply()<<endl;
 cout<<endl;
 }
 
 cout<<"temperature dependent resistance : "<<endl;
 
 
 TempCircuit circuit2;
 circuit2.setCircuit(9,1,8,configuration);
 circuit2.setParameter(1,1,1);
 //circuit2.setVDD(0.1);
 //circuit2.resetVoltage();
 for (double vdd = 0.0; vdd<=5.0; vdd+=0.0001) {
 circuit2.setVDD(vdd);
 //circuit2.resetVoltage();
 v = circuit2.eval(0.00000001, 100000);
 cout<<fixed << setw( 2 ) << setprecision( 4 ) <<" for vdd = "<<vdd<<" (Volt) , total current supplied = "
 <<fixed << setw( 7 ) << setprecision( 6 ) <<circuit2.totalCurrentSupply()<<"  , ";
 cout<<"T(r2) = "<<circuit2.temperature(4, 7)<<" , T(r7) = "<<circuit2.temperature(5, 8)<<"  , T(r12) = "<<circuit2.temperature(6, 9);
 //cout<<" Req =  "<<vdd/circuit2.totalCurrentSupply()<<endl;
 cout<<endl;
 }
 
*/



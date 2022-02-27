//
//  Circuit.cpp
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/3.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "Circuit.h"

Circuit::Circuit(){}
double Circuit::I(double V){return NAN;}
double Circuit::dIdV(double V){return NAN;}
void Circuit::setCircuit(int num,int vdd,int gnd,MAT &adjecentM){
    node_num = num;
    vdd_index = vdd;
    gnd_index = gnd;
    this->adjecentM = &adjecentM;
    this->voltages = new VEC(num);
}

void Circuit::resetVoltage(){
    cout<<"rst"<<endl;
    this->voltages = new VEC(node_num);
//    for (int i = 0; i<node_num; i++) {
//        (*this->voltages)[i] = rand()%10000*VDD/100000000.;
//    }
//    (*this->voltages)[0] = VDD/1000;
    VEC sequence(6);
    for (int i = 0; i<6; i++) {
        sequence[i] = (rand()%10000)*VDD/100000.;
    }
    sequence.sort();
    (*this->voltages)[0] = sequence[0];
    (*this->voltages)[1] = sequence[1];
    (*this->voltages)[2] = sequence[2];
    (*this->voltages)[5] = sequence[3];
    (*this->voltages)[8] = sequence[4];
    (*this->voltages)[7] = sequence[5];
    //initialize the restricted nodes
    (*this->voltages)[4] = rand()%10000*((*this->voltages)[5]-(*this->voltages)[1])/10000. + (*this->voltages)[1];
    (*this->voltages)[3] = rand()%10000*((*this->voltages)[4]-(*this->voltages)[0])/10000. + (*this->voltages)[0];
    (*this->voltages)[6] = rand()%10000*((*this->voltages)[7]-(*this->voltages)[3])/10000. + (*this->voltages)[3];
    //VEC &v = (*this->voltages);
    //if (v[0]>v[1]&&v[1]>v[2]&&v[2]>v[5]&&v[5]>v[8]&&v[8]>v[7]&&v[1]>v[4]&&v[4]>v[5]&&v[0]>v[3]&&v[3]>v[4]&&v[3]>v[6]&&v[6]>v[7]) {cout<<"ok";}
    //else{cout<<"bad"<<endl;}
}


double Circuit::F(int node,VEC &Vs){
    double result;
    if (node == vdd_index) {
        return Vs[vdd_index-1]-VDD;
    }
    if (node == gnd_index) {
        return Vs[gnd_index-1];
    }
    for (int i = 0; i<node_num; i++) {
        if ((*this->adjecentM)[node-1][i]!=0.) {
            result+=I(Vs[node-1]-Vs[i]);
        }
    }
    return result;
}
double Circuit::dFdV(int fi,int vj,VEC &vs){
    double result=0.;
    if (fi == vdd_index) {
        if (vj == vdd_index) {
            return 1.;
        }else{
            return 0.;
        }
    }
    if (fi == gnd_index) {
        if (vj == gnd_index) {
            return 1.;
        }else{
            return 0.;
        }
    }
    if (fi == vj) {
        for (int i = 0; i<node_num; i++) {
            if ((*this->adjecentM)[fi-1][i]==1.) {
                result+=dIdV(vs[fi-1]-vs[i]);
            }
        }
    }else{
        if ((*this->adjecentM)[fi-1][vj-1]==1.) {
            result-=dIdV(vs[fi-1]-vs[vj-1]);
        }
        
    }
    return result;
}
VEC Circuit::Fs(VEC &vs){
    VEC result(vs.len());
    for (int i = 0; i<vs.len(); i++) {
        result[i] = F(i+1, vs);
    }
    return result;
}

MAT Circuit::Jacobian(VEC &vs){
    MAT result(vs.len()); //result[row][col]
    
    for (int i = 0; i<vs.len(); i++) {
        for (int j = 0 ; j<vs.len(); j++) {
            result[i][j] = dFdV(i+1, j+1, vs);
        }
    }
    return result;
}
VEC Circuit::eval(double e,int maxIter){
    VEC &x = *voltages;
    for (int j = 0; j<node_num; j++) {
        x[j] = x[j] + rand()%100000*0.000001/100000.;
    }
    int k = 0;
    double err = 1+e;
    VEC F(x.len());
    MAT JF(x.len());
    VEC dx(x.len());
    while (err>e&&k<maxIter) {
        F=Fs(x);
        JF = Jacobian(x);
        luFact(JF);
        dx = bckSubs(JF, fwdSubs(JF, -F));
        x = x + dx;
        //if the voltages are out of bound, set to random number
//        for (int i = 0; i< node_num; i++) {
//            if (x[i]>VDD||x[i]<0) {
//                for (int j = 0; j<node_num; j++) {
//                    x[j] = rand()%100000*VDD/100000.;
//                }
//                i = node_num;
//            }
//        }
        if (x[7]>0&&VDD>x[0]&&x[0]>x[1]&&x[1]>x[2]&&x[2]>x[5]&&x[5]>x[8]&&x[8]>x[7]&&x[1]>x[4]&&x[4]>x[5]&&x[0]>x[3]&&x[3]>x[4]&&x[3]>x[6]&&x[6]>x[7]){
            
        }
        else{
            //resetVoltage();
            //x = *voltages;
        }
        
        k++;
        err = F.norm(0);
    }
    if (k>=maxIter-1) {
        cout<<"iteration count too large!"<<endl;
    }
    cout<<"k = "<<k;
    this->voltages = &x;
    return x;
}
double Circuit::totalCurrentSupply(){
    double current;
    for (int i = 0; i<node_num; i++) {
        if ((*this->adjecentM)[gnd_index-1][i]!=0.) {
            current+=I((*voltages)[i]);
        }
    }
    return current;
}
double Circuit::current(int node1,int node2){
    if (node1==node2||(*this->adjecentM)[node1-1][node2-1]==0.) {
        cout<<"no divice between node "<<node1<<" and node "<<node2 <<endl;
        return NAN;
    }
    return fabs(I((*voltages)[node1-1]-(*voltages)[node2-1]));
}


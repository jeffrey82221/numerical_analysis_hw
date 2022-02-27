//
//  100061212 林奕勳 HW6
//  Numerical_Analysis
//
//  Created by 林奕勳 on 2015/4/13.
//  Copyright (c) 2015年 林奕勳. All rights reserved.
//

#include <iostream>
#include <fstream>
//#include "Complex.h"
#include "VEC.h"
#include "MAT.h"
#include <math.h>
#include <string>
#include <sstream>
using namespace std;

// with this function, we can find if the nodes is on the corner of the circuit given the index of a node.
bool isVirtex(int node,int size){
    if (node == size*size) { // right-buttom corner
        return true;
    }else if(node == 1){ // left-top corner
        return true;
    }else if(node == size*size - size + 1){ //left-buttom corner
        return true;
    }else if(node == size){ // right-top corner
        return true;
    }else{
        return false;
    }
}
// with this function, we can find if the nodes is on the edge of the circuit given the index of a node.
bool isEdge(int node, int size){
    if(node > 1&& node <size){ //top edge
        return true;
        
    }else if(node > size*size - size + 1 && node < size*size){ //buttom edge
        return true;
    }else if(node%size == 0 && isVirtex(node, size) == false){ //right edge
        return true;
        
    }else if(node%size == 1 && isVirtex(node, size) == false){ //left edge
        return true;
        
    }else
        return false;
}
//with this function, we can find the corresponsed connected nodes' index of a indexed nodes and put them in a Vector.
VEC neighbor(int node,int size){
    VEC neighbor(4);
    if (node == 1) { //left-top corner
        neighbor[0] = node+1;
        neighbor[1] = node+size;
    }else if(node == size){ //right-top corner
        neighbor[0] = node-1;
        neighbor[1] = node+size;
    }else if(node == size*size - size + 1){ //left-buttom corner
        neighbor[0] = node+1;
        neighbor[1] = node-size;
    }else if(node == size*size){ // right-buttom corner
        neighbor[0] = node-1;
        neighbor[1] = node-size;
    }else if(node > 1&& node <size){ // top edge
        neighbor[0] = node-1;
        neighbor[1] = node+1;
        neighbor[2] = node+size;
    }else if(node > size*size - size + 1 && node < size*size){ //buttom edge
        neighbor[0] = node-1;
        neighbor[1] = node+1;
        neighbor[2] = node-size;
    }else if(node%size == 0 && isVirtex(node, size) == false){ //right edge
        neighbor[0] = node-size;
        neighbor[1] = node+size;
        neighbor[2] = node-1;
    }else if(node%size == 1 && isVirtex(node, size) == false){ //left edge
        neighbor[0] = node-size;
        neighbor[1] = node+size;
        neighbor[2] = node+1;
    }else{
        neighbor[0] = node-size;
        neighbor[1] = node+size;
        neighbor[2] = node-1;
        neighbor[3] = node+1;
    }
    return neighbor;
    
}
//with this function , we can find the index of the GND node corresponds to the size of the circuit
int findGound(int size){
    return (2*size*size - size + 1)/2;
}
// with this function, we can find the submatrix of a Matrix which elimate the iith colume and the jjth row
MAT findMinor(int ii,int jj,MAT in){
    MAT out(in.dim()-1);
    for (int i = 0; i<ii; i++) {
        for (int j = 0; j<jj; j++) {
            out[i][j] = in[i][j];
        }
    }
    for (int i = ii; i<out.dim(); i++) {
        for (int j = 0; j<jj; j++) {
            out[i][j] = in[i+1][j];
        }
    }
    for (int j = jj; j<out.dim(); j++) {
        for (int i = 0; i<ii; i++) {
            out[i][j] = in[i][j+1];
        }
    }
    for (int i = ii; i<out.dim(); i++) {
        for (int j = jj; j<out.dim(); j++) {
            out[i][j] = in[i+1][j+1];
        }
    }
    return out;
}

// with this function , we can find Linear System of A,b of Ax = b with Stamping Approach
MAT findSMatrixRN(int size,double R,VEC &b){
    MAT circuit(size*size);
    //set up the Matrix which represent the K'F equations with unknown [v1,v2,v3,...,vm]
    //where the nodes in the squared-shaped circuit is indexed in the order from left to right and then from top to buttom
    for (int i = 0; i<size*size ; i++) {
        if (isVirtex(i+1, size)) {
            circuit[i][i] = 2;
            circuit[i][neighbor(i+1, size)[0]-1] = -1;
            circuit[i][neighbor(i+1, size)[1]-1] = -1;
        }else if(isEdge(i+1, size)){
            circuit[i][i] = 3;
            circuit[i][neighbor(i+1, size)[0]-1] = -1;
            circuit[i][neighbor(i+1, size)[1]-1] = -1;
            circuit[i][neighbor(i+1, size)[2]-1] = -1;
        }else{
            circuit[i][i] = 4;
            circuit[i][neighbor(i+1, size)[0]-1] = -1;
            circuit[i][neighbor(i+1, size)[1]-1] = -1;
            circuit[i][neighbor(i+1, size)[2]-1] = -1;
            circuit[i][neighbor(i+1, size)[3]-1] = -1;
        }
    }
    
    //find the node represent ground
    int ground = findGound(size);
    //VEC b(size*size-1); //8
    b[0] = -2;    //-2V
    b[1] = 1;     //V
    b[size] = 1;  //V
    
    MAT A_=circuit.tpose();
    
    // the first colume of Matrix A associate with the unknown IR
    A_[0] = VEC(size*size);
    A_[0][0] = -1; //-IR
    A_[0][ground-1] = 1; //IR
    
    //combine the ground node K'f equation to the first equation, so that the matrix will be banded
    A_ = A_.tpose();
    A_[0] -= A_[ground-1];
    
    //after combine the ground node's K'f equation, we can erase the row and colume associate with the Ground nodes, because we already know that the voltage of GND is 0
    return findMinor(ground-1, ground-1, A_);
    
}

bool isFixVoltage(int size,int node){
    if (node == 1||node == findGound(size)) {
        return true;
    }else return false;
}
MAT findMatrixRN(int size,double R,VEC &b){
    
    MAT A(size*size);
    
    for (int i = 0; i<size*size ; i++) {
        if (isVirtex(i+1, size)) {
            A[i][i] = 2;
            A[i][neighbor(i+1, size)[0]-1] = -1;
            A[i][neighbor(i+1, size)[1]-1] = -1;
        }else if(isEdge(i+1, size)){
            A[i][i] = 3;
            A[i][neighbor(i+1, size)[0]-1] = -1;
            A[i][neighbor(i+1, size)[1]-1] = -1;
            A[i][neighbor(i+1, size)[2]-1] = -1;
        }else{
            A[i][i] = 4;
            A[i][neighbor(i+1, size)[0]-1] = -1;
            A[i][neighbor(i+1, size)[1]-1] = -1;
            A[i][neighbor(i+1, size)[2]-1] = -1;
            A[i][neighbor(i+1, size)[3]-1] = -1;
        }
    }
    for (int i = 0; i<size*size; i++) {
        for (int j = 0; j<size*size; j++) {
            if (isFixVoltage(size,i+1)) {
                b[i]=1;
                if (i==j) {
                    A[i][j]=1;
                }else{
                    A[i][j]=0;
                    
                }
            }else{
                if (isFixVoltage(size, j+1)) {
                    A[i][j]=0;
                }
            }
        }
    }
    //for VDD put all neighbor index b to 1
    for (int i = 0; i<2; i++) {
        b[neighbor(1, size)[i]-1]=1;
    }
    //for GND put all neighbor index b to 1
    for (int i = 0; i<3; i++) {
        b[neighbor(8, size)[i]-1]=0;
    }
    
    return A;
    
}
int main(int argc, const char * argv[]) {
    
    int resistor_length = 40;
    if (argc>1) {
        istringstream iss(argv[1]);
        iss>>resistor_length;
    }
    cout<<"resistor per side = "<<resistor_length<<endl;
    
    int size = resistor_length+1;
    double R = 2000/resistor_length;
    //VEC b(size*size-1);
    //MAT A = findMatrixRN(size, R, b);
    VEC b(size*size);
    MAT A = findMatrixRN(size,R,b);
    cout<<"formulation done!"<<endl;
    //cout<<"A = \n"<<A<<endl;
    A = A/R;
//    
//    MAT LU(A);
//    LU = luFact(LU);
//    VEC x_lu = bckSubs(LU, fwdSubs(LU, b));
//    cout<<"By LU Decomposition:"<<endl;
//    cout<<"Equivalent Resistance = "<<1/(2-x_lu[1]-x_lu[size])*R<<"ohm"<<endl;
//    cout<<"Vse = "<<x_lu[size*size-1]<<"Volt"<<endl;
//    cout<<"Vne = "<<x_lu[size-1]<<"Volt"<<endl;
//    cout<<"Vsw = "<<x_lu[size*size-size]<<"Volt"<<endl<<endl;
//    
//    double tor = 0.1;
//    for (int i = 0; i<4; i++) {
//        cout << "tor = "<<tor<<endl;
    
//        VEC x(b.len());
//        int iter = cg(A, b, x, 10000,0.0000001);
//        cout<<"By CG:"<<endl;
//        cout<<"Iteration count = "<<iter<<endl;
//        cout<<"Equivalent Resistance = "<<1/(2-x[1]-x[size])*R<<"ohm"<<endl;
//        cout<<"Vse = "<<x[size*size-1]<<"Volt"<<endl;
//        cout<<"Vne = "<<x[size-1]<<"Volt"<<endl;
//        cout<<"Vsw = "<<x[size*size-size]<<"Volt"<<endl;
    
//        cout<<"Comparison:"<<endl;
//        cout<<"dReq = "<<(1/(2-x[1]-x[size])*R)-(1/(2-x_lu[1]-x_lu[size])*R)<<"ohm"<<endl;
//        cout<<"dVse = "<<x[size*size-1]-x_lu[size*size-1]<<"Volt"<<endl;
//        cout<<"dVne = "<<x[size-1]-x_lu[size-1]<<"Volt"<<endl;
//        cout<<"dVsw = "<<x[size*size-size]-x_lu[size*size-size]<<"Volt"<<endl<<endl;
//        tor = tor*0.01;
//    }
    
    double Lmax=0.0,Lmin=0.0,shiftedLmax,shiftedLmin;
    int iterP,iterIP;
    
    double tol = 0.01;
//    for (double tol = 0.0000001; tol<=1; tol*=10) {

        cout<<"\ntol="<<tol<<endl;
        iterP =powerMethod(A, Lmax, *new VEC(A.dim()), 100000, tol);
        iterIP =invPowerMethod(A, Lmin, *new VEC(A.dim()), 100000, tol);
        cout<<"iteration of power method = "<<iterP<<endl;
        cout<<"iteration of inverse power method = "<<iterIP<<endl;
        cout<<"condition number = Lmax/Lmin = "<<Lmax<<"/"<<Lmin<<"="<<Lmax/Lmin<<endl;

//    }
    
//    for (double delta = 0.0000001; delta<=1; delta*=10) {

        cout<<"after using shifted power method to find Lmax: /n iteration = "<<iterP+sinvPowerMethod(A, shiftedLmax, *new VEC(A.dim()), 1000000, 0.000000001,Lmax)<<endl;
        cout<<"after using shifted power method to find Lmin: /n iteration = "<<iterIP+sinvPowerMethod(A, shiftedLmin, *new VEC(A.dim()), 1000000, 0.000000001,Lmin)<<endl;
//    }

        cout<<"condition number = Lmax/Lmin = "<<shiftedLmax<<"/"<<shiftedLmin<<"="<<shiftedLmax/shiftedLmin<<endl;
//    }
//    MAT B(3);
//    B[0][0]=3; B[0][1] = 2; B[0][2] = 1;
//    B[1][0]=2; B[1][1] = 3; B[1][2] = 2;
//    B[2][0]=1; B[2][1] = 2; B[2][2] = 3;
//    cout<<B<<endl;
//    MAT Q(3);
//    MAT R_(3);
//    cout<<Q<<endl;
//    MAT T = B;
//
//    for (int  i = 0; i<10; i++) {
//        QRDecomposition(T, Q , R_);
//        T = Q*R_;
//    }
//    cout<<"Q = \n"<<Q<<endl;
//    cout<<"T = \n"<<T<<endl;
}

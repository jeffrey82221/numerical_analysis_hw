//
//  main.cpp
//  Numerical_Analysis
//
//  Created by jeffrey_MAC on 2015/3/3.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include <iostream>
#include "Complex.h"
#include "VEC.h"
#include "MAT.h"
#include <math.h>
using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    // Gram-Schmit Process
    double m3[9] = {3,2,0,2,3,2,1,2,3};
    MAT A(3,m3);
    cout<<"A = "<<endl;
    cout<<A<<endl;
    A = A.tpose();
    MAT G(3);
    
    //first algo:
    /*
    G[0] = A[0]; // first step
    for (int k = 1; k<A.dim(); k++) {
        for (int i = 0; i<k; i++) {
            G[k] += ((A[k]*G[i])*G[i])/(G[i]*G[i]);
        }
        G[k] = A[k] - G[k];
    }
    
    //second algo
    G[0] = A[0]; // first step
    for (int k = 1; k<A.dim(); k++) {
        G[k] = A[k];
        for (int i = 0; i<k; i++) {
            G[k] -= ((A[k]*G[i])*G[i])/(G[i]*G[i]);
        }
    }
    */
    //third algo
    G[0] = A[0]; // first step
    for (int k = 1; k<A.dim(); k++) {
        G[k] = A[k];
        for (int i = 0; i<k; i++) {
            G[k] -= (A[k]*G[i])/(G[i]*G[i])*G[i];
        }
    }

    G = G.tpose();
    cout<<"G = "<<endl;
    cout<<G<<endl;
    cout<<"D = "<<endl;
    MAT D = G.tpose()*G;
    cout<<D<<endl;
    
    double sc = 0;
    
    for (int i = 0; i<3; i++) {
        for (int j = 0 ; j<3; j++) {
            if (i == j) {
                continue;
            }else{
                sc += (D[i][j])*(D[i][j]);
            }
        }
    }
    cout<<"sc = \n"<<sc<<endl;

}

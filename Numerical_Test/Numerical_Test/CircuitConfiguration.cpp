//
//  CircuitConfiguration.cpp
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/3.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#include "CircuitConfiguration.h"


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

//out of this class
int ResistorIndex(int node1,int node2,int size){
    int small = node1<node2?node1:node2;
    int big = node1>node2?node1:node2;
    
    if (big-small == 1) { //parallel resistor
        return ((small-1)%size)*(size+size-1)+(size)+small/3;
    }else{ //verticle resistor
        return ((small-1)%size)*(size+size-1)+small/3;
    }
}
int smallNode(int index,int size){
    int fold = (index-1)%(size+size-1)+1;
    if (fold<size) {
        return (fold-1)*size+1+index/(size+size-1);
    }else{
        return (fold-size)*size+1+(index-1)/(size+size-1);
    }
    
}
int bigNode(int index,int size){
    int fold = (index-1)%(size+size-1)+1;
    if (fold<size) {
        return (fold)*size+1+(index)/(size+size-1);
    }else{
        return (fold-size)*size+1+(index-1)/(size+size-1)+1;
    }
}
MAT adjMatrix(int size){
    MAT adjM(size*size);
    VEC tmp(4);
    for (int i = 0; i<size*size; i++) {
        tmp = neighbor(i+1, size);
        for (int j = 0; j<4; j++) {
            if (tmp[j]==0) {
                break;
            }
            adjM[i][tmp[j]-1] = 1.;
        }
    }
    return adjM;
}

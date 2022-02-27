//
//  CircuitConfiguration.h
//  Numerical_Test
//
//  Created by jeffrey_MAC on 2015/6/3.
//  Copyright (c) 2015年 jeffrey_MAC. All rights reserved.
//

#ifndef __Numerical_Test__CircuitConfiguration__
#define __Numerical_Test__CircuitConfiguration__

#include <stdio.h>
#include "MAT.h"
#include "VEC.h"
bool isVirtex(int node,int size); //given a square network size, find if the node is an virtex
bool isEdge(int node, int size); //given a square network size, find if the node is along the edge
VEC neighbor(int node,int size); //given a .............. size, find if node's connected neighbors and put the index in a VEC(4)
int findGound(int size);    //given a square network size, find the gnd index
int ResistorIndex(int node1,int node2,int size); //given two node index and the square network size, find the index of divice between
int smallNode(int index,int size); //given the network size and the divice's index, find the smaller index of the connected node.
int bigNode(int index,int size); //given the network size and the divice's index, find the bigger index of the connected node.
MAT adjMatrix(int size); //given the network size, find the adjecent matrix with which the element i,j equal's 1 if node i and node j are connected.

#endif /* defined(__Numerical_Test__CircuitConfiguration__) */

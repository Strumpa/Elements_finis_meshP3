# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 17:50:44 2023

@author: Loutre
"""
import numpy as np 

def MESHP3(XA,XB,NELM):
    NNODE =4
    NNOD = 3*NELM +1
    COORD = np.zeros(NNOD)
    for i in range(NNOD):
        COORD[i]=XA+(XB-XA)/NELM/3*(i) 
    CONNECD = np.zeros([NELM,NNODE])
    CONNECBA= np.zeros([1,1])
    CONNECBB= np.zeros([1,1])
    for k in range(NELM):
        CONNECD[k][0] = 3*k+1
        CONNECD[k][1] = 3*k+3
        CONNECD[k][2] = 3*k+4
        CONNECD[k][3] = 3*k+2
    CONNECBA[0][0] = 1
    CONNECBB[0][0] = NNOD
	
    return (COORD, NNOD, CONNECD, CONNECBA, CONNECBB)

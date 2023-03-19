# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 18:05:31 2023

@author: Loutre
"""
import numpy as np
from GLQ1D import *

#Diffusion et terme source avec des polynomes P2 

def DIFF_P3(NUMER,CONNEC,COORD,T,DATA_DIFFUSION,VAR):
    DIMI = 4
    DIMJ = 4
    NNODE= 3
    #On recupere les foctions Kd (diffusivité) et F(source volumique)
    fctKD = DATA_DIFFUSION[0]
    fctF  = DATA_DIFFUSION[1]
    NPQ = 4
    # Choix de la regle d'integration
    [DEG, W, Xi] = GLQ1D(NPQ)
    # Gestion du syteme elementaire
    AE = np.zeros([DIMI,DIMJ])
    RE = np.zeros(DIMI)
    XE = np.zeros(NNODE)
    TE = np.zeros(DIMJ)
    # On recupere le vecteur d'adressage de la variable primaire
    VAD = NUMER.VAD(CONNEC,VAR[0])
    if (len(VAD) != DIMJ):
        print('Verifier la longueur de adr (%d) \n', len(VAD))
        #On recupere les coordonnees 
    for i in range(NNODE):
        XE[i] = COORD[int(CONNEC[i])-1]
        #On recupere la variable primaire
    for j in range(DIMJ):
        TE[j] = T[int(VAD[j])-1]
    for j in range(NPQ):
        r=Xi[j]
        #Interpolants de la geometrie
        LG1  = (1 - r) / 2.
        LG2  = (1 + r) / 2.
        DLG1 = -1./2.
        DLG2 = 1./2.
        Xh   = LG1*XE[0] + LG2*XE[1] #Jacobien de la transformation
        J    = DLG1*XE[0] + DLG2*XE[1]
        #Evaluation de la diffusivite et du terme source 
        KD   = fctKD(Xh)
        f    = fctF(Xh)
        #Interpolants de la variiable primaire
        LT1  = (r+(1/3))*(r-(1/3))*(r-1)/(-16/9)
        LT2  = (r+1)*(r+(1/3))*(r-(1/3))/(16/9)
        LT3  = (r+1)*(r-(1/3))*(r-1)/(16/27)
        LT4  = (r+1)*(r+(1/3))*(r-1)/(-16/27)
        DLT1 = (1/16)*(-27*r**2+18*r+1)
        DLT2 = (1/16)*(27*r**2+18*r-1)
        DLT3 = (9/16)*(9*r**2-2*r-3)
        DLT4 = (-9/16)*(9*r**2+2*r-3)
        DLT1dx = DLT1/J
        DLT2dx = DLT2/J
        DLT3dx = DLT3/J
        DLT4dx = DLT4/J

        LT   = [LT1,LT2,LT3,LT4]
        dLTdx= [DLT1dx,DLT2dx,DLT3dx,DLT4dx]
        dTdx = np.dot(dLTdx,TE)

        #Residu elementaire
        for k in range(DIMI):
            RE[k]=RE[k]-(f*LT[k]*J*W[j])
            RE[k]=RE[k]+(KD*dTdx*dLTdx[k]*J*W[j])
        for k in range(DIMI):
            for l in range(DIMJ):
                AE[k][l]+=dLTdx[k]*KD*dLTdx[l] * J * W[j]
    VADI = VAD
    VADJ = VAD
    return(AE,RE,VADI,VADJ)

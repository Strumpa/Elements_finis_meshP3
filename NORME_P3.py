# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 15:00:24 2023

@author: Loutre
"""
from GLQ1D import *
import numpy as np

def NORME_P3(NUMER,CONNECD, COORD, T, VAR, FCTSOLEXACT):
        NNODE = 3
        DIMJ  = 4
        normeL2 = 0
        normeH1 = 0 
        seminormeH1 = 0
        NELM = CONNECD.shape[0]
        fctF = FCTSOLEXACT[0] #FCTSOLEXACT.F
        fctdFdx = FCTSOLEXACT[1] #FCTSOLEXACT.dFdx
        #Choix de la regle d'integration
        NPQ = 4
        [DEG,W,Xi] = GLQ1D(NPQ)
        #Gestion du systeme elementaire 
        XE = np.zeros(NNODE)
        TE = np.zeros(DIMJ)
        for  k in range(int(NELM)):
                CONNEC = CONNECD[k,:]
                VAD = NUMER.VAD(CONNEC,VAR)
                for i in range(NNODE):
                        XE[i] = COORD[int(CONNEC[i])-1]
                for i in range(DIMJ):
                        TE[i]=T[int(VAD[i])-1]
                for j in range(NPQ):
                        r=Xi[j]
                        LG1 = (1-r)/2
                        LG2 = (1+r)/2
                        DLG1 = -1/2
                        DLG2 = 1/2
                        Xh = LG1*XE[0]+LG2*XE[1]
                        J = DLG1*XE[0]+DLG2*XE[1]
                        #Evaluation de la diffusivite et du terme source
                        TExact = fctF(Xh)
                        dTExact_dx = fctdFdx(Xh)
                        # Interpolants de la variable primaire
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


                        LT =np.array([LT1, LT2, LT3,LT4])
                        dLTdx = np.array([DLT1dx, DLT2dx, DLT3dx, DLT4dx])
                        Th = np.dot(LT,TE)
                        dTdx = np.dot(dLTdx , TE)
                        # Norme L2
                        normeL2 += (TExact-Th)**2*J*W[j]
                        seminormeH1 += (dTExact_dx-dTdx)**2*J*W[j]
        normeH1= np.sqrt(normeL2+seminormeH1)
        seminormeH1 = np.sqrt(seminormeH1)
        normeL2 = np.sqrt(normeL2)
        return (normeL2,normeH1,seminormeH1)


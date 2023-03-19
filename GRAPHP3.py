# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 15:26:14 2023

@author: Loutre
"""
import numpy as np
import matplotlib.pyplot as plt

def GRAPHP3 (NUMER, CONNECD, COORD, T, VAR):
        fig2=plt.figure()
        ax2 = fig2.add_subplot(1, 1, 1)
        NNODE = 3
        DIMJ  = 4
        NELM  = CONNECD.shape[0]
        XE    = np.zeros(NNODE)
        TE    = np.zeros(DIMJ)
        for k in range(NELM):
                CONNEC = CONNECD[k][:]
                VAD= NUMER.VAD(CONNEC,VAR)
                for i in range(NNODE):
                        XE[i]= COORD[int(CONNEC[i])-1]
                for i in range(DIMJ):
                        TE[i]=T[int(VAD[i])-1]
                Npt = 20
                pt  = np.linspace(-1,1,Npt)
                X   = np.zeros(Npt)
                TX  = np.zeros(Npt)
                for i in range(Npt):
                        r = pt[i]
                        #Interpolants de la géométrie 
                        LG1 = (1-r)/2.
                        LG2 = (1+r)/2.
                        #Interpolants de la variable primaire
                        LT1  = (r+(1/3))*(r-(1/3))*(r-1)/(-16/9)
                        LT2  = (r+1)*(r+(1/3))*(r-(1/3))/(16/9)
                        LT3  = (r+1)*(r-(1/3))*(r-1)/(16/27)
                        LT4  = (r+1)*(r+(1/3))*(r-1)/(-16/27)
                        X[i] = LG1*XE[0] + LG2*XE[1]
                        TX[i] = LT1*TE[0] + LT2*TE[1] + LT3*TE[2] + LT4*TE[3]
                ax2.plot(X,TX,'k')
                ax2.set_title(VAR)
        plt.show()
        return  

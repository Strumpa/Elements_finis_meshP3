# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 17:58:07 2023

@author: Loutre
"""
#Importation des fonctions utilisées pour la méthode éléments finis
from NUMEROTATION1D import *

#Routines de maillage
from MESHP1 import *
from MESHP2 import *
#Routines des équations sur le domaine
from MASSE import *
from DIFF_P1 import *
from DIFF_P2 import *
from DIFF_P2H import *
from DIFFcp_P1 import *
from DIFFcp_P2 import *
from DIFFcp_P2H import *
from CONVNC_P1 import * 
from CONVNC_P2 import * 
from CONVNC_P2H import * 
from CONVCS_P1 import * 
from CONVCS_P2 import * 
from CONVCS_P2H import * 

#Routines des équations sur un élément (condition de neumann, robin, flux radiatif ...)
from SYSFLUXCONV import *
from FLUXDIFF import *
from FLUXRADIA import *

from CORRECTIONRESIDU import *
from RESOUDRE import *

#Routines de post-traitement
from REACTION1D import * 
from GRAPHP1 import * 
from GRAPHP2 import *
from GRAPHP2H import *
from NORME_P1 import *
from NORME_P2 import *
from NORME_P2H import * 

#Routines d'affichage
from affichage_maillage import *
from affichage_numerotation import *
from affichage_ART import *

#Bibliothèques utiles de python
import matplotlib.pyplot as plt
import numpy as np

#importation des fonctions custom
from MESHP3 import * 
from DIFF_P3 import *
from CONVNC_P3 import *
from NORME_P3 import *
from GRAPHP3 import *

NELM = 5
XA = 0.0
XB = 1.0
print("On considère maintenant un maillage P3")
(COORD, NNOD, CONNECD1, CONNECBA, CONNECBB) = MESHP3(XA,XB, NELM)
kD=70.0
u=1.0
#Affichage du maillage
affichage_maillage(COORD,NNOD)
print("Coord : ", COORD)
print("N Nodes : ", NNOD)
print("ConnecBA : ", CONNECBA)
print("ConnecBB : ", CONNECBB)
print("ConnecD1 : ", CONNECD1)
#Proprietes pour les equations et membres de droite
####################################################
#           Diffusion: -d(KD*dT/dx)/dx=F           #
####################################################

DATA_DIFFUSION = [lambda x:kD, lambda x:0] #KD et F resp. 

####################################################
#           Convection: U*dT/dx=0                  #
####################################################
DATA_CONVECTION = [lambda x:u] #U

####################################################
#           Masse: M*T=0                           #
####################################################
DATA_MASSE = [lambda x:0] #M

####################################################
#Flux diffusif: k*dT/dn=F+B*T (condition de robin) #
####################################################
DATA_FLUXDIFF = [lambda x:0, lambda x:0] #F et B resp.

####################################################
#      Flux radiatif:  k*dT/dn=S*(T0**4 - T**4)      #
####################################################
DATA_FLUXRADIA = [lambda x:0, lambda x:0] #T0 et S resp.



# EQUATIONS POUR LES DIFFERENTES FORMES FAIBLES

EQS = [[DIFF_P3, DATA_DIFFUSION, CONNECD1, ['U'], [[1,1,1,1]]],
       [CONVNC_P3, DATA_CONVECTION, CONNECD1, ['U'], [[1,1,1,1]]],
       ]

# CONDITIONS LIMITES ESSENTIELLES                  
CLE = [[['U'], [1],0, CONNECBA],
[['U'], [1], 1, CONNECBB],
]

#Numerotation                                      
NUMEROTATION = NUMEROTATION1D(NNOD, ['U'],EQS, CLE)

#Affichage de la numérotation
affichage_numerotation(NUMEROTATION)


#Le systeme d'equations
NINC = NUMEROTATION.NINC
NDDL = NUMEROTATION.NDDL
A = np.zeros([NINC,NINC])
R = np.zeros(int(NDDL))
T = np.zeros(int(NDDL))


#Copie des conditions limites connues dans T(vecteur solution)
T=NUMEROTATION.SOLUTION(T,T)
T

#Resolution en residu-correction
tol = 1.0e-9
maxiter = 5
[T,R,A]=CORRECTIONRESIDU(tol,maxiter, RESOUDRE,EQS, NUMEROTATION, COORD, T, A, R)
#affichage_ART(A, R, T)


#Solutio exacte
FCTSOLEXACTZERO= [lambda x: (np.exp(u*x/kD)-1)/(np.exp(u/kD)-1), 
                  lambda x: (u*np.exp(u*x/kD))/(kD*(np.exp(u/kD)-1))]
                    #SOLUTION et DERIVEE resp.
#Affichage

GRAPHP3(NUMEROTATION,CONNECD1,COORD,T,'U')
L2,H1,semiH1 = NORME_P3(NUMEROTATION, CONNECD1, COORD, T, 'U', FCTSOLEXACTZERO)
print("Diff en L2:", L2)
print("Diff en H1 :", H1)
print("Diff en semiH1 :", semiH1)
print(REACTION1D(NUMEROTATION, CLE[0], R))
print(REACTION1D(NUMEROTATION, CLE[1], R))
erreur_0 = -kD*FCTSOLEXACTZERO[1](XA)-REACTION1D(NUMEROTATION, CLE[0], R)
erreur_1 = kD*FCTSOLEXACTZERO[1](XB)-REACTION1D(NUMEROTATION, CLE[1], R)

print("Diff entre flux et reaction à x=0 :", erreur_0)
print("Diff entre flux et reaction à x=1 :", erreur_1)
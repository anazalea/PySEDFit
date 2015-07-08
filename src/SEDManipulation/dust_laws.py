# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 21:43:08 2015

@author: anneya

                      .-. 
                 .--.(   ).--.                
      <-.  .-.-.(.->          )_  .--.
       `-`(     )-'             `)    )
         (o  o  )                `)`-'
        (      )                ,)   
        ( ()  )                 )
         `---"\    ,    ,    ,/`  
               `--' `--' `--'
                |  |   |   | 
                |  |   |   |
                '  |   '   |  
                
Functions to apply dust reddening curves to your sad blue :( spectra
"""

from __future__ import print_function
from __future__ import division

import sys
sys.path.append('.')
import math
import numpy as np

#------------------------------------------------------------------------------ 
def Calzetti2000(lam,lnu,ebv):
    '''
    Apply the Calzetti 20000 las as defined in Calzetti et al 2000 (ApJ 533 682) Equations 2,3,4
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''    
    r = 4.
    k = 2.659 * (-1.857 + 1.040/lnu) + r
    k[lam<0.63] = 2.659 * (-2.156 + (1.509/lnu[lam<0.63]) - (0.198/lnu[lam<0.63]**2) +\
                (0.011/lnu[lam<0.63]**3)) + r
    newLnu = lnu * 10 ** (-0.4 * ebv * k)
    return(newLnu)
    
#------------------------------------------------------------------------------ 
def Calzetti1997(lam,lnu,ebv):
    '''
    Apply the Calzetti 1997 las as defined in Calzetti et al 1997 (astro-ph/9706121) Eqns 1 and 2. 
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''    
    k =  (((1.86 - 0.48/lam)/lam - 0.1) / lam) + 1.73
    k[lam<0.63] = 2.6596 * (-2.156 + (1.509/lnu[lam<0.63]) - (0.198/lnu[lam<0.63]**2) +\
                (0.011/lnu[lam<0.63]**3)) + 4.88
    newLnu = lnu * 10 ** (-0.4 * ebv * k)
    return(newLnu)
 
#------------------------------------------------------------------------------ 
def Fitzpatrick(lamInv,c1,c2,c3,c4,lambda0inv,gamma,ebv):
    '''
    Inputs:
            lamInv = NumPy array of wavelengths in microns^-1
            c1,c2,c3,c4 = Fitzpatricky constants (floats)
            lambda0inv = some other Fitzpatricky constant
            gamma = ditto that ^
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of taus
    '''   
    c4v = np.zeros(len(laminv))
    c4v[laminv>=5.9] = c4
    
    elvOverEbv = c1 + (c2*lamInv) + c3/((lamInv - lambda0inv**2/lamInv)**2 + gamma**2) +\
                c4v * (0.593*(lamInv-5.9)**2 + 0.0564*(lamInv-5.9)**3)
    avOverEbv = 3.1 # Assumes R(V)=A(V)/E(B-V)=3.1
    tau = (elvOverEbv + avOverEbv) * ebv/1.086
    return(tau)
#------------------------------------------------------------------------------
def LMC(lam,lnu,ebv):
    '''
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)    
    Output:
            NumPy array of dust attenuated lnus
    '''   
    c1,c2,c3,c4 = -0.69,0.89,2.55,0.50
    lambda0inv = 4.608
    gamma = 0.994
    return(lam*Fitzpatrick(1./lam,c1,c2,c3,c4,lambda0inv,gamma,ebv))
#------------------------------------------------------------------------------  
def MW(lam,lnu,ebv):
    '''
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''   
    c1,c2,c3,c4 = -0.38,0.74,3.96,0.26
    lambda0inv = 4.595
    gamma = 1.051
    return(lam*Fitzpatrick(1./lam,c1,c2,c3,c4,lambda0inv,gamma,ebv))
#------------------------------------------------------------------------------  
def Dor30(lam,lnu,ebv):
    '''
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
    Output:
            NumPy array of dust attenuated lnus
    '''   
    c1,c2,c3,c4 = -2.19,1.39,1.49,0.43
    lambda0inv = 4.606
    gamma = 0.894
    return(lam*Fitzpatrick(1./lam,c1,c2,c3,c4,lambda0inv,gamma,ebv))
    



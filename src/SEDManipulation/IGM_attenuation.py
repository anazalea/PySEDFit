# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:27:26 2015

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


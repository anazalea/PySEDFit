# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:13:01 2015

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
                
Model selection tools common to all selection algorithms
"""

from __future__ import print_function
from __future__ import division

import sys
import math
import numpy as np

#------------------------------------------------------------------------------   
def Scale(modelFlux,dataFlux,dDataFlux):
    '''
    Calculate scale factor of best fit according to Sawicki (2012) A5 
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of mselected odels, array of size (#objects,#filters)
            
    Output:
            Scale factors in numpy array of size (#objects)
            
    '''
    #  If there's only one object or one model, we need to reshape so 
    #  that array operations are working on the right indices
    if len(np.shape(modelFlux))==1:   
        modelFlux = np.reshape(modelFlux,(1,len(modelFlux)))
    if len(np.shape(dataFlux))==1:
        dataFlux = np.reshape(dataFlux,(1,len(dataFlux)))
    if len(np.shape(dDataFlux))==1:
        dDataFlux = np.reshape(dDataFlux,(1,len(dDataFlux)))
    s = np.sum(dataFlux*modelFlux/dDataFlux**2,axis=1) / np.sum(modelFlux**2/dDataFlux**2,axis=1)
    return(s)
#------------------------------------------------------------------------------    



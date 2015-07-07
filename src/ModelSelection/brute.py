# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:15:52 2015

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
                
Brute force chi^2 model selection tools
These functions do a brute force chi^2 minimization to find the best fit model
"""

from __future__ import print_function
from __future__ import division

import sys
sys.path.append('.')
import math
import numpy as np
import fit_tools 

#------------------------------------------------------------------------------
def MinChi2(models,data,weights):
    '''
    Returns minimum chi2 value and index of best model
    
    Inputs:
            models = observables of models, array of size (#models,#filters)
            data = observed quantities for an object, array of size (#filters)
            weights = uncertainties in data, array of size (#filters)
            
    Output:
            NumPy array of size (#objects,3)
            Columns [index of model with min chi^2, scale factor, chi^2 value]
    '''
    chi2 = np.sum(((models-data)/weights)**2,axis=1)
    minChi2 = np.min(chi2)
    n = np.argmin(chi2)
    return(np.array([minChi2,n]))
    
#------------------------------------------------------------------------------   
def BruteColorSpace(dataFluxes,dDataFluxes,modelFluxes):
    '''
    Compare objects' observed colors and propagated uncertainties in N filters to model 
    spectra (sampled in the same N filters) to identify the model that produces 
    the lowest value of chi2 for each observed object
    
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            
    Output:
            NumPy array of size (#objects,3)
            Columns [index of model with min chi^2, scale factor, chi^2 value]
    '''
    modelColors = modelFluxes[:,1:] / modelFluxes[:,:-1]
    dataColors = dataFluxes[:,1:]/dataFluxes[:,:-1]
    dDataColors = np.sqrt( (1./dataFluxes[:,:-1])**2 * (dDataFluxes[:,1:])**2 \
                + (dataFluxes[:,1:]/dataFluxes[:,:-1]**2)**2 * (dDataFluxes[:,:-1])**2)
    results = np.array([]).reshape(0,3)
    for i in range(len(dataFluxes)):
        minChi2,n = MinChi2(modelColors,dataColors[i],dDataColors[i])
        s = fit_tools.Scale(modelFluxes[int(n)],dataFluxes[i],dDataFluxes[i])
        results = np.r_[results,[[n,s,minChi2]]]
    return(results)
    
#------------------------------------------------------------------------------     
def BruteFluxSpace(dataFluxes,dDataFluxes,modelFluxes):
    '''
    Compare objects' observed fluxes and associated uncertainties in N filters to model 
    spectra (sampled in the same N filters) to identify the model that produces 
    the lowest value of chi2 for each observed object
    
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            
    Output:
            NumPy array of size (#objects,3)
            Columns [index of model with min chi^2, scale factor, chi^2 value]
    '''
    results = np.array([]).reshape(0,3)
    for i in range(len(dataFluxes)):
        scales = fit_tools.Scale(modelFluxes,dataFluxes[i],dDataFluxes[i])
        scaledModelFluxes = (modelFluxes.transpose() * scales.transpose()).transpose()
        minChi2,n = MinChi2(scaledModelFluxes,dataFluxes[i],dDataFluxes[i])
        results = np.r_[results,[[n,scales[int(n)],minChi2]]]
    return(results)
#------------------------------------------------------------------------------ 

    
    
    




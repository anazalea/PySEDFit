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
import os

#------------------------------------------------------------------------------
def MinChi2(models,data,weights,returnAll=False):
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
    if returnAll:
        return([minChi2,n,chi2])
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
    
#-------------------------------------------------------------------------------     
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
def BruteFitErrorBars(dataFluxes,dDataFluxes,modelFluxes,modelParams,dChi2):
    '''
    OldSchool SEDfit errorbar fitting by chi2 difference
    The best fit (min chi2) model is found as in BruteFluxSpace, and then the subset of
    models with chi2 values within dChi2 of its minimum value are considered to give 
    ranges in each parameter
    
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            modelParams = parameters of each model to be recorded, array of size (#models,#parameters)
            dChi2 = Difference from min chi2 value for errorbar definition
            
    Output:
            NumPy array of size (#objects,3)
            Columns [index of model with min chi^2, scale factor, chi^2 value]
    '''
    results = np.array([]).reshape(0,3+2*len(modelParams[0]))
    for i in range(len(dataFluxes)):
        scales = fit_tools.Scale(modelFluxes,dataFluxes[i],dDataFluxes[i])
        scaledModelFluxes = (modelFluxes.transpose() * scales.transpose()).transpose()
        minChi2,n,chi2 = MinChi2(scaledModelFluxes,dataFluxes[i],dDataFluxes[i],returnAll=True)
        maxChi2 = minChi2 + dChi2
        result = [n,scales[int(n)],minChi2]
        for j in range(len(modelParams[0])): #for each parameter of interest
            allVals = modelParams[:,j][chi2<maxChi2]
            result.append(np.min(allVals))
            result.append(np.max(allVals))
        results = np.r_[results,[result]]
    return(results)
#------------------------------------------------------------------------------ 

def BruteFit(params,dataFlux,dDataFlux,dataParams,dataParamInfo,modelFlux,modelParams,modelParamInfo,dchi2=None):
    # params to add BruteSpace: color/flux
    #               BruteDChi2Error: dchi2 val         
    
    if dchi2==None:
        if params['brute_space']=='color':
            fit = BruteColorSpace(dataFlux,dDataFlux,modelFlux)
        else:
            fit = BruteFluxSpace(dataFlux,dDataFlux,modelFlux)
        outModelParams = modelParams[fit[:,0].astype(int)]
        modelParamNames = modelParamInfo[:,1]
        modelParamFluxScale = modelParamInfo[:,3]

        
        # this is ugly, but my brain isn't indexing well
        scaleUs = np.arange(len(modelParamNames))[modelParamFluxScale>0]
        
        scaleUs = [scaleUs]
        print(scaleUs)
        for n in scaleUs:
            outModelParams[:,n]*=fit[:,1]
            
        outStuff = np.c_[dataParams,outModelParams,fit[:,2]]
        head = ''
        fmt = ''
        i=0
        for p in dataParamInfo:
            head+=str(i)+' '+p[1]+'\n'
            fmt+=p[2]+' '
            i+=1
        for p in modelParamInfo:
            head+=str(i)+' '+p[1]+'\n'
            fmt+=p[2]+' '
            i+=1
        head+=str(i+1)+' min(chi2)'
        fmt+='%f'
        np.savetxt(params['output_file'],outStuff,header=head,fmt=fmt)
        os.system('say "I work! I work!"')
        return()
        
        
        
    elif params[brute_space]=='flux' and params[brute_dchi2]!=None:
        fit = BruteColorSpace(dataFlux,dDataFlux,modelFlux)
        
        
    
    
    




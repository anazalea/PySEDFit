# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 19:48:58 2015

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
                
Model selection with MC iterations
"""
from __future__ import print_function
from __future__ import division

import sys
sys.path.append('.')
import math
import numpy as np
import fit_tools 
from scipy.spatial.distance import seuclidean
from sklearn.neighbors import NearestNeighbors as nn
from sklearn.neighbors import BallTree
from astroML.density_estimation import bayesian_blocks
from astroML.plotting import hist as amlHist
import matplotlib.pyplot as plt

'''
Dual Tree, MC

Use MC method instead of creating a new model tree for each object. Create X 
sets of perturbed fluxes for each object. Use a dual tree query to find closest 
model to each (euclidean dist.) Return sampled PDF of each model parameter 
(and covariance matrix?)
'''
#------------------------------------------------------------------------------
def DualTreeFlux(dataFlux,dDataFlux,modelFlux,modelParams,mcIts,columnsToScale=[]):
    '''
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            modelParams = parameters of each model to be recorded, array of size (#models,#parameters)
            mcIts = number of times to perturb fluxes for each object, int
            columnsToScale = list of column indices in modelParams of parameters that need to be multiplied by scale factor
            
    Output:
            NumPy array of size (#objects,mcIts,#params)
            e.g. the zeroth element gives you a 2d array where each row represents the
            fit parameters from one monte carlo iteration 
    '''
    
    fitParams = []
    for i in range(len(dataFlux)):
        newFlux = dataFlux[i] + dDataFlux[i] * np.random.randn(mcIts,len(dataFlux[i]))
        scales = fit_tools.Scale(modelFlux,dataFlux[i],dDataFlux[i])
        scaledModelFlux = (modelFlux.transpose() * scales.transpose()).transpose()
        tree = BallTree(scaledModelFlux)
        query = tree.query(newFlux,k=1,dualtree=True)
        s = scales[query[1][:,0]]
        myParams = s
        for j in range(len(modelParams[0])):
            if j in columnsToScale:
                myParams = np.c_[myParams,np.multiply(s,modelParams[query[1][:,0]][:,j])] 
            else:
                myParams = np.c_[myParams,modelParams[query[1][:,0]][:,j]]
        fitParams.append(myParams)
    return(np.array(fitParams))
                

#------------------------------------------------------------------------------ 
def DualTree(dataFlux,dDataFlux,modelFlux,modelParams,mcIts,columnsToScale=[]):
    '''
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            modelParams = parameters of each model to be recorded, array of size (#models,#parameters)
            mcIts = number of times to perturb fluxes for each object, int
            columnsToScale = list of column indices in modelParams of parameters that need to be multiplied by scale factor
            
    Output:
            NumPy array of size (#objects,mcIts,#params)
            e.g. the zeroth element gives you a 2d array where each row represents the
            fit parameters from one monte carlo iteration 
    '''
    modelColors = modelFlux[:,1:] / modelFlux[:,:-1]
    tree = BallTree(modelColors)
    fitParams = []
    for i in range(len(dataFlux)):
        newFlux = dataFlux[i] + dDataFlux[i] * np.random.randn(mcIts,len(dataFlux[i]))
        newColors = newFlux[:,1:] / newFlux[:,:-1]
        query = tree.query(newColors,k=1,dualtree=True)
        s = fit_tools.Scale(modelFlux[query[1][:,0]],newFlux,np.ones(np.shape(newFlux)))
        myParams = s
        for j in range(len(modelParams[0])):
            if j in columnsToScale:
                myParams = np.c_[myParams,np.multiply(s,modelParams[query[1][:,0]][:,j])]                
            else:
                myParams = np.c_[myParams,modelParams[query[1][:,0]][:,j]]
        fitParams.append(myParams)
    return(np.array(fitParams))
#------------------------------------------------------------------------------ 
def DaisyChainDualTree(dataFlux,dDataFlux,modelFlux,modelParams,mcIts,columnsToScale=[]):
    '''
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            modelParams = parameters of each model to be recorded, array of size (#models,#parameters)
            mcIts = number of times to perturb fluxes for each object, int
            columnsToScale = list of column indices in modelParams of parameters that need to be multiplied by scale factor
            
    Output:
            NumPy array of size (#objects,mcIts,#params)
            e.g. the zeroth element gives you a 2d array where each row represents the
            fit parameters from one monte carlo iteration 
    '''
    modelColors = modelFlux[:,1:] / modelFlux[:,:-1]
    modelColors = np.c_[modelColors,modelFlux[:,0]/modelFlux[:,-1]]
    tree = BallTree(modelColors)
    fitParams = []
    for i in range(len(dataFlux)):
        newFlux = dataFlux[i] + dDataFlux[i] * np.random.randn(mcIts,len(dataFlux[i]))
        newColors = newFlux[:,1:] / newFlux[:,:-1]
        newColors = np.c_[newColors,newFlux[:,0]/newFlux[:,-1]]
        query = tree.query(newColors,k=1,dualtree=True)
        s = fit_tools.Scale(modelFlux[query[1][:,0]],newFlux,np.ones(np.shape(newFlux)))
        myParams = s
        for j in range(len(modelParams[0])):
            if j in columnsToScale:
                myParams = np.c_[myParams,np.multiply(s,modelParams[query[1][:,0]][:,j])]                
            else:
                myParams = np.c_[myParams,modelParams[query[1][:,0]][:,j]]
        fitParams.append(myParams)
    return(np.array(fitParams))
#------------------------------------------------------------------------------ 
        
def DualTreePeakProbs(data,flagMultimodal=False,saveBBlocks=False): 
    '''
    Creates a Bayesian blocks histogram of the set of values found for each parameter
    for each object. The peak probability value is taken to be the centre of the block
    with highest value.
    Inputs:
        DualTree Output array of size (#objects,mcIts,#params)
    Output:
        Peak probability parameter values for each object, a NumPy array of size(#objects,#parameters)
    '''
    allPeakLocs = [] #for all objects
    multimo = []
    bBlocks = []
    for i in range(len(data)):
        peakLocs = [] #for this object, each parameter
        myMultimo = []
        myBBlocks = []
        for j in range(len(data[0][0])):
            
            bins = bayesian_blocks(data[i][:,j],fitness='events',p0=0.05)
            histo = np.histogram(data[i][:,j],bins)
            # Optional Bayesian Block Histogram storage
            if saveBBlocks:
                myBBlocks.append([bins,histo])
            try:
                nMax = np.argmax(histo[0])
            except:
                print(i,j)
                return(histo)
            loc = (histo[1][nMax]+histo[1][nMax+1])/2.
            peakLocs.append(loc)
            # Optional check for possible multimodality. Not remotely rigorous, but I haven't seen it fail yet
            if flagMultimodal:
                left = histo[0][1:-1] > histo[0][:-2]
                right = histo[0][1:-1] > histo[0][2:]
                nPeaks = np.sum(left*right)
                if nPeaks>1:
                    myMultimo.append(True)
                else:
                    myMultimo.append(False)
        if flagMultimodal:
            multimo.append(myMultimo)
        if saveBBlocks:
            bBlocks.append(myBBlocks)
        allPeakLocs.append(peakLocs)
    if flagMultimodal and not saveBBlocks:
        return([allPeakLocs,multimo])
    if not flagMultimodal and saveBBlocks:
        return([allPeakLocs,bBlocks])
    if flagMultimodal and saveBBlocks:
        return([allPeakLocs,multimo,bBlocks])
    else:
        return(np.array(allPeakLocs))
#------------------------------------------------------------------------------ 
        
def DualTreePeakProbsBins(data,binType='knuth'): 
    '''
    Creates a histogram of the set of values found for each parameter
    for each object. The peak probability value is taken to be the centre of the block
    with highest value.
    Inputs:
        DualTree Output array of size (#objects,mcIts,#params)
    Output:
        Peak probability parameter values for each object, a NumPy array of size(#objects,#parameters)
    '''
    allPeakLocs = [] #for all objects
    for i in range(len(data)):
        peakLocs = [] #for this object, each parameter
        for j in range(len(data[0][0])):
            histo = amlHist(data[i][:,j],bins=binType)
            plt.clf()
            try:
                nMax = np.argmax(histo[0])
            except:
                print(i,j)
                return(histo)
            loc = (histo[1][nMax]+histo[1][nMax+1])/2.
            peakLocs.append(loc)
        allPeakLocs.append(peakLocs)

    return(np.array(allPeakLocs))
#------------------------------------------------------------------------------ 
def DualTreePercentiles(data,percentiles=[25,50,75]):
    '''
    Input:
        DualTree array
    Output:
        Numpy array of size (#percentiles,#objects,#parameters)
        e.g. stats[0] gives the value of the zeroth parameter at the zeroth percentile
        for each object
    '''
    stats = np.percentile(data,percentiles,axis=1)
    return(stats)
#------------------------------------------------------------------------------    

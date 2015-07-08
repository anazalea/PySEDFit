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

'''
Dual Tree, MC

Use MC method instead of creating a new model tree for each object. Create X 
sets of perturbed fluxes for each object. Use a dual tree query to find closest 
model to each (euclidean dist.) Return sampled PDF of each model parameter 
(and covariance matrix?)
'''

#------------------------------------------------------------------------------ 
def DualTree(dataFlux,dDataFlux,modelFlux,modelParams,mcIts):
    '''
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            modelParams = parameters of each model to be recorded, array of size (#models,#parameters)
            mcIts = number of times to perturb fluxes for each object, int
            
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
            nMax = np.argmax(histo[0])
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
    if flagMultimodal:
        return(np.array([allPeakLocs,multimo]))
    else:
        return(np.array(allPeakLocs))
#------------------------------------------------------------------------------ 
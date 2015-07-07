# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 18:47:47 2015

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
"""  
from __future__ import print_function
from __future__ import division

import sys
import math
import numpy as np
from scipy.spatial.distance import seuclidean
from sklearn.neighbors import NearestNeighbors as nn
from sklearn.neighbors import BallTree
#from astroML.plotting import hist as aMLhist
#from astroML.density_estimation import knuth_bin_width
from astroML.density_estimation import bayesian_blocks

###############################################################################
###############################################################################
###############################################################################
'''
Straight Up

These functions do a brute force chi^2 minimization to find the best fit model
'''
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
    if len(np.shape(modelFlux))==1: # there's only one object, need to reshape so that array operations are working on the right indices
        modelFlux = np.reshape(modelFlux,(1,len(modelFlux)))
        dataFlux = np.reshape(dataFlux,(1,len(dataFlux)))
        dDataFlux = np.reshape(dDataFlux,(1,len(dDataFlux)))
    s = np.sum(dataFlux*modelFlux/dDataFlux**2,axis=1) / np.sum(modelFlux**2/dDataFlux**2,axis=1)
    return(s)

     
#------------------------------------------------------------------------------    
def SimpleBrute(dataFluxes,dDataFluxes,modelFluxes):
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
    modelColors = modelFluxes[:,1:] / modelFluxes[:,:-1]
    dataColors = dataFluxes[:,1:]/dataFluxes[:,:-1]
    dDataColors = np.sqrt( (1./dataFluxes[:,:-1])**2 * (dDataFluxes[:,1:])**2 \
                + (dataFluxes[:,1:]/dataFluxes[:,:-1]**2)**2 * (dDataFluxes[:,:-1])**2)
    results = np.array([]).reshape(0,3)
    for i in range(len(dataFluxes)):
        minChi2,n = MinChi2(modelColors,dataColors[i],dDataColors[i])
        s = Scale(modelFluxes[n],dataFluxes[i],dDataFluxes[i])
        results = np.r_[results,[[n,s,minChi2]]]
    return(results)
        
###############################################################################
###############################################################################
###############################################################################
'''
Chi^2 tree
Uses SciKitLearn's nearest neighbor objects with ball tree algorithm and seuclidean metric
(which is equivalent to chi^2 distance when the variance of each observed flux is used as a weight)
A new tree is created for each observation. Probably wouldn't want to use this for real, but 
I'm putting it here for easy sanity checks.
'''

def ChiTree(dataFlux,dDataFlux,modelFlux):
    '''
    Finds the model that minimizes chi^2 distance for each object using a ball_tree search
    with seuclidean metric (weighting each dimension by the variance in that flux)
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            
    Output:
            NumPy array of size (#objects,3)
            Columns [index of model with min chi^2, scale factor, chi^2 value]
    '''
    modelColors = modelFlux[:,1:] / modelFlux[:,:-1]
    dataColors = dataFluxes[:,1:]/dataFluxes[:,:-1]
    dDataColors = np.sqrt( (1./dataFluxes[:,:-1])**2 * (dDataFluxes[:,1:])**2 \
                + (dataFluxes[:,1:]/dataFluxes[:,:-1]**2)**2 * (dDataFluxes[:,:-1])**2)
    results = np.array([]).reshape(0,3)
    for i in range(len(dataFlux)):
        tree = nn(n_neighbors=1,algorithm='ball_tree',metric='seuclidean',metric_params={'V':dDataColors[i]**2})
        tree.fit(modelColors)
        query = tree.kneighbors(dataColors[i],1)
        n,chi2 = query[1][0][0],query[0][0][0]**2.
        s = Scale(modelFlux[n],dataFlux[i],dDataFlux[i])
        results = np.r_[results,[[n,s,chi2]]]
    return(results)
        
        
###############################################################################
###############################################################################
###############################################################################
'''
Dual Tree, MC

Use MC method instead of creating a new model tree for each object. Create X 
sets of perturbed fluxes for each object. Use a dual tree query to find closest 
model to each (euclidean dist.) Return sampled PDF of each model parameter 
(and covariance matrix?)
'''
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
    stats = np.percentile(data,percentiles,axis=2)
    return(stats)
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
        s = Scale(modelFlux[query[1][:,0]],newFlux,np.ones(np.shape(newFlux)))
        myParams = s
        for j in range(len(modelParams[0])):
            myParams = np.c_[myParams,modelParams[query[1][:,0]][:,j]]
        fitParams.append(myParams)
    return(np.array(fitParams))

        
###############################################################################
###############################################################################
###############################################################################
'''
Upper Limits
''' 

def DualTreeWithUpperLimits(dataFlux,dDataFlux,modelFlux,modelParams,mcIts,upperLimFlagVal,upperLimProb):
    '''
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            modelParams = parameters of each model to be recorded, array of size (#models,#parameters)
            mcIts = number of times to perturb fluxes for each object, int
            upperLimFlagVal = where dDataFlux == this value, corresponding elements of dataFlux will br treated as upper limits
            upperLimProb = the probability that upper
            
    Output:
            NumPy array of size (#objects,mcIts,#params)
            e.g. the zeroth element gives you a 2d array where each row represents the
            fit parameters from one monte carlo iteration 
    '''
        
###############################################################################
###############################################################################
###############################################################################
'''
MPI
'''
###############################################################################
###############################################################################
###############################################################################

        
        
<<<<<<< HEAD
modelfile = np.genfromtxt('../../TestFiles/salp_tau1_m20.bbsed')
modelFluxes = 10**(modelfile[:,5:]/-2.5)
modelParams = modelfile[:,:5]
data = np.genfromtxt('../../TestFiles/elliptical.mockobs.cat')
=======
modelfile = np.genfromtxt('../TestFiles/salp_tau1_m20.bbsed')
modelFluxes = 10**(modelfile[:,5:]/-2.5)
modelParams = modelfile[:,:5]
data = np.genfromtxt('../TestFiles/elliptical.mockobs.cat')
>>>>>>> FETCH_HEAD
dataFluxes = 10**(data[:,np.arange(0,len(data[0]),2)]/-2.5)
dDataFluxes = dataFluxes*data[:,1+np.arange(0,len(data[0]),2)]/1.086
dataColors = dataFluxes[:,1:]/dataFluxes[:,:-1]
dDataColors = np.sqrt( (1./dataFluxes[:,:-1])**2 * (dDataFluxes[:,1:])**2 \
                + (dataFluxes[:,1:]/dataFluxes[:,:-1]**2)**2 * (dDataFluxes[:,:-1])**2)


fits = SimpleBrute(dataFluxes,dDataFluxes,modelFluxes)

dtFits = DualTree(dataFluxes,dDataFluxes,modelFluxes,modelParams,1000)


# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:04:22 2015

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
                
Top level module for model fitting
"""

from __future__ import print_function
from __future__ import division

import sys
import os
sys.path.append('.')
sys.path.append('../../')
sys.path.append('../IO/')
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.signal as sig
from scipy import integrate
from scipy.interpolate import interp1d
import MC_dual_tree as mcdt 
from astroML.density_estimation import bayesian_blocks
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

def weighted_percentile(data, percents, weights=None):
    ''' percents in units of 1%
    weights specifies the freqency (count) of data.
    '''
    if weights is None:
        return np.percentile(data, percents)
    ind=np.argsort(data)
    d=data[ind]
    w=weights[ind]
    p=1.*w.cumsum()/w.sum()*100
    y=np.interp(percents, p, d)
    return y

def makePriorGrid(zs,rs,dr,outfile,rmin=15,rmax=30):
    rMids = np.arange(16.,28.,0.1)
    drs = 0.5*np.ones(len(rMids))
    drs[rMids<23]=0.5
    drs[rMids<20]=1.0
    drs[rMids<19]=2.0
    zEdges=np.arange(-0.075,5.075,0.1)
    zMids = (zEdges[1:]+zEdges[:-1])/2.
    allH = []
    for i in range(len(rMids)):
        #print(rMids[i])
        rMsk = np.ma.masked_outside(rs,rMids[i]-dr,rMids[i]+dr)
        zMsk = np.ma.masked_array(zs,mask=np.ma.getmask(rMsk)).compressed()

        h = np.histogram(zMsk,bins=zEdges)[0]
        kernel=np.ones(5)*(1./5.)
        h2=sig.convolve(h,kernel,mode='same')
        h3=sig.convolve(h2,kernel,mode='same')
        g = interp1d(zMids,h3,bounds_error=False, fill_value=0.0)
        tot = integrate.quad(g,0.,7.)
        h3 = h3/tot[0]
        
        if i%5==0:
            plt.plot(zMids,h3,lw=3,alpha=0.75,color=cm.jet(i/len(rMids)),label='r='+str(rMids[i]))
        else:
            plt.plot(zMids,h3,lw=3,alpha=0.5,color=cm.jet(i/len(rMids)))
        allH.append(h3)
    return([rMids,zMids,np.array(allH)])
        
def ProbRemoval(data,h,rMags):
    # data = output from DCDT
    # h = output from makePrior
    # rMags = rMags from catalog, same length as data
    z = []
    bBlocks = []
    for i in range(len(data)):
        rN = np.searchsorted(rMids,rMags[i]) # index of nearest rMag with defined zDist
        if rN > 239:
            print(rMags[i])
            rN = 239
        ws = np.interp(data[i][:,1],zMids,h[rN]) # values of zDist at zFits
        ws = ws/np.max(ws) # Normalize to make probabalistic removal possible
        ran = np.random.rand(len(ws))
        msk = ran > ws # True where ran exceeds prior prob
        zs = np.ma.masked_array(data[i][:,1],mask=msk).compressed()
        bins = bayesian_blocks(zs,fitness='events',p0=0.25)
        histo = np.histogram(zs,bins)
        bBlocks.append([bins,histo])
        try:
            nMax = np.argmax(histo[0])
            loc = (histo[1][nMax]+histo[1][nMax+1])/2.
            z.append(loc)
        except:
            loc = np.percentile(zs,50)
            z.append(loc)
    return(z)
        
        
# Make prior
cat = np.genfromtxt('/Users/anneya/PySEDFit/TestFiles/cmcgals_ugrizy_randomized.cat',usecols=(1,6))
rMids,zMids,h=makePriorGrid(cat[:,0],cat[:,1],0.5,'')
'''
# Read data 
magSoftening = 0.0
cat = np.genfromtxt('/Users/anneya/PySEDFit/TestFiles/LSST_OneYear_simObs_noisier.cat')
ns = np.random.randint(0,len(cat),5000)
cat=cat[ns]
dataMags = cat[:,[4,6,8,10,12,14]]
dDataMags = cat[:,[5,7,9,11,13,15]]
dataFlux = 10.**(dataMags/-2.5)
#dDataFlux += magSoftening
dDataFlux = dataFlux*dDataMags/1.086

# Read models
models=np.genfromtxt('/Users/anneya/PySEDFit/TestFiles/cmc_castor.bbsed')
modelMags = models[:,[12,13,14,15,16,17]]
modelFlux = 10.**(modelMags/-2.5)
modelParams = models[:,[0,1]]


dcdt=mcdt.DaisyChainDualTree(dataFlux,dDataFlux,modelFlux,modelParams,200,columnsToScale=[])

#cols = cm.jet((rMags-17.0041)/8.)

#plt.scatter(z,zC,color='blue',alpha=0.5,label='Corrected',s=20)
'''
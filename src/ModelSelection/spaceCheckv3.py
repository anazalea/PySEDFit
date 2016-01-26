# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 18:15:08 2016

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
sys.path.append('.')
sys.path.append('../IO/')
import numpy as np
import spectrum
from astropy import units as u
from astropy import cosmology 
import param5
import filterv2
from copy import deepcopy 
import matplotlib.pyplot as plt
from skimage import measure
import matplotlib.path as mplPath
import matplotlib.cm as cm
import scipy.signal as sig

plt.rc('text',usetex=False)
plt.rc('font',family='serif')

params=param5.SetSpaceCheckParams('/Users/anneya/PySEDFit/Testfiles/spacecheck.param')

#def main(params):
nModels = params['nmodels']
nModels = 200000
nColours = len(params['filter_names'])-1
colourNames = []
for i in range(nColours):
    colourNames.append(params['filter_names'][i]+' - '+params['filter_names'][i+1])
    
# Read model flux columns
modelFlux = np.genfromtxt(params['model_file'],usecols=(params['model_flux_columns']))

# If we have mags
modelColours = modelFlux[:,1:] - modelFlux[:,:-1]
modelColours = modelColours[range(0,len(modelColours),5)]

# Read data flux columns
dataFlux = np.genfromtxt(params['data_file'],usecols=(params['data_flux_columns']))
# If we have mags
dataColours = dataFlux[:,1:] - dataFlux[:,:-1]
#dataColours = modelFlux[:,1:] - modelFlux[:,:-1]

# Get ranges
#minColour = np.min(np.array([np.min(modelColours),np.min(dataColours)]))
minColour = np.percentile(modelColours,1)#np.min(modelColours)
#maxColour = np.max(np.array([np.max(modelColours),np.max(dataColours)]))
maxColour = np.percentile(modelColours,99)#np.max(modelColours)
pad = 0.1*(maxColour - minColour)

# Set up plot
f, axarr = plt.subplots(nColours-1,nColours-1,sharex='col',sharey='row',figsize=(9,9))
f.subplots_adjust(hspace=0,wspace=0)
for i in range(nColours):
    if 1==1:#i > 0:
        axarr[i-1][0].set_ylabel(colourNames[i],fontsize=8, labelpad=5)
        axarr[i-1][0].set_ylim(minColour-pad,maxColour+pad)
        axarr[i-1][0].set_xlim(minColour-pad,maxColour+pad)
    if i !=  nColours-1:
        axarr[nColours-2][i].set_xlabel(colourNames[i],fontsize=8, labelpad=5)
        axarr[nColours-2][i].set_ylim(minColour-pad,maxColour+pad)
        axarr[nColours-2][i].set_xlim(minColour-pad,maxColour+pad)
f2, axarr2 = plt.subplots(nColours-1,nColours-1,sharex='col',sharey='row',figsize=(9,9))
f2.subplots_adjust(hspace=0,wspace=0)
for i in range(nColours):
    if 1==1:#i > 0:
        axarr2[i-1][0].set_ylabel(colourNames[i],fontsize=8, labelpad=5)
        axarr2[i-1][0].set_ylim(minColour-pad,maxColour+pad)
        axarr2[i-1][0].set_xlim(minColour-pad,maxColour+pad)
    if i !=  nColours-1:
        axarr2[nColours-2][i].set_xlabel(colourNames[i],fontsize=8, labelpad=5)
        axarr2[nColours-2][i].set_ylim(minColour-pad,maxColour+pad)
        axarr2[nColours-2][i].set_xlim(minColour-pad,maxColour+pad)


# Get Nnumber of bins for histogram
if params['nbins']==0:
    nbins = 200
else:
    nbins = params['nbins']


colourChecks = []
# For each pair of colours
for i in range(nColours):
    for j in range(nColours-1):
        if i>j:
            print(i,j)
            h1,xedges,yedges = np.histogram2d(modelColours[:,i],modelColours[:,j],\
                    range=[[minColour-pad,maxColour+pad],[minColour-pad,maxColour+pad]],bins=[nbins,nbins])
            kernel = (1./9.)*np.ones((3,3))
            his = sig.convolve(h1,kernel,mode='same')
            axarr[i-1][j].hist2d(modelColours[:,i],modelColours[:,j],cmap=cm.gray_r,\
                    range=[[minColour-pad,maxColour+pad],[minColour-pad,maxColour+pad]],bins=[nbins,nbins]) 
            contours = measure.find_contours(his, 0.)
            contourLengths = []
            for c in contours:
                contourLengths.append(np.shape(c)[0])
            contourLengths = np.array(contourLengths)
            n = np.argmax(contourLengths)
            contour = contours[n]
            scaledContour = np.c_[(contour[:,1]*(maxColour+pad-(minColour-pad))/nbins)+(minColour-pad),
                                  (contour[:,0]*(maxColour+pad-(minColour-pad))/nbins)+(minColour-pad)]
            axarr[i-1][j].plot(scaledContour[:,1],scaledContour[:,0],color='#5200cc',alpha=0.75,lw=1.5)
            axarr2[i-1][j].plot(scaledContour[:,1],scaledContour[:,0],color='#5200cc',alpha=0.75,lw=1.5)
            axarr2[i-1][j].scatter(modelColours[:,i],modelColours[:,j],color='black',alpha=0.01)
            axarr[i-1][j].scatter(dataColours[:,i],dataColours[:,j],color='#5200cc',alpha=0.3,s=1)
            path = mplPath.Path(scaledContour)
            checkPoints = zip(dataColours[:,j],dataColours[:,i])
            contains = path.contains_points(checkPoints)
            colourChecks.append(contains)
            

            
colourChecks = np.array(colourChecks).astype(int)
checkSum = np.sum(colourChecks,axis=0)/len(colourChecks[:,0])
nOutsiders = len(checkSum[checkSum<1.])

print(str(nOutsiders)+' of '+str(len(dataColours[:,0]))+' galaxies ('+\
        str(round(100.*nOutsiders/len(dataColours[:,0]),2))+'%) outside model space.')
        
if params['data_plot_filename']!=None:
    f.savefig(params['data_plot_filename'])
if params['model_plot_filename']!=None:
    f2.savefig(params['model_plot_filename'])
        





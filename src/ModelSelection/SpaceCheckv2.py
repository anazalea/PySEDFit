# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 13:52:13 2015

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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from skimage import measure
import matplotlib.path as mplPath
import matplotlib.cm as cm
import scipy.signal as sig



plt.rc('text',usetex=False)
plt.rc('font',family='serif')

params=param5.SetSpaceCheckParams('/Users/anneya/PySEDFit/Testfiles/spacecheck.param')
#def main(params):
nColours = len(params['filter_names'])-1
colourNames = []
for i in range(nColours):
    colourNames.append(params['filter_names'][i]+' - '+params['filter_names'][i+1])
    
# Read model flux columns
modelFlux = np.genfromtxt(params['model_file'],usecols=(params['model_flux_columns'])) 
# If we have mags
modelColours = modelFlux[:,1:] - modelFlux[:,:-1]

# Read data flux columns
dataFlux = np.genfromtxt(params['data_file'],usecols=(params['data_flux_columns']))
# If we have mags
dataColours = modelFlux[:,1:] - modelFlux[:,:-1]

minColourX = np.min(np.array([np.min(modelColours[:,:-1]),np.min(dataColours[:,:-1])])) 
maxColourX = np.max(np.array([np.max(modelColours[:,:-1]),np.max(dataColours[:,:-1])]))
minColourY = np.min(np.array([np.min(modelColours[:,1:]),np.min(dataColours[:,1:])]))
maxColourY = np.max(np.array([np.max(modelColours[:,1:]),np.max(dataColours[:,1:])]))

# Set up plot
#f, axarr = plt.subplots(nColours-1,nColours-1,sharex='col',sharey='row',figsize=(9,9))
#f.subplots_adjust(hspace=0,wspace=0)

###############################################################################
nbins = 100
i,j=2,3
padX,padY = 0.1*(maxColourX - minColourX),0.1*(maxColourY - minColourY)


h2,xedges,yedges = np.histogram2d(modelColours[:,i],modelColours[:,j],\
                    range=[[minColourX-padX,maxColourX+padX],[minColourY-padY,maxColourY+padY]],bins=[nbins,nbins])
kernel = (1./9.)*np.ones((3,3))
his = sig.convolve(h2,kernel,mode='same')

plt.pcolormesh(np.arange(len(his)+1),np.arange(len(his)+1),his,cmap=cm.gray_r)
contours = measure.find_contours(his, 0.5)
plt.plot(contours[0][:,1],contours[0][:,0])
if len(contours)>1:
    # Won;t work.
    pass

contourLengths = []
for c in contours:
    contourLengths.append(np.shape(c)[0])
contourLengths = np.array(contourLengths)
n = np.argmax(contourLengths)
contour = contours[n]

path = mplPath.Path(contour)
checkPoints = zip(nbins*(dataColours[:,i]-(minColourX-padX))/(maxColourX+padX-(minColourX-padX)),\
                nbins*(dataColours[:,j]-(minColourY-padY))/(maxColourY+padY-(minColourY-padY)))
#checkPoints = zip(50.*(modelColours[:,j]-(minColourY-padY))/(maxColourY+padY-(minColourY-padY)),\
#                                50.*(modelColours[:,i]-(minColourX-padX))/(maxColourX+padX-(minColourX-padX)))
contains = path.contains_points(checkPoints)
checkPoints = np.array(checkPoints)

plt.scatter(checkPoints[:,1],checkPoints[:,0],color='red',alpha=0.05)

colourChecks = []
###############################################################################
'''
for i in range(nColours):
    for j in range(nColours):
        if i!=j and i>j:
            his,xedges,yedges = np.histogram2d(modelColours[:,i],modelColours[:,j],\
                    range=[[minColourX-padX,maxColourX+padX],[minColourY-padY,maxColourY+padY]],bins=25)
            contours = find_contours(his, 0)
            path = mplPath.Path(contour)
            checkPoints = zip(25.*(dataColours[:,i]-(minColourX-padX))/(maxColourX+padX-(minColourX-padX)),\
                                25.*(dataColours[:,j]-(minColourY-padY))/(maxColourY+padY-(minColourY-padY)))
            contains = path.contains_points(checkPoints)
            colourChecks.append(contains)
            
colourChecks = np.array(colourChecks)
'''    

axarr[1][0].set_xlabel('$\mathrm{M_*}$', fontsize=18, labelpad=5)
axarr[1][1].set_xlabel(r'$\Phi_*\times\,10^4$', fontsize=18, labelpad=5)
axarr[0][0].set_ylabel(r'$\Phi_*\times\,10^4$', fontsize=18, labelpad=5)
axarr[1][0].set_ylabel(r'$\alpha$', fontsize=18, labelpad=5)


'''
i,j=1,2

plt.scatter(modelColours[:,j],modelColours[:,i],color='purple',alpha=0.75,s=1)
plt.xlim(minColourX,maxColourX)
plt.ylim(minColourY,maxColourY)
plt.plot(checkPoints)
'''
'''
for i in range(nColours):
    for j in range(nColours):
        if i!=j and i>j:
            #axarr[i-1][j].scatter(modelColours[:,j],modelColours[:,i],color='black',alpha=0.5,s=2)
            axarr[i-1][j].scatter(dataColours[:,j],dataColours[:,i],color='red',alpha=0.75,s=1)
            axarr[i-1][j].set_xlim(minColourX,maxColourX)
            axarr[i-1][j].set_ylim(minColourY,maxColourY)
            axarr[i-1][0].set_ylabel(colourNames[i],fontsize=10, labelpad=5)
            axarr[nColours-2][j].set_xlabel(colourNames[j],fontsize=10, labelpad=5)

'''
    
'''
if __name__ == "__main__":
    pfile = sys.argv[1]
    args = sys.argv[2:]
    params = param5.SetSpaceCheckParams(pfile,args)
    main(params)
'''
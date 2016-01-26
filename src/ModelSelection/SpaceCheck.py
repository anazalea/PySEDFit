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
f, axarr = plt.subplots(nColours,nColours,sharex='col',sharey='row',figsize=(9,9))
f.subplots_adjust(hspace=0,wspace=0)
'''
axarr[1][0].set_xlabel('$\mathrm{M_*}$', fontsize=18, labelpad=5)
axarr[1][1].set_xlabel(r'$\Phi_*\times\,10^4$', fontsize=18, labelpad=5)
axarr[0][0].set_ylabel(r'$\Phi_*\times\,10^4$', fontsize=18, labelpad=5)
axarr[1][0].set_ylabel(r'$\alpha$', fontsize=18, labelpad=5)
'''

for i in range(nColours):
    for j in range(nColours):
        if i!=j and i>=j:
            axarr[i][j].scatter(modelColours[:,j],modelColours[:,i],color='black',alpha=0.5,s=2)
            axarr[i][j].scatter(dataColours[:,j],dataColours[:,i],color='red',alpha=0.75,s=1)
            axarr[i][j].set_xlim(minColourX,maxColourX)
            axarr[i][j].set_ylim(minColourY,maxColourY)
            axarr[i][0].set_ylabel(colourNames[i],fontsize=10, labelpad=5)
            axarr[nColours-1][j].set_xlabel(colourNames[j],fontsize=10, labelpad=5)


    
'''
if __name__ == "__main__":
    pfile = sys.argv[1]
    args = sys.argv[2:]
    params = param5.SetSpaceCheckParams(pfile,args)
    main(params)
'''
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 13:24:04 2016

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

plt.rc('text',usetex=False)
plt.rc('font',family='serif')

x,y = 2.*np.random.randn(20000),np.random.randn(20000)
his,xedges,yedges = np.histogram2d(x,y,range=[[-10,10],[-10,10]],bins=[20,20])
#plt.hist2d(x+20,y+20,range=[[0,100],[0,100]],bins=[100,100],alpha=0.3)

a = np.zeros((40,40))
a[0:20][:,0:20] = his
#a[a<=1] = 0.
contours = measure.find_contours(a, 0)
contourLengths = []

for c in contours:
    contourLengths.append(np.shape(c)[0])
contourLengths = np.array(contourLengths)
n = np.argmax(contourLengths)
plt.plot(contours[n][:,0],contours[n][:,1],lw=3,alpha=0.4)
plt.xlim(0,40)
plt.ylim(0,40)
plt.gca().set_aspect('equal')
path = mplPath.Path(contours[n])


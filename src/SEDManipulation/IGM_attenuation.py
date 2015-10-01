# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:27:26 2015

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
                
Functions to apply dust reddening curves to your sad blue :( spectra
"""

from __future__ import print_function
from __future__ import division

import sys
sys.path.append('.')
import math
import numpy as np
from astropy import units as u
from copy import deepcopy 

def IGMAttenuateSpectrum(spec,igmLaw,igmOpacity,z):
    '''
    Apply IGM opacity to spectrum
    
    Inputs:
        spec = spectrum object
        dustlaw = string, one of ['Calzetti2000','Calzetti1997','LMC','SMC','MW','Dor30']
        ebv = float in (0,1)
    Ouput:
        Astropy quantity of reddened spectrum of length len(spec.spec), with units of spec.spec
    '''
    #assert isinstance(spec,spectrum.Spectrum)
    assert isinstance(igmOpacity,float)
    igmLaws = ['Calzetti2000','Calzetti1997','LMC','SMC','MW','Dor30']
    igmLawFuncs = [Calzetti2000,Calzetti1997,LMC,SMC,MW,Dor30]
    if dustlaw not in dustLaws:
        raise ValueError(dustlaw,' is not a valid IGM Attenuation Law.')

    nLaw = dustLaws.index(dustlaw)
    specUnit = spec.spec.unit
    reddenedSpec = igmLawFuncs[nLaw](spec.wavelengths.to('micron').value,spec.spec.value,ebv) * u.Unit(specUnit)
    newSpec = spectrum.Spectrum(deepcopy(spec.wavelengths).value,reddenedSpec.value,spec.wavelengths.unit,u.Unit(specUnit),params=deepcopy(spec.params)) 
    newSpec.params['ebv']=ebv
    return(newSpec)    
#------------------------------------------------------------------------------ 
def Madau(lam,lnu,z):
    '''
    Apply the Calzetti 20000 las as defined in Calzetti et al 2000 (ApJ 533 682) Equations 2,3,4
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in A
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''    
'''
lam = spec.wavelengths.value
z=5.0
lymanLimit = 912 * u.Unit('AA')
lyLines = np.array([1216, 1026, 973, 950]).reshape((4,1))
Aj = np.array([0.0036, 0.0017, 0.0012, 0.00093]).reshape((4,1))

lamObs = lam * (1.+z)

# Line Blanketing
tauLines = Aj * (np.divide(np.vstack([lamObs,lamObs,lamObs,lamObs]),lyLines)) ** 3.46
mask = np.vstack([lam<lyLines[0][0],lam<lyLines[1][0],lam<lyLines[2][0],lam<lyLines[3][0]]).astype(float)
cosmicTrans = np.product(np.exp(-tauLines*mask),axis=0)
    
# Photoelectric
    
xEm = 1.+z
zC = lamObs / lymanLimit -1.
mask = lamObs > lymanLimit
'''
#########################################################
lam = spec.wavelengths.value
z = 3.5
wavEm = lam
wavObs = lam * (1.+z)
cosmicTrans = np.ones(len(lam))
wavL = 912.
aJ = [0.0036, 0.0017, 0.0012, 0.00093]
wavJ = [1216., 1026., 973., 950.]

for i in [0,1,2,3]:
    tauThisLine = aJ[i] * ((wavObs)/wavJ[i])**3.46
    mask = wavEm < wavJ[i]
    mask = wavObs < wavJ[i]*(1.+z)
    tauThisLine*=mask.astype(float)
    cosmicTrans*=np.exp(-1.*tauThisLine)
#plt.plot(wavObs[:300],cosmicTrans[:300])

xEm = 1. + z
zC = wavObs/wavL - 1.
xC = 1. + zC

tauPhot = 0.25*xC**3*(xEm**0.46 - xC**0.46) \
        + 9.4*xC**1.5*(xEm**0.18 - xC**0.18) \
        - 0.7*xC**3*(xC**-1.32-xEm**-1.32) \
        - 0.023*(xEm**1.68 - xC**1.68)
        
mask = wavObs<wavL*(1.+z)
tauPhot *= mask.astype(float)
#mask = wavObs<wavL
#tauPhot *= mask.astype(float)
cosmicTrans *= np.exp(-1.*tauPhot)
plt.ylim(0,1)
plt.plot(wavObs[:300],cosmicTrans[:300])


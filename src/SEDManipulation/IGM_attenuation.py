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
def madau(lam,lnu,ebv):
    '''
    Apply the Calzetti 20000 las as defined in Calzetti et al 2000 (ApJ 533 682) Equations 2,3,4
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''    
    r = 4.
    k = 2.659 * (-1.857 + 1.040/lnu) + r
    k[lam<0.63] = 2.659 * (-2.156 + (1.509/lnu[lam<0.63]) - (0.198/lnu[lam<0.63]**2) +\
                (0.011/lnu[lam<0.63]**3)) + r
    newLnu = lnu * 10 ** (-0.4 * ebv * k)
    return(newLnu)
#------------------------------------------------------------------------------ 
def Calzetti2000(lam,lnu,ebv):
    '''
    Apply the Calzetti 20000 las as defined in Calzetti et al 2000 (ApJ 533 682) Equations 2,3,4
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''    


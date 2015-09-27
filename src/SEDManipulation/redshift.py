# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 17:44:08 2015

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
import math
import numpy as np
import spectrum
from astropy import units as u
from astropy import cosmology 
import dust_laws 
'''
To redshift spectrum
    1 - Interstellar dust reddening
    2 - IGM attenuation
    3 - Cosmological redshifting
'''


def RedshiftSpectrum(spec,z,cosmo):
    '''
    Redshift restframe spectrum to observed frame, apply cosmological dimming
    Inputs:
        spec
        z
        cosmology
    Output:
        z_obs = Observed Spectrum 
    '''
    dL = cosmo.luminosity_distance(z).to('m')
    fNuNew = (spec.spec*u.Unit('m')**2)/(4*np.pi*dL**2/(1+z))
    newWavelengths = spec.wavelengths * (1+z)
    return(spectrum.Spectrum(newWavelengths.value,fNuNew.value,newWavelengths.unit,fNuNew.unit))
    
def 
    

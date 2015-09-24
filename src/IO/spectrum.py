# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 14:02:46 2015

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

from __future__ import division, print_function
from astropy import constants as const
from astropy import units as u
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os


class Spectrum:
    def __init__(self,specFile,xUnit,yUnit,columns=[0,1]):
        # check file
        try:
            f = np.genfromtxt(specfile)
        except:
            raise ValueError(specFile,' could not be read.')
        
        # check wavelength units
        try:
            xUnitType = u.get_physical_type(x.Unit(unit))
        except:
            raise TypeError(xUnit," is not a recognized unit. Don't waste my time.")
        if xUnitType not in ['length','frequency']:
            raise TypeError(xUnit,' is neither a unit of frequency nor a unit of length. Get it together.')
        
        # check flux density units
        okFluxUnits = ['Lsun/A']
        if yUnit not in okFluxUnits:
            raise TypeError(yUnit,' is not an accepted unit of flux density.')

        if xUnitType == 'length':    
            self.wavelengths = f[:,0]*u.Unit(xUnit)
            self.wavelengths = self.wavelengths.to('micron')
            
        


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
from filterv2 import Filter 


class Spectrum:
    '''
    Spectrum with wavelength and associated spectral flux densities, stored in
    units of microns and erg/s/cm^2/Hz for ease of use with SED manipulation 
    functions and convolving with filters to get AB mags
    '''
    def __init__(self,x,y,xUnit,yUnit):
        '''
        Inputs:
                x = wavelength/frequency values, 1d array-like
                y = spectral flux density, 1d array-like
                xUnit = astropy.units.Unit of physical type length or frequency
                yUnit = astropy.units.Unit of physical type luminosity/flux density /length or frequency

        '''
        # Validate x input
        wavelengths = np.asarray(x)
        spec = np.asarray(y)
        if wavelengths.ndim!=1:
            raise ValueError('Wavelength/Frequency must be given as a 1d array')
        if wavelengths.shape!=spec.shape:
            raise ValueError('Wavelength/Frequency array and flux density array have different shapes.')
        if u.get_physical_type(xUnit) not in ['frequency','length']:
            raise TypeError(xUnit,' is neither a unit of frequency nor a unit of length. Get it together.')
        wavelengths = wavelengths * xUnit
        wavelengths.to('micron',equivalencies=u.spectral())
        if not np.all(np.ediff1d(wavelengths)>0.):
            if u.get_physical_type(xUnit) == 'frequency':
                raise ValueError('Frequencies must be monotonically decreasing.')
            else:
                raise ValueError('Wavelengths must be monotonically increasing.')
        self.wavelengths = wavelengths
        
        # Validate y input
        spec = np.asarray(y)
        if spec.ndim!=1:
            raise ValueError('Spectrum must be given as a 1d array')
            
        # Define desired flux units from base units
        uFnu = u.Unit('erg')/u.Unit('s')/u.Unit('Hz') 
        uFnuN = u.Unit('erg')/u.Unit('s')/u.Unit('Hz')/u.Unit('m')**2 # astropy equivalency only works when area uncluded
        uFlam = u.Unit('erg')/u.Unit('s')/u.Unit('Angstrom')
        uFlamN = u.Unit('erg')/u.Unit('s')/u.Unit('Angstrom')/u.Unit('m')**2
        
        # Need F_nu, but if there's no area, need to add because too lazy to add equivalency 
        if yUnit.is_equivalent(uFnu):
            spec = spec * yUnit / u.Unit('m')**2
        elif yUnit.is_equivalent(uFlam):
            spec = spec * yUnit/ u.Unit('m')**2
        else:
            spec = spec * yUnit
        
        if not spec.unit.is_equivalent(uFlamN) and not spec.unit.is_equivalent(uFnuN):
            raise ValueError(spec.unit,' not recognized as a unit of spectral flux density.')
        spec = spec.to(uFnuN,equivalencies=u.spectral_density(wavelengths))
        
        if wavelengths[-1]<wavelengths[0]: # then we originally had flux units
            wavelengths = np.flipud(wavelengths)
            spec = np.flipud(spec)
        self.spec = spec
        
    def plot(self):
        x = self.wavelengths.value
        y = self.spec.value
        plt.plot(x,-2.5*np.log10(y))
            
    def convolve(self,filt):
        if not isinstance(filt,Filter):
            raise TypeError('The "filter" you specified is invalid.')
        # Find the part of the spectrum where the FTC is defined
        startN,stopN = np.searchsorted(self.wavelengths.value,[filt.wavelengths[0].value,filt.wavelengths[-1].value])            
        if startN!=0:
            startN-=1
        if stopN!=(len(self.wavelengths)-1):
            stopN+=1
        print(startN,stopN)
        ls = self.wavelengths[startN:stopN].value
        # Resample filter
        newFtc = np.interp(ls,filt.wavelengths.value,filt.ftc)
        num = interp1d(ls,newFtc * self.spec[startN:stopN])
        denom = interp1d(ls,newFtc)
        # Integrate
        numerator = quad(num,filt.wavelengths[0].value,0.999*filt.wavelengths[-1].value)
        denominator = quad(denom,filt.wavelengths[0].value,0.999*filt.wavelengths[-1].value)
        ABmag = -2.5 * np.log10(numerator[0]/denominator[0]) - 48.60
        return(ABmag)
         
         

        
        



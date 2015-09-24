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
from Filter import Filter 


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
            if u.get_phyiscal_type(xUnit) == 'frequency':
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
        if yUnit.is_equivalent(uFlam):
            spec = spec * yUnit/ u.Unit('m')**2
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
        startN,stopN = np.searchsorted(self.wavelengths.value,[filt.ftcLambda[0].value,filt.ftcLambda[-1].value])            
        nus = self.wavelengths[startN-1:stopN+1] # this will break for filters at ends
        # Resample filter
        newFtc = np.interp(nus,filt.ftcNu.value,filt.ftcTransNu)
        num = interp1d(nus,newFtc * specFnu[startN-1:stopN+1])
        denom = interp1d(nus,newFtc)
         
         
                def convolve(self,specNu,specFnu):
        '''
        Convolves filter with a spectrum, returns AB mag
        '''
        # Find the part of the spectrum where the FTC is defined
        startN,stopN = np.searchsorted(specNu.value,[self.ftcNu[0].value,self.ftcNu[-1].value])     
        nus = specNu[startN-1:stopN+1]
        # Resample Filter 
        newFtc = np.interp(nus,self.ftcNu.value,self.ftcTransNu)
        # Define functions to integrate
        num = interp1d(nus,newFtc * specFnu[startN-1:stopN+1])
        denom = interp1d(nus,newFtc)
        # Integrate
        numerator = quad(num,filt.ftcNu[0].value,0.999*filt.ftcNu[-1].value)
        denominator = quad(denom,filt.ftcNu[0].value,0.999*filt.ftcNu[-1].value)
        ABmag = -2.5 * np.log10(numerator[0]/denominator[0]) - 48.60
        return(ABmag)
        
        



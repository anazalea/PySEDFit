# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:42:51 2015

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
 

class Filter:
    '''
    A filter.
    '''
    def __init__(self,ftcFile,unit):
        # Validate unit 
        try:
            unitType = u.get_physical_type(u.Unit(unit))
        except:
            os.system("say '"+unit+" is not a recognized unit. Don't waste my time "+name+"'")
            raise TypeError(unit," is not a recognized unit. Don't waste my time.")
        if unitType not in ['length','frequency']:
            os.system("say '"+unit+" is neither a unit of frequency nor a unit of length. Get it together "+str(name)+"'")
            raise TypeError(unit,' is neither a unit of frequency nor a unit of length. Get it together.')
            
        # Load FTC
        if unitType == 'length':
            self.ftcLambda, self.ftcTransLam = np.loadtxt(ftcFile,unpack=True)
            self.ftcLambda = self.ftcLambda * u.Unit(unit)
            self.ftcLambda = self.ftcLambda.to(u.AA)
            self.ftcNu = const.c / self.ftcLambda
            self.ftcNu = self.ftcNu.to(u.Hz)
            self.ftcNu = self.ftcNu[::-1]
            self.ftcTransLam[self.ftcTransLam<0.] = 0.
            self.ftcTransNu = self.ftcTransLam[::-1]
            
        if unitType == 'frequency':
            self.ftcNu, self.ftcTransNu = np.loadtxt(ftcFile,unpack=True)
            self.ftcNu = self.ftcNu * u.Unit(unit)
            self.ftcNu = self.ftcNu.to(u.Hz)
            self.ftcLambda = const.c / self.ftcNu
            self.ftcLambda = self.ftcLambda.to(u.AA)
            self.ftcTransNu[self.ftcTransNu<0.] = 0.
            self.ftcTransLam = self.ftcTransNu[::-1]
            
        
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
        numerator = quad(num,self.ftcNu[0].value,0.999*self.ftcNu[-1].value)
        denominator = quad(denom,self.ftcNu[0].value,0.999*self.ftcNu[-1].value)
        ABmag = -2.5 * np.log10(numerator[0]/denominator[0]) - 48.60
        return(ABmag)
        
    def plotFilt(self):
        pass

os.system('whoami > namedsflsdfjhdsfjkhdsfjkhdsfjka.txt')
f=open('namedsflsdfjhdsfjkhdsfjkhdsfjka.txt','r')
name=f.readline()
os.system('rm namedsflsdfjhdsfjkhdsfjkhdsfjka.txt')
spec = np.genfromtxt('/Users/anneya/science/GNIRS/2015A/A5V.txt')
ls = spec[:,0] * u.Unit('AA')# AA?
mas=5.
dist = 1/(mas*1.e-3) * u.Unit('parsec')
radius = 3.83174e-8 * u.Unit('parsec')#pc
radius*=2/1.7
flambdas = (radius**2/dist**2)*((((spec[:,1] * u.Unit('erg') )/ u.Unit('s')) / u.Unit('cm') )/u.Unit('cm')/ u.Unit('AA'))

fnus = flambdas * ls**2 / const.c
fnus=fnus[::-1]
nus = const.c / ls
nus = nus.to(u.Hz)[::-1]
    
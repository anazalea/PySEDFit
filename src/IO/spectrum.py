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
    def __init__(self,specFile,xUnit,yUnit):
        # Validate units 
        try:
            xUnitType = u.get_physical_type(x.Unit(unit))
        except:
            os.system("say '"+xUnit+" is not a recognized unit. Don't waste my time "+name+"'")
            raise TypeError(xUnit," is not a recognized unit. Don't waste my time.")
        if xUnitType not in ['length','frequency']:
            os.system("say '"+xUnit+" is neither a unit of frequency nor a unit of length. Get it together "+str(name)+"'")
            raise TypeError(xUnit,' is neither a unit of frequency nor a unit of length. Get it together.')
        
        # How to check flux units?
            
        


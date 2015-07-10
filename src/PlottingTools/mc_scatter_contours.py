# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 20:22:20 2015

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
import math
import numpy as np
import matplotlib.pyplot as plt

def MakeScatterContours(data):
    '''
    
    '''
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')
    nParams = len(data[0][0])
    f, axarr = plt.subplots(nParams,nParams, sharex='col', sharey='row',figsize=(9,9))
    bigax=f.add_subplot(111,zorder=-1)
    bigax.spines['top'].set_color('none')
    bigax.spines['bottom'].set_color('none')
    bigax.spines['left'].set_color('none')
    bigax.spines['right'].set_color('none')
    bigax.tick_params(top='off', bottom='off', left='off', right='off')
    

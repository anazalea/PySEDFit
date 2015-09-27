# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 10:10:47 2015

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
                
Top level module for model fitting
"""

from __future__ import print_function
from __future__ import division

import sys
import os
sys.path.append('../../')
sys.path.append('../IO/')
import math
import numpy as np
import param5

def main(params):
    '''
    Given a param dictionary parsed by param.SetParams, carry out specified model selection procedure
    '''
    print("PySEDFit v0.1\nBlahBlah\nBlah.")
    
    
    #---------------------------------------------------------------------------
    # Read data
    dataFlux = np.genfromtxt(params['data_file'],usecols=(params['data_flux_columns']))
    dDataFlux = np.genfromtxt(params['data_file'],usecols=(params['data_error_columns']))
    dataParamInfo = np.array(params['data_param'])
    if len(np.shape(dataParamInfo))==1:
        dataParamInfo.reshape((1,len(dataParamInfo)))
    dataParams = np.genfromtxt(params['data_file'],usecols=([dataParamInfo[:,0].astype(int)]))

    
    # Apply mag_softening and convert mags to fluxes (if required)
    magSoftening = params['mag_softening']
    if params['data_flux_unit']=='mag':
        dataFlux = 10.**(dataFlux/-2.5)
        dDataFlux += magSoftening
        dDataFlux = dataFlux*dDataFlux/1.086
    # @@@TODO If you have data in flux units, but specify mag_softening in mags, need to add
    # @@@@TODO: if data_flux_unit = an undesireable flux unit, convert with astropy   
    # @@@@TODO: treat upper limits
        
    #---------------------------------------------------------------------------
    # Read models
    modelFlux = np.genfromtxt(params['model_file'],usecols=(params['model_flux_columns']))
    if params['model_flux_unit']=='mag':
        modelFlux = 10.**(modelFlux/-2.5)
    modelParamInfo = np.array(params['model_param'])
    if len(np.shape(modelParamInfo))==1:
        modelParamInfo.reshape((1,len(modelParamInfo)))    
    
    
    modelParams = np.genfromtxt(params['model_file'],usecols=modelParamInfo[:,0].astype(int))
    
    outputOverwrite = params['output_overwrite']
    if outputOverwrite == 'n':
        outputOverwrite = False
    if not outputOverwrite and os.path.isfile(params['output_file']):
        reprimand = 'You told me not to overwrite output but your output file exists. Get it together. Seriously'
        raise IOError(reprimand)    
    #os.system('touch '+params['output_file'])
        
    #---------------------------------------------------------------------------
    # Good to go
        
    #---BRUTE----------------------------------------------------------------------
    if params['fitting_method']=='brute':
        import brute
        brute.BruteFit(params,dataFlux,dDataFlux,dataParams,dataParamInfo,modelFlux,modelParams,modelParamInfo)
      
    
        
if __name__ == "__main__":
    
    pfile = sys.argv[1]
    args = sys.argv[2:]
    params = param5.SetParams(pfile,args)
    main(params)



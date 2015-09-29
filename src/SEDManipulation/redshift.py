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
import numpy as np
import spectrum
from astropy import units as u
from astropy import cosmology 
import dust_laws 
import filterv2
from copy import deepcopy 
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
def make_bbsed():
    pass

#------------------------------------------------------------------------------
def ProcessSpectrum(spec,dustLaw=None,ebvs=None,igmLaw=None,igmOpacities=None,cosmology=None,zs=None,filters=None):
    '''
    Takes a single Spectrum instance and produces all combinations of reddening/redshifting requested, returns spectra or magnitudes
    Inputs:
        spec = base spectrum to be modified, instance of Spectrum
        dustLaw = name of dust law to be applied (or None), one of ['Calzetti2000','Calzetti1997','LMC','SMC','MW','Dor30']
        ebvs = list of E(B-V) values to be applied to base spectrum
        igmLaw = name of IGM attenuation law to be applied , one of []
        igmOpacities = list of IGM opacity values to be applied 
        cosmology = cosmology to use, instance of astropy.cosmology
        zs = list of redshifts to which base spectrum will be shifted
        filters = list of filters to be convolved with modified spectra, if None specified, full spectra are returned
        
    '''    
    # Validate input
        
    
    spectra = [deepcopy(spec)] # spectra will store all of the modified spectra in a list
    # If we're going to add dust/IGM attenuation/redshifting to spectrum and these quantities aren't defined
    # in the original spectrum's params, add 0.0 values to make things easy later    
    if dustLaw!=None and 'ebv' not in spectra[0].params.keys():
        spectra[0].params['ebv']=0.0
    if zs!=None and 'z' not in spectra[0].params.keys():
        spectra[0].params['z']=0.0
        
    # Add requested ISM reddening to base spectrum, each reddened spectrum is appended to spectra
    if dustLaw!=None:
        ebvs = np.asarray(ebvs)
        ebvs = ebvs.reshape((len(ebvs))) 
        for ebv in ebvs:
            newSpec = dust_laws.dustReddenSpectrum(spectra[0],dustLaw,ebv)
            spectra.append(newSpec)
    
    if igmLaw==None and cosmology==None: # Then we just want to return what we already have
        if filters==None: # Then we're returning spectra
            return(spectra)
        else:
            bbseds = []
            for spect in spectra:
                mybbsed = spectrum.BBsed(params=deepcopy(spect.params))
                for filt in filters:
                    mag = spect.convolve(filt)
                    mybbsed.addMag(filt.name,mag)
                bbseds.append(mybbsed)
            return(bbseds)
            
    # If we got here, we either want to redshift and IGM attenuate or just redshift
    zs = np.asarray(zs)
    zs = zs.reshape(len(zs))
    # Calculate dLs?
    if igmLaw!=None:
        igmOpacities = np.asarray(igmOpacities)
        igmOpacities = igmOpacities.reshape((len(igmOpacities)))
    bbseds = []
    for z in zs:
        for spect in spectra:
            newSpec = RedshiftSpectrum(spectra[i],z,cosmology)
            if igmLaw!=None:
                for igmO in igmOpacities:
                    pass
            else:
                pass
                
                    
            
    
    

    return(spectra)
#------------------------------------------------------------------------------    
def OProcessSpectrum(spec,returnSpec=False,zs=None,cosmology=None,ebvs=None,dustLaw=None,igmLaw=None,igmOpacities=None,filters=None):
    '''
    Takes a single Spectrum instance and produces all combinations of reddening/redshifting requested, returns spectra or magnitudes
    '''    
    # Validate input

    # If we're going to add dust/IGM attenuation/redshifting to spectrum and these quantities aren't defined
    # in the original spectrum's params, add 0.0 values to make things easy later    
    if dustLaw!=None and 'ebv' not in spec.params.keys():
        spec.params['ebv']=0.0
        # should also take values of e(b-v)=0 out of ebvs array
    if zs!=None and 'z' not in spec.params.keys():
        spec.params['z']=0.0
    print(spec.params)
    spectra = [deepcopy(spec)]    
    
    # Add requested ISM reddening to base spectrum, each reddened spectrum is appended to spectra
    if dustLaw!=None:
        ebvs = np.asarray(ebvs)
        ebvs = ebvs.reshape((len(ebvs)))
        for ebv in ebvs:
            newSpec = dust_laws.dustReddenSpectrum(spec,dustLaw,ebv)
            spectra.append(newSpec)
            print(newSpec.params['ebv'])

            
    # Add requested IGM attenuation to each dust reddened spectrum in spectra, each attenuated spectrum
    # will also be appended to spectra
    if igmLaw!=None:
        igmOpacities = np.asarray(igmOpacities).reshape((1))
        # Do stuff
    
    # Redshift each spectrum in [spectra] to all requested redshifts and apply cosmological dimming
    # If we're supposed to return full spectra, we'll keep them all in memory, otherwise unnecessary
    if returnSpec: 
        if cosmology!=None:
            zs = np.asarray(zs)
            zs = zs.reshape((len(zs)))
            for z in zs:
                for i in range(len(spectra)):
                    newSpec = RedshiftSpectrum(spectra[i],z,cosmology)
                    spectra.append(newSpec)
        return(spectra)
    
    # If not spectra, we're returning bbSED objects
    bbseds = []
    if cosmology==None: # Then we just wanted restframe spectra convolved
        for spect in spectra:
            mybbsed = spectrum.BBsed(params=spect.params)
            for filt in filters:
                mag = spect.convolve(filt)
                print(mag)
                mybbsed.addMag(filt.name,mag)
            bbseds.append(mybbsed)
        return(bbseds)
    print(len(spectra))
    # We do need to redshift stuff before convolving, but we won't store redshifted spectra
    zs = np.asarray(zs)
    zs = zs.reshape((len(zs)))
    print(zs)
    for z in zs:
        for i in range(len(spectra)):
            newSpec = RedshiftSpectrum(spectra[i],z,cosmology)
            mybbsed = spectrum.BBsed(params=newSpec.params)
            print(newSpec.params)
            for filt in filters:
                mag = newSpec.convolve(filt)
                mybbsed.addMag(filt.name,mag)
            bbseds.append(mybbsed)
    return(bbseds)

#------------------------------------------------------------------------------    
def ProcessSpectrumValidator(spec,returnSpec,zs,cosmology,ebvs,dustLaw,igmLaw,igmOpacities,filters):
    # z=0 not allowed
    # if IGM law specified, cosmology and zs must be specified
    pass
#------------------------------------------------------------------------------
def RedshiftSpectrum(spect,z,cosmo):
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
    print(dL)
    fNuNew = (spect.spec*u.Unit('m')**2)/(4*np.pi*dL**2/(1+z))
    newWavelengths = deepcopy(spect.wavelengths * (1.+z))
    newParams = deepcopy(spect.params)
    newParams['z']=z
    return(spectrum.Spectrum(newWavelengths.value,fNuNew.value,newWavelengths.unit,fNuNew.unit,params=newParams))
    

#------------------------------------------------------------------------------
#Test

cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

uf=np.genfromtxt('../../Testfiles/FTCs/uSDSS.ftc')
uSDSS = filterv2.Filter(uf[:,0],u.Unit('AA'),uf[:,1],'uSDSS')
gf=np.genfromtxt('../../Testfiles/FTCs/gSDSS.ftc')
gSDSS = filterv2.Filter(gf[:,0],u.Unit('AA'),gf[:,1],'gSDSS')
rf=np.genfromtxt('../../Testfiles/FTCs/rSDSS.ftc')
rSDSS = filterv2.Filter(rf[:,0],u.Unit('AA'),rf[:,1],'rSDSS')
iF=np.genfromtxt('../../Testfiles/FTCs/iSDSS.ftc')
iSDSS = filterv2.Filter(iF[:,0],u.Unit('AA'),iF[:,1],'iSDSS')
zf=np.genfromtxt('../../Testfiles/FTCs/zSDSS.ftc')
zSDSS = filterv2.Filter(zf[:,0],u.Unit('AA'),zf[:,1],'zSDSS')
filters = [uSDSS,gSDSS,rSDSS,iSDSS,zSDSS]

f = np.genfromtxt('../../Testfiles/chab_tau0.1_m20.out.cols',usecols=(0,50))
spec = spectrum.Spectrum(f[:,0],f[:,1],u.Unit('Angstrom'),u.Unit('solLum')/u.Unit('Angstrom'))
spec.params['Age']=7.
#newspex=ProcessSpectrum(spec,filters=filters,cosmology=cosmo,zs=np.arange(0,1.,0.2)[1:],dustLaw='Calzetti2000',ebvs=np.arange(0,1.,0.2)[1:])
newspex=ProcessSpectrum(spec,dustLaw='Calzetti2000',ebvs=np.arange(0.1,1.,0.2))

import matplotlib.cm as cm
i=0
for spex in newspex:
    plt.plot(spex.wavelengths.value,-2.5*np.log10(spex.spec.value)-48.6,color=cm.jet(i/len(newspex)))
    i+=1

plt.gca().invert_yaxis()
plt.xlim(0,2000)
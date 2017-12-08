#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 14:45:35 2017

@author: wwaldron
"""

from   __future__ import division

import pysynphot
import numpy      as     np
import pandas     as     pd
from   os         import path as p

# --- Create a pysynphot source from a pandas dataframe -----------------------
def sb99specsrc(fileName, distToSrc=np.sqrt(3/4/np.pi), redshift=0):
    '''Returns the Starburst99 spectrum file as a pysynphot list of sources
    
    Starburst99 spectrum files are not convenient to work with by themselves.
    This function converts the data in the file to a pysynphot list of sources
    to make data crunching easier. The sources can be combined with a filter
    within the pysynphot package to get photometry measurements. The simulation
    years of these sources are also returned as a numpy array for reference.
    
    Parameters
    ----------
    fileName : str
        The Starburst99 spectrum file
    distToSrc : float, optional
        The distance from the observer to the source in [cm]. Do not pass in an
        argument if total flux is desired.
    
    Returns
    -------
    srcs : list
        The each element in the list is a
        pysynphot.spectrum.ArraySourceSpectrum object that corresponds by index
        to the years in yrs. The units of the object are in
        flam=[erg/s/cm/cm/A].
    yrs : ndarray
        The years each Starburst99 simulation calculation was made
    
    '''
    
    # Import data first as pandas dataframe
    # I don't like hard coding these col widths in, but Starburst99 doesn't
    # have a simpler way right now to parse the data
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'
    specdf  = pd.read_fwf(fileName,header=5,widths=[11,17,12,13,26])
    
    # Convert to Array
    specArr = specdf.as_matrix()
    
    # Get the length of one spectrum
    specLen = np.flatnonzero(np.diff(specdf['TIME [YR]']))[0] + 1
    nSpec   = len(np.flatnonzero(np.diff(specdf['TIME [YR]']))) + 1
    nCols   = len(specdf.columns.values)
    
    # Reshape the array to be 3D based on spectrum length
    specArr = specArr.reshape((nSpec,specLen,nCols))
    specArr = np.moveaxis(specArr,(0,1,2),(2,0,1))
    
    # Get Years and Wavelengths since they are common throughout
    yrs   = np.unique(specdf['TIME [YR]'])
    waves = np.unique(specdf['WAVELENGTH [A]'])
    
    # Create a list of sources
    srcs        = nSpec*[None]
    absFluxCor  = (4/3)*np.pi*distToSrc*distToSrc
    for i in np.arange(nSpec):
        absFlux = 10**specArr[:,2,i] # Total Flux in erg / sec / A
        flux    = absFlux / absFluxCor # Total Flux in erg / sec / A / cm / cm
        srcs[i] = pysynphot.ArraySpectrum(wave=waves, flux=flux,
            waveunits='angstrom',fluxunits='flam')
        srcs[i] = srcs[i].redshift(redshift)
    
    return srcs, yrs


# --- Read Equivalent Width Files from SB99 -----------------------------------
def sb99ewidth(fileName,outscale='linear'):
    '''
    '''
    
    # Import data first as pandas dataframe
    # I don't like hard coding these col widths in, but Starburst99 doesn't
    # have a simpler way right now to parse the data
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'
    widthDF = pd.read_fwf(fileName,header=5,widths=([14] + 4*[10,9,9]))
    
    # Change Scales
    if outscale.lower() == 'linear':
        widthDF.iloc[:,1:] = 10**widthDF.iloc[:,1:]
    
    return widthDF
    


# --- From BPASSv2 SED File ---------------------------------------------------
def bpasssedsrc(fileName, distToSrc=np.sqrt(3/4/np.pi), redshift=0):
    '''
    '''
    
    # Check File First
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'
    
    # Declare Years
    n   = np.arange(2,43)
    yrs = 10**(6 + 0.1*(n-2))
    
    # Get Dataframe
    names   = ['wave'] + yrs.tolist()
    specdf  = pd.read_table(fileName, delim_whitespace=True, names=names)
    
    # Get the sources
    srcs       = []
    absFluxCor = (4/3)*np.pi*distToSrc*distToSrc
    for yr in yrs:
        absFlux  = specdf[yr].as_matrix()
        flux     = (absFlux / absFluxCor) * 3.826e33
        srcs.append(pysynphot.ArraySpectrum(wave=specdf['wave'].as_matrix(),
            flux=flux, waveunits='angstrom', fluxunits='flam'))
        srcs[-1] = srcs[-1].redshift(redshift)
    
    return srcs, yrs


# --- From GALEV spec File ----------------------------------------------------
def galevspecsrc(fileName,distToSrc=np.sqrt(3/4/np.pi), redshift=0):
    '''
    '''
    
    # Check File First
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'
    
    # Get commented parameters
    lnNum = 0
    with open(fileName,'r') as f:
        for line in f:
            if lnNum == 1:
                lineList = line.split()[1:]
                yrs      = [float(x) for x in lineList]
                break
            else:
                lnNum += 1
    
    # Get Dataframe
    names   = ['wave'] + yrs
    yrs     = np.array(yrs)
    specdf  = pd.read_table(fileName, delim_whitespace=True, names=names,
                            comment='#', index_col=0)
    
    # Get the sources
    srcs       = []
    absFluxCor = (4/3)*np.pi*distToSrc*distToSrc
    for yr in yrs:
        absFlux  = specdf[yr].as_matrix()
        flux     = absFlux / absFluxCor
        srcs.append(pysynphot.ArraySpectrum(wave=specdf['wave'].as_matrix(),
            flux=flux, waveunits='angstrom', fluxunits='flam'))
        srcs[-1] = srcs[-1].redshift(redshift)
    
    return srcs, yrs
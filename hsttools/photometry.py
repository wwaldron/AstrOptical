#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 08:54:01 2017

@author: wwaldron
"""

# Future Imports
from __future__ import division

# Numerical Imports
import numpy as np

# Astro Imports
from astropy             import units as u
from astropy.wcs         import WCS
from astropy.io          import fits
from astropy.coordinates import SkyCoord
from photutils           import CircularAperture,       CircularAnnulus
from photutils           import SkyCircularAperture, SkyCircularAnnulus
from photutils           import aperture_photometry


# --- Magnitude ---------------------------------------------------------------
def magnitude(flux, zeroPoint, distToSrc=10, galExt=0, apCor=0):
    '''
    '''
    
    # If flux is an array, convert values below zero to NaN
    if isinstance(flux,np.ndarray):
        flux[flux < 0] = np.NaN
    
    # Calculate Magnitude
    mag    = -2.5*np.log10(flux) + zeroPoint - galExt - apCor
    absMag = mag - 5*np.log10(distToSrc/10)
    
    return mag, absMag


# --- Zeropoint Calculations --------------------------------------------------
def zeropoint(photFlam, photPlam):
    '''
    '''
    
    # Calculate Zeropoints
    stMagZPt = -2.5*np.log10(photFlam) - 21.1
    abMagZpt = -2.5*np.log10(photFlam) - 21.1 - 5*np.log10(photPlam) + 18.692
    
    return stMagZPt, abMagZpt


# --- Simple Aperature Photometery --------------------------------------------
def simpleapphot(fileName,pos,r,rI,rO,frame='image'):
    '''
    '''
    
    # Create Aperture
    frame = frame.lower()
    if frame == 'image':
        
        # Create the Aperatures
        cirAps = CircularAperture(pos, r)
        annAps = CircularAnnulus( pos, rI, rO)
        
    elif frame == 'fk5':
        
        # Create Sky Apertures
        pos = SkyCoord(frame=frame, ra=pos[:,0]*u.deg, dec=pos[:,1]*u.deg)
        cirAps = SkyCircularAperture(pos, r *u.arcsec)
        annAps = SkyCircularAnnulus( pos, rI*u.arcsec, rO*u.arcsec)
        
    else:
        raise ValueError('Unsupported coordinate system.')
    
    # Load in the files and do photometry
    hdu        = fits.open(fileName)
    cirPhotTab = aperture_photometry(hdu,cirAps)
    annPhotTab = aperture_photometry(hdu,annAps )
    if frame == 'fk5':
        cirAps = cirAps.to_pixel(WCS(fobj=hdu))
        annAps = annAps.to_pixel(WCS(fobj=hdu))
    hdu.close()
    
    # Get Photometry as ndarray
    phot = cirPhotTab['aperture_sum'].data - \
           (cirAps.area()/annAps.area())*annPhotTab['aperture_sum'].data
    return phot
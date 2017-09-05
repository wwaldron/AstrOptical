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
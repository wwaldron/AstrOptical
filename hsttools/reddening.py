#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 13:41:08 2017

@author: wwaldron
"""

# Future Imports
from __future__ import division

# Meta
__all__ = ['cal97','cal00','car89']

# Python
from warnings import warn

# Numerical Imports
import numpy as np
from   numpy import frompyfunc as vec


# frompyfuncd Functions
def cal97(wave,ebv=0.0):
    '''Calculates reddening according to Calzetti's 1997 Paper
    '''
    
    # Function
    vcal97 = vec(_cal97,2,2)
    
    return vcal97(wave,ebv)

def cal00(wave,ebv=0.0,rvp=4.05):
    '''Calculates reddening according to Calzetti's 2000 Paper
    '''
    
    # Function
    vcal00 = vec(_cal00,3,2)
    
    return vcal00(wave,ebv,rvp)

def car89(wave,ebv=0.0,rvp=3.1):
    '''Calculates reddening according to Cardelli's 1989 Paper
    '''
    
    # Function
    vcar89 = vec(_car89,3,2)
    
    return vcar89(wave,ebv,rvp)


# --- Calzetti 1997 -----------------------------------------------------------
def _cal97(wave, ebv=0.0):
    '''Calculates reddening according to Calzetti's 1997 Paper
    '''
    
    if wave < 0.63:
        if wave < 0.12:
            warn('Reddening is ill defined below 0.12 microns.', RuntimeWarning)
        
        # Calculate k value
        k = 4.88 + 2.656*(-2.156 + 1.509/wave - 0.198/wave/wave +
                             0.011/wave/wave/wave)
    
    else:
        if wave > 1.0:
            warn('Reddening is ill defined above 1.00 microns.', RuntimeWarning)
        
        # Calculate k value
        k = 1.73 - 0.1/wave + 1.86/wave/wave - 0.48/wave/wave/wave
    
    red = k * ebv
    
    return k, red


# --- Calzetti 2000 -----------------------------------------------------------
def _cal00(wave, ebv=0.0, rvp=4.05):
    '''Calculates reddening according to Calzetti's 2000 Paper
    '''
    
    if wave < 0.63:
        if wave < 0.12:
            warn('Reddening is ill defined below 0.12 microns.', RuntimeWarning)
        
        # Calculate k value
        k = rvp + 2.659*(-2.156 + 1.509/wave - 0.198/wave/wave +
                             0.011/wave/wave/wave)
    
    else:
        if wave > 2.20:
            warn('Reddening is ill defined above 2.20 microns.', RuntimeWarning)
        
        # k Value
        k = rvp + 2.659*(-1.857 + 1.04/wave)
        
    red = k * ebv
    
    return k, red


# --- Cardelli 1989 -----------------------------------------------------------
def _car89(wave, ebv=0.0, rvp=3.1):
    '''Cardelli Milky Way reddening (1989)
    '''
    
    # First, get the inverse params
    x    = 1/wave
    
    # Calculate a(x) and b(x)
    if x < 1.1:
        if x < 0.3:
            warn('Reddening is ill defined below 0.3 inverse microns.', RuntimeWarning)
        
        # Get a and b
        a =  0.574*np.power(x,1.61)
        b = -0.527*np.power(x,1.61)
    
    elif x < 3.3:
        
        # Get a and b
        y = x - 1.82
        a = 1 + 0.17699*y - 0.50447*y*y - 0.02427*y*y*y + 0.72085*y*y*y*y + \
            0.01979*y*y*y*y*y - 0.7753*y*y*y*y*y*y + 0.32999*y*y*y*y*y*y*y
        b = 1.41338*y + 2.28305*y*y + 1.07233*y*y*y - 5.38434*y*y*y*y - \
            0.62251*y*y*y*y*y + 5.3026*y*y*y*y*y*y - 2.09002*y*y*y*y*y*y*y
    
    elif x < 8:
        
        # Get Fa and Fb
        if x > 5.9:
            y = x - 5.9
            fA = -0.0447*y*y - 0.009779*y*y*y
            fB =  0.2130*y*y + 0.1207*y*y*y
        else:
            fA = fB = 0
        
        # Get a and b
        a = fA + 1.752 - 0.316*x - 0.104/(0.341 + (x - 4.67)**2)
        b = fB - 3.090 + 1.825*x + 1.206/(0.263 + (x - 4.62)**2)
    
    else:
        if x > 10:
            warn('Reddening is ill defined above 10 inverse microns.', RuntimeWarning)
        
        # Get a and b
        y = x - 8
        a = -1.073 - 0.628*y + 0.137*y*y - 0.070*y*y*y
        b = 13.670 + 4.257*y - 0.420*y*y + 0.374*y*y*y
            
    # Calculate K
    k   = rvp*a + b
    red = k * ebv
    
    return k, red

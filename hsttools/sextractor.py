#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 16:06:35 2017

@author: wwaldron
"""

# Future Imports
from __future__ import division

# Numerics
import numpy as np

# Data Tables
import pandas as pd
import xarray as xr

# DS9 Regions
import pyregion

# Warnings
import warnings


# --- Read Catalog ------------------------------------------------------------
def readcatalog(fileName,ignoreVal=[99,'inf'],returnType='pandas'):
    '''
    '''
    
    # Get the column names
    nVars = prevVarNum = 0
    varNames = []
    with open(fileName) as f:
        for line in f:
            if line[0] == '#':
                col     = line.split()
                varNum  = int(col[1])
                varName = col[2]
                for i in range(varNum - prevVarNum - 1):
                    varNames.append(varNames[-1])
                    nVars += 1
                varNames.append(varName)
                nVars += 1
                prevVarNum = varNum
            else:
                nCols = len(line.split())
                for i in range(nCols - prevVarNum):
                    varNames.append(varNames[-1])
                    nVars += 1
                break
    
    # Load in the table
    cat = pd.read_table(fileName, delim_whitespace=True, header=None,
                        names=varNames, index_col=False, na_values=ignoreVal,
                        comment='#')
    
    # If a numpy array is desired
    if returnType.lower() in ['numpy','array','ndarray']:
        cat = cat.as_matrix()
    elif returnType.lower() == 'dataset':
        cat = cat.to_xarray()
    elif returnType.lower() == 'dataarray':
        cat = xr.DataArray(cat)
        cat = cat.rename({'dim_0':'observation','dim_1':'measurement'})
    
    # Return the Catalog
    return cat


# --- Find sources in regions -------------------------------------------------
def indentifysrcinreg(xSrc,ySrc,regFile):
    '''Returns source mask array for sources that appear in the regions of
    regFile (in image coordinates)
    '''
    
    # Initialize Mask Array
    msk = np.zeros(xSrc.size,dtype=np.bool)
    
    # Read in Regions
    regs = pyregion.open(regFile)
    if regs[0].coord_format != 'image':
        raise ValueError('"image" format only supported region format.')
    
    # Identify Sources
    for reg in regs:
        if reg.name == 'circle':
            x, y, r = reg.coord_list
            inReg = ((xSrc - x)**2 + (ySrc - y)**2 <= r**2)
            msk[inReg] = True
        elif reg.name == 'ellipse':
            x, y, a, b, t = reg.coord_list
            X = (xSrc-x)*np.cos(np.radians(t)) + (ySrc-y)*np.sin(np.radians(t))
            Y = (xSrc-x)*np.sin(np.radians(t)) - (ySrc-y)*np.cos(np.radians(t))
            inReg = ((X/a)**2 + (Y/b)**2 <= 1);
            msk[inReg] = True
        else:
            warnings.warn('Region Type "{}" not recognized.'.format(reg.name),
                          RuntimeWarning)
    
    # Return
    return msk
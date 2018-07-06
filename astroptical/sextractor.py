#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 16:06:35 2017

@author: wwaldron
"""

# Future Imports
from __future__ import division

# Metadata
__all__ = ['readcatalog', 'identifysrcinreg', 'finduncommonsrcs']

# Numerics
import numpy as np

# Data Tables
import pandas as pd
import xarray as xr

# DS9 Regions
import pyregion

# Python
import warnings
from   itertools import combinations as combins
from   os        import path


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
                    if i == 0:
                        dupVarName = varNames[-1]
                        varNames[-1] = dupVarName + '.' +str(i)
                    varNames.append(dupVarName + '.' + str(i + 1))
                    nVars += 1
                varNames.append(varName)
                nVars += 1
                prevVarNum = varNum
            else:
                nCols = len(line.split())
                for i in range(nCols - prevVarNum):
                    if i == 0:
                        dupVarName = varNames[-1]
                        varNames[-1] = dupVarName + '.' +str(i)
                    varNames.append(dupVarName + '.' + str(i + 1))
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
def identifysrcinreg(xSrc,ySrc,regFile):
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
            inReg = ((X/a)**2 + (Y/b)**2 <= 1)
            msk[inReg] = True
        else:
            warnings.warn('Region Type "{}" not recognized.'.format(reg.name),
                          RuntimeWarning)

    # Return
    return msk


# --- Find Uncommon Sources - ImageWise ---------------------------------------
def finduncommonsrcs(catFiles, roiFile=None, maxPixSep=7, writeRegions=True,
                     outDir=None, circRegRadius=5):
    '''Used to identify sources that only occur in one image'''

    # Read in the catalogs
    cats  = []
    nCats = len(catFiles)
    for fileName in catFiles:
        cats.append(readcatalog(fileName))

    # Only Keep sources in ROI
    if roiFile is not None:
        for i, cat in enumerate(cats):
            msk = identifysrcinreg(cat['X_IMAGE'], cat['Y_IMAGE'], roiFile)
            cats[i] = cat.loc[msk, :]

    # Go About identifying spurrious sources
    goodInds = [np.zeros(cat.shape[0], dtype='bool') for cat in cats]
    for i, j in combins(range(nCats), 2):

        # Get XY from each Catalog
        xI = cats[i]['X_IMAGE'].values.reshape(1, -1)
        yI = cats[i]['Y_IMAGE'].values.reshape(1, -1)
        xJ = cats[j]['X_IMAGE'].values.reshape(-1, 1)
        yJ = cats[j]['Y_IMAGE'].values.reshape(-1, 1)

        # Get Pairings
        goodSrcPairs = ((xI - xJ)**2 + (yI - yJ)**2 <= maxPixSep**2)
        goodI, goodJ = np.any(goodSrcPairs, 0), np.any(goodSrcPairs, 1)
        goodInds[i] = np.logical_or(goodInds[i], goodI)
        goodInds[j] = np.logical_or(goodInds[j], goodJ)

    # Spurrious Mask and Region Output
    spurMsk = tuple(np.logical_not(msk) for msk in goodInds)
    for i, fileName in enumerate(catFiles):
        dirName = outDir if outDir is not None else path.dirname(fileName)
        outName = path.join(dirName, path.splitext(path.basename(fileName))[0]
                            + '_uncommon.reg')
        cat  = cats[i].loc[spurMsk[i], :]
        with open(outName, 'w') as f:
            f.write('image\n')
            for x, y in zip(cat['X_IMAGE'], cat['Y_IMAGE']):
                f.write('circle({}, {}, {})\n'.format(x, y, circRegRadius))

    # Return
    return cats, spurMsk

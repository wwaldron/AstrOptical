#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 16:06:35 2017

@author: wwaldron
"""

# Future Imports
from __future__ import division

# Data Tables
import pandas as pd


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
    
    # Return the Catalog
    return cat
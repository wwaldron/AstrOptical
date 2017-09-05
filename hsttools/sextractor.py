#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 16:06:35 2017

@author: wwaldron
"""

# Future Imports
from __future__ import division

# Numerical Imports
import numpy as np

# Data Tables
import pandas as pd


# --- Read Catalog ------------------------------------------------------------
def readcatalog(fileName,ignoreVal=[99,'inf'],returnType='pandas'):
    '''
    '''
    
    # Get the column names
    nVars = varNum = 0
    varNames = []
    with open(fileName) as f:
        for line in f:
            if line[0] == '#':
                # NEED TO RE-THINK
                cols    = line.split()
                varName = cols[2]
                for i in range(cols[1] - varNum):
                    varNames.append(varName)
                    nVars += 1
                varNum = cols[1]
            else:
                
                break
    
    cat = pd.read_table(fileName, delim_whitespace=True, header=None,
                        names=varNames, index_col=0, na_values=ignoreVal,
                        comment='#')
    
    if returnType.lower() == 'numpy' or returnType.lower() == 'ndarray':
        cat = cat.as_matrix()
    
    return cat


if __name__ == '__main__':
    
    tmp = readcatalog('//home/wwaldron/DoctoralResearch/Images/ESO_137-001/'
                'MAST_DATA/WaldronPipeline/SExtractor/Catalogs/'
                'ESO_F275WxF160W.cat')
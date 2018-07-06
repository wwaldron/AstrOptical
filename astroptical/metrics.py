#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Provides various metrics for Astrophysics calculations

@author: wwaldron
"""

# Imports
import numpy as np


# --- Angular Pairwise Distances ----------------------------------------------
def angpdist(X, Y=None):
    '''Gives the pairwise distances for RA/DEC coordinates

    Parameters
    ----------
    X : ndarray
        A 2D numpy array with shape (M,K) where M is the number of observations
        and K is the number of features/dimensions. Note that only the 1st
        columns are used as RA/DEC, respectively.
    Y : ndarray
        A 2D numpy array with shape (N,K) where N is the number of observations
        and K is the number of features/dimensions. Note that only the 1st
        columns are used as RA/DEC, respectively.

    Returns
    -------
    D : ndarray
        A 2D numpy array of the weighted distances of shape (M,M) if Y is None
        or (M,N) if Y is not None.

    '''

    # Add Y if none given = self pairwise distances
    if Y is None: Y = X

    # Get RA and Dec for each
    ra  = {'X' : X[:, 0].reshape((X.shape[0], -1)),
           'Y' : Y[:, 0].reshape((-1, Y.shape[0]))}
    dec = {'X' : X[:, 1].reshape((X.shape[0], -1)),
           'Y' : Y[:, 1].reshape((-1, Y.shape[0]))}

    # Get differences and averages
    dRa  = 3600 * (ra[ 'X'] - ra[ 'Y'])  # ArcSec
    dDec = 3600 * (dec['X'] - dec['Y'])  # ArcSec
    decBar = np.deg2rad( (dec['X'] + dec['Y'])/2 )

    # Get the Result
    return np.sqrt( (dRa * np.cos(decBar))**2 + dDec**2 )


# --- Weighted Euclidean Distances --------------------------------------------
def weuclidean(X, Y=None, wX=None, wY=None):
    '''Computes the weighted Euclidean Values

    Parameters
    ----------
    X : ndarray
        A 2D numpy array with shape (M,K) where M is the number of observations
        and K is the number of features/dimensions
    Y : ndarray
        A 2D numpy array with shape (N,K) where N is the number of observations
        and K is the number of features/dimensions
    wX : ndarray
        A 2D numpy array with shape (M,K) where M is the number of observations
        and K is the number of features/dimensions that establish the weights
        on X
    wX : ndarray
        A 2D numpy array with shape (N,K) where N is the number of observations
        and K is the number of features/dimensions that establish the weights
        on Y

    Returns
    -------
    D : ndarray
        A 2D numpy array of the weighted distances of shape (M,M) if Y is None
        or (M,N) if Y is not None.

    '''

    # None Checks
    if Y  is None: Y  = X.copy()
    if wX is None: wX = np.zeros_like(X)
    if wY is None: wY = np.zeros_like(Y)

    # Assertions
    assert X.ndim == 2, 'X must have two dimensions.'
    assert Y.ndim == 2, 'Y must have two dimensions.'
    assert X.shape == wX.shape, 'X and wX must have the same shape.'
    assert Y.shape == wY.shape, 'Y and wY must have the same shape.'
    assert X.shape[1] == Y.shape[1], 'X and Y must have the same number of features.'

    # Perform the Distance Calculations
    w = 1 + wX[:, np.newaxis, :] + wY[np.newaxis, :, :]
    D = (X[:, np.newaxis, :] - Y[np.newaxis, :, :])/w
    return np.sqrt(np.sum(D*D, 2))

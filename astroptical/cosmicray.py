#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 20:23:19 2017

@author: wwaldron
"""

# Futures
from __future__ import division

# Metadata
__all__ = ['createmask']

# Imports
from os import path as p
from astroscrappy import detect_cosmics
from astropy.io import fits
from numpy import uint8

# create a mask file of the detected cr pixels
def createmask(fileName, sigclip=4.5, sigfrac=0.3, objlim=5.0, gain=1.0,
               readnoise=6.5, satlevel=65536.0, pssl=0.0, niter=4, sepmed=True,
               cleantype='meanmask', fsmode='median', psfmodel='gauss',
               psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False,
               exposure=1.0, fitsExten='SCI'):
    """
    Wraps astroscrappy
    See astroscrappy.pyx for documentation on the parameters
    exposure is if the image is in counts/sec
    """

    # Check file existance
    if not p.exists(fileName):
        raise IOError("File " + fileName + " does not exist")

    # Get the path to the image
    fileWOExt, exten = p.splitext(fileName)
    cleanName = fileWOExt + "_cleaned" + exten

    # Read in the FITS image
    hduListIn   = fits.open(fileName)

    # Change data to e instead of e/s
    inImg  = hduListIn[fitsExten].data*exposure + pssl

    # Create the lacosmic object
    c = detect_cosmics(inImg, sigclip=sigclip, sigfrac=sigfrac,
                       objlim=objlim, gain=gain, readnoise=readnoise,
                       pssl=0., niter=niter,
                       sepmed=sepmed, cleantype=cleantype, fsmode=fsmode,
                       psfmodel=psfmodel, psffwhm=psffwhm, psfsize=psfsize,
                       psfk=psfk, psfbeta=psfbeta, verbose=verbose)

    # Remove the added sky value
    outImg = (c[1] - pssl)/exposure

    # Write out the mask and cleaned images
    primHdr = fits.PrimaryHDU(header=hduListIn['PRIMARY'].header)
    imgHdr     = fits.ImageHDU(data=outImg, header=hduListIn[fitsExten].header,
                               name='SCI')
    mskHdr     = fits.ImageHDU(data=uint8(c[0]), header=hduListIn['SCI'].header,
                               name='MSK')
    hduListOut = fits.HDUList([primHdr,imgHdr,mskHdr])
    hduListOut.writeto(cleanName, output_verify='warn', overwrite=True)
    hduListIn.close()

    return

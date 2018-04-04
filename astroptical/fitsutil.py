#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 20:27:51 2017

@author: wwaldron
"""

# Futures
from __future__ import division

# Metadata
__all__ = ['fixirafcrfix', 'creatermsimage', 'copyheadervalue']

# Imports
from os         import path as p
from astropy.io import fits
from numpy      import sqrt


# --- Fix IRAF's CRFix --------------------------------------------------------
def fixirafcrfix(origFile,crFile):
    """
    IRAF's crfix task only copies the original image's SCI ext to the Primary
    extension. This function moves the headers around and puts things in their
    proper place.
    """

    # Check original file existance
    if not p.exists(origFile):
        raise IOError("File " + origFile + " does not exist")

    # Check cr file existance
    if not p.exists(crFile):
        raise IOError("File " + crFile + " does not exist")

    # Open Files
    origHDUList = fits.open(origFile)
    crHDUList   = fits.open(crFile)

    # Create the Headers
    priHdr = fits.PrimaryHDU(header=origHDUList['PRIMARY'].header)
    imgHdr = fits.ImageHDU(data=crHDUList['PRIMARY'].data,
                           header=crHDUList['PRIMARY'].header,
                           name='SCI')

    # Close FITS
    origHDUList.close()
    crHDUList.close()

    # Create Header List
    outHDUList = fits.HDUList([priHdr,imgHdr])
    outHDUList.writeto(crFile, output_verify='warn', overwrite=True)

    return


# --- Create RMS Image --------------------------------------------------------
def creatermsimage(fileName, hduExt='WHT', scaleFact=1.):
    """
    Takes the square root inverse of a weight image to get the RMS
    """

    # Check file existance
    if not p.exists(fileName):
        print("File " + fileName + " does not exist")
        return

    # Get the path to the image
    fileWOExt, exten = p.splitext(fileName)
    rmsName = fileWOExt + "_rms" + exten

    # Read in the FITS image
    hduListIn   = fits.open(fileName)

    # Get the weight data
    whtData = hduListIn[hduExt].data

    # Create the RMS Image
    rmsData = 1/sqrt(whtData * scaleFact)

    # Create the RMS Image
    primHdr = fits.PrimaryHDU(header=hduListIn['PRIMARY'].header)
    imgHdr  = fits.ImageHDU(data=rmsData, header=hduListIn[hduExt].header,
                               name='RMS')
    hduListOut = fits.HDUList([primHdr,imgHdr])
    hduListOut.writeto(rmsName, output_verify='warn', overwrite=True)
    hduListIn.close()

    return


# --- Copy Header Keyword/Value -----------------------------------------------
def copyheadervalue(fromFile, toFile, extension, keywords):
    '''Copies the keyword/value pair in specified extension from first file to
    second'''

    # Open Files
    ext = extension
    with fits.open(fromFile) as fromHDU:
        with fits.open(toFile, 'update') as toHDU:
            if isinstance(keywords, str):
                toHDU[ext].header[keywords] = fromHDU[ext].header[keywords]
            else:
                for key in keywords:
                    toHDU[ext].header[key] = fromHDU[ext].header[key]

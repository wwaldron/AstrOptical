#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Runs AstroDrizzle for each HST visit to get CR corrected images, then runs
SExtractor to get the sources, then aligns images with TweakReg
See: http://www.stsci.edu/hst/acs/documents/isrs/isr1504.pdf
"""
# TODO: Make work from anywhere

# Python Imports
import os
import sys
import glob
import shutil
import logging
import subprocess as sbpr
from argparse import ArgumentParser
from multiprocessing import cpu_count

# 3rd Party Imports
from stsci.tools import teal
from astropy.io import fits
from drizzlepac.tweakreg import TweakReg
from drizzlepac.astrodrizzle import AstroDrizzle

# Setup
crDir = 'CRClean'


# --- Reset TweakReg and AstroDrizzle -----------------------------------------
teal.unlearn('tweakreg', True)
teal.unlearn('astrodrizzle', True)


# --- Functions ---------------------------------------------------------------
# --- File/ID Association -----------------------------------------------------
def fileassoc(fileList):
    """Associates files according to ASN_ID"""

    assocs = {}

    for fileName in fileList:
        with fits.open(fileName) as fid:
            asnId = fid[0].header['ASN_ID'].upper()
            if asnId not in assocs:
                assocs[asnId] = []
            assocs[asnId].append(fileName)

    return assocs


# --- Drizzle Associated Files ------------------------------------------------
def drizassoc(fAssocs):
    """Drizzles associated files together to get crfix images"""

    # Create the AstroDrizzle Options we want the same
    asOpts = {
        'context' : False,
        'preserve' : False,
        'in_memory' : True,
        'clean' : True,
        'driz_cr_corr' : True,
        'driz_combine' : False,
        'num_cores' : cpu_count()//2
    }

    # Loop Through file Associations and Drizzle
    crFiles = {}
    for key, fileList in fAssocs.items():

        # Run AstroDrizzle
        AstroDrizzle(fileList, **asOpts)

        # Get File Names
        crFiles[key] = [f.replace('flc', 'crclean') for f in fileList]

    # Move the Files to a Directory
    if not os.path.isdir(crDir): os.makedirs(crDir)
    for fileName in glob.iglob('*crclean*'):
        shutil.move(fileName, os.path.join(crDir, fileName))
    shutil.move('astrodrizzle.log', os.path.join(crDir, 'astrodrizzle.log'))

    return crFiles


# --- Get the Image Science HDUs in a Fits File -------------------------------
def _imgsci(fileName):
    """Finds which extensions of a fits file are a SCI Image"""

    # Open Fits Image
    ext, x, nm = [], 'XTENSION', 'EXTNAME'
    with fits.open(fileName) as fid:
        for i, hdu in enumerate(fid):
            hdr = hdu.header
            if (x in hdr and hdr[x] == 'IMAGE' and
                    nm in hdr and hdr[nm] == 'SCI'):
                ext.append(i)
    return ext


# --- Get the SExtractor Command ----------------------------------------------
def _sexcom():

    # Check the Possibilities
    if not sbpr.call('which sextractor', shell=True):
        sCom = 'sextractor'
    elif not sbpr.call('which sex', shell=True):
        sCom = 'sex'
    else:
        raise OSError('Could not find SExtractor on the PATH.')

    return sCom


# --- Parse SExtractor Args ---------------------------------------------------
def _parsesex(argString):
    """Parses the Arguments to be passed to SExtractor"""
    argList = argString.split()

    # Add Dashes
    for i, a in enumerate(argList):
        if not i % 2: argList[i] = '-' + a

    return ' '.join(argList)


# --- Build SExtractor Catalogs -----------------------------------------------
def buildcats(fAssocs, crFiles, sexArgs, shrLoc):
    """Runs SExtractor for all the Files in the glob pattern"""

    # Get the SExtractor Command
    sCom = _sexcom()

    # Create the Parameters file
    pFile = 'tweak.param'
    with open(os.path.join(crDir, pFile), 'w') as fid:
        fid.write('X_WORLD\nY_WORLD\n')

    # Copy Appropriate Files
    for fExt in ['conv', 'nnw', 'sex']:
        if not os.path.isfile(os.path.join(crDir, 'default.' + fExt)):
            shutil.copyfile(os.path.join(shrLoc, 'default.' + fExt),
                            os.path.join(crDir, 'default.' + fExt))

    # Delete the Master File if Catalog Assoc File
    mCat = 'cat_assoc.cat'
    if os.path.isfile(mCat): os.remove(mCat)

    # Loop through Files and Run SExtractor
    for key in crFiles:
        for fName, cName in zip(fAssocs[key], crFiles[key]):

            # Write File Assoc
            with open(mCat, 'a') as fid: fid.write(fName)

            # Get the Base Name of the Catalogs
            catBase = 'tweak_' + os.path.splitext(fName)[0]

            # Get Valid Extensions and SExtractor on Them
            ext = _imgsci(os.path.join(crDir, cName))
            for i in ext:
                catName = catBase + str(i) + '.cat'
                totalCommand = ' '.join([sCom, sexArgs, '-catalog_name',
                                         catName, '-parameters_name', pFile,
                                         "{}'[{:d}]'".format(cName, i)])
                _ = sbpr.run(totalCommand, check=True, cwd=crDir, shell=True)
                with open(mCat, 'a') as fid: fid.write('  ' +
                                                       os.path.join(
                                                           'SExtractor',
                                                           catName))

            # Add a New Line for the Next File to Run
            with open(mCat, 'a') as fid: fid.write('\n')

    # Move Files
    sDir = 'SExtractor'
    if not os.path.isdir(sDir): os.makedirs(sDir)
    for fName in glob.iglob(os.path.join(crDir, 'tweak*.cat')):
        shutil.move(fName, os.path.join(sDir, os.path.basename(fName)))
    shutil.move('cat_assoc.cat', os.path.join(sDir, 'cat_assoc.cat'))


# --- Align Images ------------------------------------------------------------
def alignimgs(fileList, refCat, updateHdr=False, refCols=(1, 2),
              searchRad=1.0):
    """Aligns Images"""

    # Create the AstroDrizzle Options we want the same
    tOpts = {
        'updatehdr' : updateHdr,
        'writecat' : False,
        'interactive' : False,
        'clean' : True,
        'reusename' : True,
        'catfile' : os.path.join('SExtractor', 'cat_assoc.cat'),
        'xyunits' : 'degrees',
        'refcat' : refCat,
        'refxcol' : refCols[0], 'refycol' : refCols[1],
        'searchrad' : searchRad
    }

    # Loop Through file Associations and Drizzle
    # Run AstroDrizzle
    TweakReg(fileList, **tOpts)

    # Clean up PNGs
    for fName in glob.iglob('hist2d*.png'): os.remove(fName)
    for fName in glob.iglob('vector*.png'): os.remove(fName)
    for fName in glob.iglob('residuals*.png'): os.remove(fName)

# --- Final Drizzle -----------------------------------------------------------
def finaldriz(fileList):
    """Drizzles all similar files together"""

    # Set the Options
    asOpts = {
        'context' : False,
        'preserve' : False,
        'clean' : True,
        'build' : True,
        'num_cores' : cpu_count()//2
    }

    # Get File Associations
    fAssocs = {}
    for fName in fileList:
        with fits.open(fName) as fid:
            hdr = fid[0].header
            idStr = '-'.join([hdr['INSTRUME'], hdr['FILTER']])
        if idStr not in fAssocs: fAssocs[idStr] = []
        fAssocs[idStr].append(fName)

    # Run AstroDrizzle
    for idStr in fAssocs:
        AstroDrizzle(fAssocs[idStr], output=idStr, **asOpts)


# --- Main Script -------------------------------------------------------------
if __name__ == '__main__':

    # --- Read in Command Line ------------------------------------------------
    prsr = ArgumentParser(description=__doc__)

    # Fits Files to Read in
    prsr.add_argument('files', nargs='*', type=str,
                      help='A list of file globs to be processed')
    prsr.add_argument('-c', '--clean', action='store_true',
                      help='Flags whether the script should clean up files')
    driz = prsr.add_mutually_exclusive_group()
    driz.add_argument('-d', '--final-drizzle', action='store_true',
                      help='Determines whether to run a final drizzle without '
                           'prompting later')
    driz.add_argument('-D', '--no-final-drizzle', action='store_true')
    prsr.add_argument('-f', '--flat', type=str, default='flc',
                      help='The flat field type (default: flc)')
    prsr.add_argument('-l', '--log', metavar='FILE', type=str,
                      help='Save the log to a file (does not change STDOUT '
                           'verbosity)')
    verb = prsr.add_mutually_exclusive_group()
    verb.add_argument('-v', '--verbose', action='count', default=0,
                      help='Set the verbosity level')
    verb.add_argument('-q', '--quiet', action='store_true',
                      help='Print no status messages')

    # SExtractor
    sex = prsr.add_argument_group('SExtractor')
    sex.add_argument('-a', '--args', metavar='ARGS', type=str, dest='sexArgs',
                     default='',
                     help='Extra SExtractor configuration parameters (in '
                          'quotes; without flag dashes). Example: "filter Y"')
    sex.add_argument('-s', '--share', metavar='DIR',
                     default='/usr/share/sextractor', type=str,
                     help='The location of the default SExtractor helper '
                          'files')

    # TweakReg
    tr = prsr.add_argument_group('TweakReg')
    tr.add_argument('-r', '--ref', metavar='FILE', type=str,
                    help='A catalog according to which each file should be '
                         'aligned (Assumed to be in FK5 WCS)')
    tr.add_argument('-N', '--catcols', metavar='C', type=int, nargs=2,
                    default=[1, 2],
                    help='The columns containing the alignment sources')
    tr.add_argument('-R', '--rad', metavar='X', type=float, default=1.0,
                    dest='searchRad',
                    help='The search radius in arcsec')
    trYN = tr.add_mutually_exclusive_group()
    trYN.add_argument('-y', action='store_true', dest='updateHdr',
                      help='Update WCS without prompting')
    trYN.add_argument('-n', action='store_true', dest='noUpdateHdr',
                      help='Do NOT update WCS and do NOT prompt')

    # Parse
    args = prsr.parse_args(sys.argv[1:])


    # --- Create Loggers and Handlers -----------------------------------------

    # Loggers
    log = logging.getLogger('sextweak')
    log.setLevel(logging.DEBUG)

    # Formatter
    frmt = logging.Formatter('%(asctime)s   %(levelno)-2s   %(message)s')

    # Stream Logger
    if not args.quiet:
        sh = logging.StreamHandler()
        sh.setLevel(logging.INFO - 5*args.verbose)
        sh.setFormatter(frmt)
        log.addHandler(sh)

    # File Logger
    if args.log is not None:
        fh = logging.FileHandler(args.log)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(frmt)
        log.addHandler(fh)


    # --- Get File Associations -----------------------------------------------
    log.info('Associating Files to Visit ID\n')
    fileAs = fileassoc(args.files)

    # --- Drizzle Associations ------------------------------------------------
    log.info('Running AstroDrizzle for Each Association to Mute CRs\n')
    crFileDict = drizassoc(fileAs)

    # --- Run SExtractor ------------------------------------------------------
    log.info('Running SExtractor to get Sources for Alignment\n')
    buildcats(fileAs, crFileDict, _parsesex(args.sexArgs), args.share)

    # --- Run TweakReg --------------------------------------------------------
    log.info('Aligning Images with TweakReg\n')
    updateWCS = args.updateHdr
    alignimgs(args.files, args.ref, updateWCS, args.catcols, args.searchRad)

    # Prompt if Necessary
    if not (args.updateHdr or args.noUpdateHdr):
        ans = input('Would you like to rerun and update the WCS [y,n]? ')
        if ans.lower() == 'y':
            alignimgs(args.files, args.ref, True, args.catcols, args.searchRad)

    # --- Run Final Drizzle ---------------------------------------------------
    if not args.no_final_drizzle:
        if (args.final_drizzle or
                input('Run final drizzle [y,n]? ').lower() == 'y'):
            log.info('Running Final Drizzle\n')
            finaldriz(args.files)

    # --- Run Clean Up --------------------------------------------------------
    if args.clean:
        shutil.rmtree(crDir)
        shutil.rmtree('SExtractor')

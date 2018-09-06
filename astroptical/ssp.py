#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 14:45:35 2017

@author: wwaldron
"""

# Future
from   __future__ import division

# Metadata
__all__ = ['Starburst99Spectrum', 'BpassSpectrum', 'GalevSpectrum']

# Classes
from abc import abstractmethod, ABC

# Other Imports
import pysynphot  as     psp
import numpy      as     np
import scipy      as     sp
import pandas     as     pd
from   os         import path as p
import matplotlib.pyplot as plt
from glob2 import iglob

# Local Imports
from astroptical.metrics import weuclidean


# --- Module Variables --------------------------------------------------------
DEFAULT_DIST = np.sqrt(3/4/np.pi)


# --- Base Class for Spectrums ------------------------------------------------
class SpectrumEvolution(ABC):
    '''Base class to define a generic spectrum evolution'''

    # --- Python Builtins -----------------------------------------------------
    def __init__(self, fileName=None, distToSrc=DEFAULT_DIST,
                 redshift=0):
        #Private
        self.__redshift = None
        self.__distToSrc = None

        # Other Setup
        self.years = None
        self.spectrumList = None

        # Public
        self.fileName  = fileName
        self.distToSrc = distToSrc
        self.redshift  = redshift

        # Read file if necessary
        if self.fileName is not None:
            self.read(self.fileName)

    def __len__(self):
        return len(self.years)

    def __getitem__(self, i):
        return (self.spectrumList[i], self.years[i])

    def __add__(self, other):

        # Error Checking
        if (self.redshift != other.redshift or
                self.distToSrc != other.distToSrc):
            raise ValueError('Cannot concatonate two Spectrums with different '
                             'redshifts or distances to the source.')

        # Concatonate the Years and the Spectrums
        concatYears = np.concatenate((self.years, other.years))
        concatSpecs = self.spectrumList + other.spectrumList

        # Make a new object
        newSpec = type(self)(distToSrc=self.distToSrc, redshift=self.redshift)
        newSpec.fileName = (self.fileName, other.fileName)
        newSpec.spectrumList, newSpec.years = concatSpecs, concatYears
        newSpec.sort()
        return newSpec

    # --- Set/Get Methods -----------------------------------------------------
    @property
    def redshift(self):
        return self.__redshift

    @redshift.setter
    def redshift(self, value):

        # Set the Redshift if not the Same as Current
        if self.__redshift != value:
            self.__redshift = value

            # Set the Redshift of the Spectrum
            if self.spectrumList is not None:
                for i, spec in enumerate(self.spectrumList):
                    self.spectrumList[i] = spec.redshift(self.redshift)

    @property
    def distToSrc(self):
        return self.__distToSrc

    @distToSrc.setter
    def distToSrc(self, value):

        # Set the Distance if not the same as current
        if self.__distToSrc != value:
            oldVal = DEFAULT_DIST if self.distToSrc is None else self.distToSrc
            self.__distToSrc = value

            # Set the Distance to the Spectrum
            if self.spectrumList is not None:
                oldVal **= 2
                oldVal  *= 4*np.pi/3
                newVal   = (4*np.pi/3)*value**2
                for i, spec in enumerate(self.spectrumList):
                    waves, flux = spec.getArrays()
                    flux *= oldVal/newVal
                    spec = psp.ArraySpectrum(wave=waves, flux=flux,
                                             waveunits=spec.waveunits,
                                             fluxunits=spec.fluxunits)
                    self.spectrumList[i] = spec.redshift(self.__redshift)


    # --- Abstract Methods ----------------------------------------------------
    @abstractmethod
    def read(self, fileName):
        '''To be implemented by child classes to read in necessary file

        Needs to add `spectrumList` and `years` attributes to the object

        '''
        pass

    # --- Public Methods ------------------------------------------------------
    def sort(self):

        # Get the Sorting Necessities
        srtInds = np.argsort(self.years)

        # Sort
        self.years = self.years[srtInds]
        outSpecs = []
        for ind in srtInds:
            outSpecs.append(self.spectrumList[ind])
        self.spectrumList = outSpecs

        return

    def plottrack(self, cMagX, cMagY, ax=None, outUnit='ABMag', **kwargs):
        '''Plots the color/magnitude diagram track over the years

        Parameters
        ----------
        cMagX : iterable
            A Pysynphot Bandpass object or iterable of 2 Bandpass objects in
            which x will be a color between the first and second elements (i.e.
            cMagX[0] - cMagX[1])
        cMagY : iterable
            A Pysynphot Bandpass object or iterable of 2 Bandpass objects in
            which y will be a color between the first and second elements (i.e.
            cMagY[0] - cMagY[1])

        Returns
        -------
        ax : pyplot.Axes
            The pointer to the axes object

        '''

        # Get measurements
        if isinstance(cMagX, (tuple, list)):
            measX = self._getcolor(cMagX, outUnit)
        else:
            measX = self.makeobservation(cMagX, outUnit)
        if isinstance(cMagY, (tuple, list)):
            measY = self._getcolor(cMagY)
        else:
            measY = self.makeobservation(cMagY, outUnit)

        # Make plot
        if ax is None:
            plt.figure()
            ax = plt.axes()

        return ax.plot(measX, measY, **kwargs)

    def makeobservation(self, observer, outUnit='ABMag'):
        '''Makes an observation of the spectrum through the years.'''

        # Check
        if not (isinstance(observer, psp.spectrum.SpectralElement) or
                hasattr(observer, 'throughput')):
            raise TypeError('observer must be a pysynphot observer')
        outUnit = outUnit.lower()

        # Make Observation measurments
        meas = []
        for spec, _ in self:
            meas.append(psp.Observation(spec, observer).effstim(outUnit))
        return np.array(meas)

    def _getcolor(self, colorFilts, outUnit='ABMag'):
        return self.makeobservation(colorFilts[0], outUnit) - \
               self.makeobservation(colorFilts[1], outUnit)

    def calcobsage(self, colorFilts, obs, err=None, nInterpYrs=1000):
        '''Calculates the age of the observations based on the model

        Parameters
        ----------
        colorFilts : iterable
            An iterable of arbitrary length that contains 2-tuples of
            PySynPhot Bandpass filters.
        obs : ndarray
            An M x N array where M is the number of color observations and N is
            the same as the length of color filt pairs.
        err : ndarray
            An M x N array where M is the number of error observations and N is
            the same as the length of color filt pairs.

        Returns
        -------
        ages : ndarray

        '''

        # Checks
        if len(colorFilts) != obs.shape[1]:
            raise ValueError('The number of columns (colors) in obs must '
                             'match the number of filter pairs in colorFilts.')
        if err is not None:
            if len(colorFilts) != obs.shape[1]:
                raise ValueError('The number of columns (colors) in err must '
                                 'match the number of filter pairs in '
                                 'colorFilts.')

        # Get the Model colors
        modColors = np.array([self._getcolor(filts) for filts in colorFilts]).T

        # Interpolate the model colors for finer resolution
        intYrs = np.logspace(np.log10(self.years.min()),
                             np.log10(self.years.max()), nInterpYrs)
        interpolator = sp.interpolate.interp1d(self.years, modColors, axis=0,
                                               bounds_error=False,
                                               fill_value='extrapolate')
        modColors = interpolator(intYrs)

        # Now find the matches
        dists = weuclidean(modColors, obs, wY=err)
        minDistInds = np.nanargmin(dists, 0)
        return intYrs[minDistInds]



# --- Starburst99 Spectrum ----------------------------------------------------
class Starburst99Spectrum(SpectrumEvolution):
    '''Class to hold Starburst99 Runs'''

    def read(self, fileName):
        '''Reads the Starburst99 File'''

        self.fileName = fileName
        self.spectrumList, self.years = sb99specsrc(self.fileName,
                                                    self.distToSrc,
                                                    self.redshift)


# --- BPASS Spectrum ----------------------------------------------------------
class BpassSpectrum(SpectrumEvolution):
    '''Class to hold BPASS runs'''

    def read(self, fileName):
        '''Reads the BPASS file'''

        self.fileName = fileName
        self.spectrumList, self.years = bpasssedsrc(self.fileName,
                                                    self.distToSrc,
                                                    self.redshift)


# --- GALEV Spectrum ----------------------------------------------------------
class GalevSpectrum(SpectrumEvolution):
    '''Class to hold BPASS runs'''

    def read(self, fileName):
        '''Reads the BPASS file'''

        self.fileName = fileName
        self.spectrumList, self.years = galevspecsrc(self.fileName,
                                                     self.distToSrc,
                                                     self.redshift)


# --- MappingsV Spectrum ------------------------------------------------------
class MappingsVSpectrum(SpectrumEvolution):
    '''Class to hold MappingsV Output'''

    # --- Python Builtins -----------------------------------------------------
    def __init__(self, fileName=None, distToSrc=DEFAULT_DIST,
                 redshift=0, years=None, sbSpec=None):

        # Run SpecEv init without a fileName to load all the meta-data
        # then add the years then load file if exists.
        super(MappingsVSpectrum, self).__init__(None, distToSrc, redshift)
        self.years = years
        self.sbSpec = sbSpec
        if fileName is not None:
            self.read(fileName, sbSpec)

    # -- Concretions ----------------------------------------------------------
    def read(self, fileName, sbSpec=None):
        '''Reads a **directory** of MappingsV files to get the spectrum
        evolution

        Parameters
        ----------
        fileName : str
            A glob pattern for the NFN files

        '''

        self.fileName = fileName
        self.sbSpec = sbSpec
        self.spectrumList = mapspecsrc(fileName, self.distToSrc,self.redshift,
                                       sbSpec)


# --- Create a pysynphot source from a pandas dataframe -----------------------
def sb99specsrc(fileName, distToSrc=DEFAULT_DIST, redshift=0):
    '''Returns the Starburst99 spectrum file as a pysynphot list of sources

    Starburst99 spectrum files are not convenient to work with by themselves.
    This function converts the data in the file to a pysynphot list of sources
    to make data crunching easier. The sources can be combined with a filter
    within the pysynphot package to get photometry measurements. The simulation
    years of these sources are also returned as a numpy array for reference.

    Parameters
    ----------
    fileName : str
        The Starburst99 spectrum file
    distToSrc : float, optional
        The distance from the observer to the source in [cm]. Do not pass in an
        argument if total flux is desired.

    Returns
    -------
    srcs : list
        The each element in the list is a
        pysynphot.spectrum.ArraySourceSpectrum object that corresponds by index
        to the years in yrs. The units of the object are in
        flam=[erg/s/cm/cm/A].
    yrs : ndarray
        The years each Starburst99 simulation calculation was made

    '''

    # Import data first as pandas dataframe
    # I don't like hard coding these col widths in, but Starburst99 doesn't
    # have a simpler way right now to parse the data
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'
    specdf  = pd.read_fwf(fileName,header=5,widths=[11,17,12,13,26])

    # Convert to Array
    specArr = specdf.values

    # Get the length of one spectrum
    specLen = np.flatnonzero(np.diff(specdf['TIME [YR]']))[0] + 1
    nSpec   = len(np.flatnonzero(np.diff(specdf['TIME [YR]']))) + 1
    nCols   = len(specdf.columns.values)

    # Reshape the array to be 3D based on spectrum length
    specArr = specArr.reshape((nSpec,specLen,nCols))
    specArr = np.moveaxis(specArr,(0,1,2),(2,0,1))

    # Get Years and Wavelengths since they are common throughout
    yrs   = np.unique(specdf['TIME [YR]'])
    waves = np.unique(specdf['WAVELENGTH [A]'])

    # Create a list of sources
    srcs        = nSpec*[None]
    absFluxCor  = (4/3)*np.pi*distToSrc*distToSrc
    for i in np.arange(nSpec):
        absFlux = 10**specArr[:,2,i] # Total Flux in erg / sec / A
        flux    = absFlux / absFluxCor # Total Flux in erg / sec / A / cm / cm
        srcs[i] = psp.ArraySpectrum(wave=waves, flux=flux,
                                    waveunits='angstrom',
                                    fluxunits=psp.units.Flam)
        srcs[i] = srcs[i].redshift(redshift)

    return srcs, yrs


# --- Read Equivalent Width Files from SB99 -----------------------------------
def sb99ewidth(fileName,outscale='linear'):
    '''
    '''

    # Import data first as pandas dataframe
    # I don't like hard coding these col widths in, but Starburst99 doesn't
    # have a simpler way right now to parse the data
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'
    widthDF = pd.read_fwf(fileName,header=5,widths=([14] + 4*[10,9,9]))

    # Change Scales
    if outscale.lower() == 'linear':
        widthDF.iloc[:,1:] = 10**widthDF.iloc[:,1:]

    return widthDF


# --- From BPASSv2 SED File ---------------------------------------------------
def bpasssedsrc(fileName, distToSrc=DEFAULT_DIST, redshift=0):
    '''Returns the BPASS spectrum file as a pysynphot list of sources

    Parameters
    ----------
    fileName : str
        The Starburst99 spectrum file
    distToSrc : float, optional
        The distance from the observer to the source in [cm]. Do not pass in an
        argument if total flux is desired.
    redshift : float
        Redhsift to the source (z)

    Returns
    -------
    srcs : list
        The each element in the list is a
        pysynphot.spectrum.ArraySourceSpectrum object that corresponds by index
        to the years in yrs. The units of the object are in
        flam=[erg/s/cm/cm/A].
    yrs : ndarray
        The years each Starburst99 simulation calculation was made

    '''

    # Check File First
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'

    # Declare Years
    n   = np.arange(2,43)
    yrs = 10**(6 + 0.1*(n-2))

    # Get Dataframe
    names   = ['wave'] + yrs.tolist()
    specdf  = pd.read_table(fileName, delim_whitespace=True, names=names)

    # Get the sources
    srcs       = []
    absFluxCor = (4/3)*np.pi*distToSrc*distToSrc
    for yr in yrs:
        absFlux  = specdf[yr].values
        flux     = (absFlux / absFluxCor) * 3.826e33
        srcs.append(psp.ArraySpectrum(wave=specdf['wave'].values,
            flux=flux, waveunits='angstrom', fluxunits=psp.units.Flam))
        srcs[-1] = srcs[-1].redshift(redshift)

    return srcs, yrs


# --- From GALEV spec File ----------------------------------------------------
def galevspecsrc(fileName,distToSrc=DEFAULT_DIST, redshift=0):
    '''Returns the GALEV spectrum file as a pysynphot list of sources

    Parameters
    ----------
    fileName : str
        The Starburst99 spectrum file
    distToSrc : float, optional
        The distance from the observer to the source in [cm]. Do not pass in an
        argument if total flux is desired.

    Returns
    -------
    srcs : list
        The each element in the list is a
        pysynphot.spectrum.ArraySourceSpectrum object that corresponds by index
        to the years in yrs. The units of the object are in
        flam=[erg/s/cm/cm/A].
    yrs : ndarray
        The years each Starburst99 simulation calculation was made

    '''

    # Check File First
    assert isinstance(fileName, str), 'fileName must be a string.'
    assert p.exists(fileName),        'File ' + fileName + 'does not exist.'

    # Get commented parameters
    lnNum = 0
    with open(fileName,'r') as f:
        for line in f:
            if lnNum == 1:
                lineList = line.split()[1:]
                yrs      = [float(x) for x in lineList]
                break
            else:
                lnNum += 1

    # Get Dataframe
    names   = ['wave'] + yrs
    yrs     = np.array(yrs)
    specdf  = pd.read_table(fileName, delim_whitespace=True, names=names,
                            comment='#', index_col=0)

    # Get the sources
    srcs       = []
    absFluxCor = (4/3)*np.pi*distToSrc*distToSrc
    for yr in yrs:
        absFlux  = specdf[yr].values
        flux     = absFlux / absFluxCor
        srcs.append(psp.ArraySpectrum(wave=specdf['wave'].values,
            flux=flux, waveunits='angstrom', fluxunits=psp.units.Flam))
        srcs[-1] = srcs[-1].redshift(redshift)

    return srcs, yrs


# --- Read the MappingsV Files ------------------------------------------------
def mapspecsrc(fileName, distToSrc=DEFAULT_DIST, redshift=0, sbSpec=None):
    '''Returns the MappingsV spectrum files as a pysynphot list of sources

    Parameters
    ----------
    dirName : str
        The glob pattern for the NFN files.
    distToSrc : float, optional
        The distance from the observer to the source in [cm]. Do not pass in an
        argument if total flux is desired.

    Returns
    -------
    srcs : list
        The each element in the list is a
        pysynphot.spectrum.ArraySourceSpectrum object that corresponds by index
        to the years in yrs. The units of the object are in
        flam=[erg/s/cm/cm/A].

    '''

    # Get the File Iterator
    fileIt = iglob(fileName)

    # Read in the Data
    srcs = []
    absFluxCor = (4/3)*np.pi*distToSrc*distToSrc
    for fileName in sorted(fileIt):

        # Load Mappings File
        mapDF = pd.read_table(fileName, delim_whitespace=True, header=None,
                              names=['wave', 'flux'], skiprows=9, skipfooter=1,
                              engine='python')
        mpFnu = mapDF['flux'].values/mapDF['wave'].values  # nuFnu --> Fnu
        mpFnu /= absFluxCor  # May have to remove since MP already defines
        mpSpec = psp.ArraySpectrum(mapDF['wave'].values, mpFnu,
                                   waveunits=psp.units.Hz,
                                   fluxunits=psp.units.Fnu).redshift(redshift)
        mpSpec.convert(psp.units.Angstrom)
        mpSpec.convert(psp.units.Flam)
        srcs.append(mpSpec)

    # MappingsV changes the scale of the spectrum (compared to  Starburst99)
    # for some reason. Therefore, this is an attempt to put them on the same
    # scale if they should be.
    if sbSpec is not None:
        if sbSpec.redshift != redshift:
            raise ValueError('Starburst99 spectrum must have the same redshift')
        pivotFreq = 4.9782e+14
        for i, (mpSrc, sbSrc) in enumerate(zip(srcs, sbSpec.spectrumList)):

            # Convert Both
            mpSrc.convert(psp.units.Hz); mpSrc.convert(psp.units.Fnu)
            sbSrc.convert(psp.units.Hz); sbSrc.convert(psp.units.Fnu)

            # Fix the Flux
            mFact = sbSrc.sample(pivotFreq)/mpSrc.sample(pivotFreq)
            mpWave, mpFlux = mpSrc.getArrays()
            mpFlux *= mFact
            mpSrc = psp.ArraySpectrum(mpWave, mpFlux,
                                      waveunits=psp.units.Hz,
                                      fluxunits=psp.units.Fnu)

            # Convert Back
            mpSrc.convert(psp.units.Angstrom); mpSrc.convert(psp.units.Flam)
            sbSrc.convert(psp.units.Angstrom); sbSrc.convert(psp.units.Flam)
            srcs[i] = mpSrc

    return srcs

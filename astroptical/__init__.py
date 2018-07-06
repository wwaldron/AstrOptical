#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Metadata
__all__ = ['cosmicray','fitsutil','photometry','sextractor','ssp',
           'reddening']
import pkg_resources
__version__ = pkg_resources.get_distribution('astroptical').version
del pkg_resources

# Local Impots
from . import cosmicray
from . import fitsutil
from . import photometry
from . import sextractor
from . import ssp
from . import reddening
from . import metrics

# Astropy Imports
#from astropy.io import fits as fitsio  # Not sure whether ethical or good practice

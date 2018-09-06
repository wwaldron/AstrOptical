#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''A toolkit to assist in visible astronomy'''

# Metadata
__all__ = ['cosmicray','fitsutil','photometry','sextractor','ssp',
           'reddening']
import pkg_resources
__version__ = pkg_resources.get_distribution('astroptical').version
del pkg_resources

# Local Impots
from astroptical import metrics
from astroptical import cosmicray
from astroptical import fitsutil
from astroptical import photometry
from astroptical import sextractor
from astroptical import ssp
from astroptical import reddening

# Astropy Imports
#from astropy.io import fits as fitsio  # Not sure whether ethical or good practice

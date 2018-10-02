#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Install script to install astroptical'''

# Imports
from setuptools import setup, find_packages

# Setup command
setup(name='astroptical',
      version='0.1.0',
      description='Tools useful in Optical Astronomy data processing.',
      author='Will Waldron',
      author_email='wvw0330@uah.edu',
      packages=find_packages(),
      scripts=['bin/sextweak.py'],
      zip_safe=True)

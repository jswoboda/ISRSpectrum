#!/usr/bin/env python
"""
setup.py
This is the setup file for the RadarDataSim python package

@author: John Swoboda
"""
import os, inspect
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Creates ISR Spectrums',
    'author': 'John Swoboda',
    'url': 'github.com/jswoboda/ISRSpectrum',
    'author_email': 'swoboj@bu.edu',
    'version': '1.0',
    'install_requires': ['numpy', 'scipy', 'tables'],
    'packages': ['ISRSpectrum'],
    'scripts': [],
    'name': 'ISRSpectrum'
}

curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
testpath = os.path.join(curpath,'Test')
if not os.path.exists(testpath):
    os.mkdir(testpath)
    print "Making a path for testing at "+testpath
setup(**config)

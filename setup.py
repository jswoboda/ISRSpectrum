#!/usr/bin/env python
"""
setup.py
This is the setup file for the RadarDataSim python package

@author: John Swoboda
"""
import os,subprocess
from setuptools import setup

try:
    subprocess.call(['conda','install','--yes','--file','requirements.txt'])
except Exception as e:
    pass

config = {
    'description': 'Creates ISR Spectrums',
    'author': 'John Swoboda',
    'url': 'github.com/jswoboda/ISRSpectrum',
    'install_requires': ['PythonISRUtilties'],
    'dependency_links': ['https://github.com/jswoboda/PythonISRUtilties/tarball/master#egg=PythonISRUtilties'],
    'version': '1.0',
    'packages': ['ISRSpectrum'],
    'name': 'ISRSpectrum'
}

curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
testpath = os.path.join(curpath,'Test')
if not os.path.exists(testpath):
    os.mkdir(testpath)
    print("Making a path for testing at "+testpath)
setup(**config)

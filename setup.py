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
    'install_requires': ['isrutilities'],
    'dependency_links': ['https://github.com/jswoboda/PythonISRUtilities/tarball/master#egg=PythonISRUtilities'],
    'version': '1.0',
    'packages': ['ISRSpectrum'],
    'name': 'ISRSpectrum'
}

curpath = os.path.dirname(__file__)
testpath = os.path.join(curpath,'Test')
try:
    os.mkdir(testpath)
    print("created {}".format(testpath))
except OSError:
    pass
    
setup(**config)

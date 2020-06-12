#!/usr/bin/env python
"""
setup.py
This is the setup file for the RadarDataSim python package

@author: John Swoboda
"""
req = ['nose','six','numpy','scipy','pathlib2','pandas']

import os
from setuptools import setup,find_packages

config = dict(
    description='Creates ISR Spectrums',
    author='John Swoboda',
    url='https://github.com/jswoboda/ISRSpectrum',
    install_requires=req,
    python_requires='>=2.7',
    extras_require={'plot':['matplotlib','jupyter','seaborn'],},
    version='2.0.1',
    packages= find_packages(),
    name= 'ISRSpectrum')

curpath = os.path.dirname(__file__)
testpath = os.path.join(curpath,'Test')
try:
    os.mkdir(testpath)
    print("created {}".format(testpath))
except OSError:
    pass

setup(**config)

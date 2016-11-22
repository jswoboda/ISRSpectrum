#!/usr/bin/env python
import os
from setuptools import setup

config = {
    'description': 'Creates ISR Spectrums',
    'author': 'John Swoboda',
    'url': 'github.com/jswoboda/ISRSpectrum',
    'install_requires': ['PythonISRUtilties'],
    'dependency_links': ['https://github.com/jswoboda/PythonISRUtilties/tarball/master#egg=PythonISRUtilties'],
    'version': '1.0',
    'packages': ['ISRSpectrum'],
    'scripts': [],
    'name': 'ISRSpectrum'
}

curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
testpath = os.path.join(curpath,'Test')
if not os.path.exists(testpath):
    os.mkdir(testpath)
    print("Making a path for testing at "+testpath)
setup(**config)

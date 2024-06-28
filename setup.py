#!/usr/bin/env python
"""
setup.py
This is the setup file for the RadarDataSim python package

@author: John Swoboda
"""
with open('requirements.txt') as f:
    req = f.read().splitlines()

with open('README.md', 'r') as f:
    long_desc = f.read()

import versioneer
from setuptools import setup, find_packages

config = dict(
    description="Creates ISR Spectrums",
    long_description=long_desc,
    long_description_content_type='text/markdown',
    author="John Swoboda",
    url="https://github.com/jswoboda/ISRSpectrum",
    install_requires=req,
    setup_requires=req,
    python_requires=">=3",
    extras_require={
        "plot": ["matplotlib", "jupyter", "seaborn"],
    },
    version='4.0.8',
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    name="ISRSpectrum",
    package_data={'ISRSpectrum': ['*.csv']}
)


setup(**config)

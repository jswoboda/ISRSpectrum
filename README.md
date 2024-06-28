# ISRSpectrum

[![GitHub version](https://badge.fury.io/gh/jswoboda%2FISRSpectrum.svg)](https://badge.fury.io/gh/jswoboda%2FISRSpectrum)
[![Coverage Status](https://coveralls.io/repos/jswoboda/ISRSpectrum/badge.svg)](https://coveralls.io/r/jswoboda/ISRSpectrum)
[![Conda image](https://anaconda.org/swoboj/isrspectrum/badges/version.svg)](https://anaconda.org/swoboj/isrspectrum)
[![Documentation Status](https://readthedocs.org/projects/isrspectrum/badge/?version=latest)](https://isrspectrum.readthedocs.io/en/latest/?badge=latest)

by John Swoboda

![alt text](https://raw.github.com/jswoboda/ISRSpectrum/master/logofig.png "ISR Spectrum")

## Overview

This is a Python module to calculate an incoherent scatter spectrum based off of Kudeki and Milla's 2011 IEEE Geophysics paper.

> Kudeki, E.; Milla, M.A., "Incoherent Scatter Spectral Theories—Part I: A General Framework and Results for Small Magnetic Aspect Angles," Geoscience and Remote Sensing, IEEE Transactions on , vol.49, no.1, pp.315,328, Jan. 2011

 Like the model covered in the paper the software can calculate a spectra given, any magnetic aspect angle not perpendicular to B, any number of ion species, and the collision frequencies associated with those ion species. As the magnetic aspect angles get closer to perpendicular to B, usually &lt; 1 degree perp to B, more and more calculations are needed for the Gordeyev to converge.

## Installation

The package can be installed through anaconda through the following command

```sh
conda install -c swoboj isrspectrum
```

The package is also avalible through pypi and using using pip

```sh
pip install isrspectrum
```

The user can also download the software and install it through pip or set up tools. It is suggested the software be installed locally as it can be updated.

```sh
git clone https://github.com/jswoboda/ISRSpectrum.git

cd ISRSpectrum

pip install -e .
```

## Code Examples

In order to create the spectrums, assuming the user installed the code, first import the ISRSpectrum class and create an instance of the class.

```python
import numpy as np
from ISRSpectrum import Specinit
import scipy.constants as spconsts

spfreq = 50e3
nspec = 512
ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=nspec, sampfreq=spfreq,dFlag=True)


ti = 1e3
te = 1e3
Ne = 1e11
Ni = 1e11
mi = 16
Ce = np.sqrt(spconst.Boltzmann*te/spconst.m_e)
Ci = np.sqrt(spconst.Boltzmann*ti/(spconst.m_p*mi))

datablock90 = np.array([[Ni,ti,0,1,mi,0],[Ne,te,0,1,1,0]])
(omega,specorig,rcs) = ISS2.getspec(datablock90, rcsflag = True)
```

Where each colum of the data block is represented as a the information for a ion/electron species. The list is in this order:

1. density(m^-3)
2. temperature(K)
3. Doppler velocity (m/s)
4. charge in elementary charges
5. mass of the species in AMU
6. collision frequency in s^-1.

The info for the electron species must be the last row also !

Alternatively, if you are only using the following species: O+,NO+,N2+,O2+,N+, H+, e-, a simpler interface is available.
Using this interface quasi-neutrality is assumed so the number of positive ions and electrons are the same.
The user can get spectrum in the following way.

```python
from ISRSpectrum import Specinit
import numpy as np

spfreq = 50e3
nspec = 512
ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=nspec, sampfreq=spfreq,dFlag=True)

ti = 1e3
te = 1e3
Ne = 1e11
Ni = 1e11
vi = 0
species = ['O+','e-']
datablock = np.array([[Ni,ti],[Ne,te]])
(omega,specorig,rcs) = ISS2.getspecsep(datablock,species,vi,90,rcsflag = True)
```

Further examples can be found in `Examples/examplespectrums.py`

## Note

I switch the default branch to main. If you have a local copy of the software you can run the following commands to deal with the renaming

```bash
git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```

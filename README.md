[![Build Status](https://travis-ci.org/jswoboda/ISRSpectrum.svg)](https://travis-ci.org/jswoboda/ISRSpectrum)
[![Coverage Status](https://coveralls.io/repos/jswoboda/ISRSpectrum/badge.svg)](https://coveralls.io/r/jswoboda/ISRSpectrum)
##ISRSpectrum
by John Swoboda

![alt text](https://raw.github.com/jswoboda/ISRSpectrum/master/logofig.png "ISR Spectrum")
## Overview
This is a Python module to calculate an ISR spectrum  based off of Kudeki and Milla's 2011 IEEE Geophysics paper.

	Kudeki, E.; Milla, M.A., "Incoherent Scatter Spectral Theories—Part I: A General Framework and Results for Small Magnetic Aspect Angles," Geoscience and Remote Sensing, IEEE Transactions on , vol.49, no.1, pp.315,328, Jan. 2011

The code has been written to be able to produce a spectrum with an arbitrary number of ion species with an arbitrary collsion frequency and magnetic aspect angle. The only issue is that I will now promise that the code will run quickly at magnetic aspect angles <1 degree perp to B.

## Requirements
Standard python package requirements such as scipy and matplotlib are stated specifially in requirements.txt. The package [isrutilities](https://github.com/jswoboda/PythonISRUtilities) also needs to be installed before hand.

## Suggestions
It is highly suggested that the [Anaconda](https://www.continuum.io/downloads) platform be used as the package manager. All of the development and testing has been done using this.
Assuming the user has installed Anaconda a [set up bash script](https://github.com/jswoboda/AnacondaEnvUtilities), which can be used in Linux or Mac environments is avalible.

The user can also take advantage of two different APIs to plot results using the SimISR. The first is in Python and is called [GeoDataPython](https://github.com/jswoboda/GeoDataPython). A MATLAB version of this API is also avalible called [GeoDataMATLAB](https://github.com/jswoboda/GeoDataMATLAB).
## Installation

	 git clone https://github.com/jswoboda/ISRSpectrum.git

	 cd ISRSpectrum
	 python setup.py develop

## Code Examples

In order to create the spectrums, assuming the user installed the code, first import the ISRSpectrum class and create an instanstance of the class.

~~~python
import numpy as np
import ISRSpectrum.ISRSpectrum as ISSnew
from isrutilities.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge, v_me, v_amu

ISS2 = ISSnew.ISRSpectrum(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=129, sampfreq=50e3,dFlag=True)

ti = 1e3
te = 1e3
Ne = 1e11
Ni = 1e11
mi = 16
Ce = np.sqrt(v_Boltz*te/v_me)
Ci = np.sqrt(v_Boltz*ti/(v_amu*mi))

datablock90 = np.array([[Ni,ti,0,1,mi,0],[Ne,te,0,1,1,0]])
(omega,specorig,rcs) = ISS2.getspec(datablock90, rcsflag = True)
~~~

Where each colum of the data block is represented as a the information for a ion/electron species. The list is in this order[density(m^-3),tempreture(K), Doppler velocity (m/s),charge in elementary charges, mass of the species in AMU, collision frequency in s^-1]. The info for the elctron species must be the last row also !

Alternatively, if you are only using the following species: O+,NO+,N2+,O2+,N+, H+, e-, a simpler interface is avalible. Using this interface quasi-neutrality is assumed so the number of positive ions and electrons are the same. The user can get spectrum in the following way.

~~~python
import ISRSpectrum.ISRSpectrum as ISSnew
ISS2 = ISSnew.ISRSpectrum(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=129, sampfreq=50e3,dFlag=True)

ti = 1e3
te = 1e3
Ne = 1e11
Ni = 1e11
vi = 0
species = ['0+','e-']
datablock = np.array([[Ni,ti],[Ne,te]])
(omega,specorig,rcs) = ISS2.getspecsep(datablock,species,vi,90,rcsflag = True)
~~~

Further examples can be found in `Examples/examplespectrums.py`

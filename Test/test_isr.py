#!/usr/bin/env python
import numpy as np
from numpy.testing import assert_allclose
import scipy.constants as spconst
#
from ISRSpectrum import Specinit
ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=129, sampfreq=50e3,dFlag=True)

ti = 1e3
te = 1e3
Ne = 1e11
Ni = 1e11
mi = 16
Ce = np.sqrt(spconst.k*te/spconst.m_e)
Ci = np.sqrt(spconst.k*ti/(spconst.m_p*mi))

datablock90 = np.array([[Ni,ti,0,1,mi,0],[Ne,te,0,1,1,0]])
(omega,specorig,rcs) = ISS2.getspec(datablock90, rcsflag = True)

#%% registration test
"""
Point of this test to simply return to console stderr=0
This is what Travis-CI needs for continuous integration testing
automatically on each git-commit.
"""
#TODO: add omega, specorig, perhaps via HDF5 file loading
# assert_allclose(Ce,123111.44138427741)
# assert_allclose(Ci,720.87235109746211)
assert_allclose(Ne,1e11)
assert_allclose(Ni,1e11)
assert_allclose(datablock90,
                np.array([[1e11,1e3,0,1,16,0],
                          [1e11,1e3,0,1,1, 0]]))
assert_allclose(mi,16)
assert_allclose(rcs,48806562035.874)
assert_allclose(te,1e3)
assert_allclose(ti,1e3)

import numpy as np
from numpy.testing import assert_allclose
#
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

#%% registration test
"""
Point of this test to simply return to console stderr=0
This is what Travis-CI needs for continuous integration testing
automatically on each git-commit.
"""
#TODO: add omega, specorig, perhaps via HDF5 file loading
assert_allclose(Ce,123111.44138427741)
assert_allclose(Ci,720.87235109746211)
assert_allclose(Ne,1e11)
assert_allclose(Ni,1e11)
assert_allclose(datablock90,
                np.array([[1e11,1e3,0,1,16,0],
                          [1e11,1e3,0,1,1, 0]]))
assert_allclose(mi,16)
assert_allclose(rcs,48806554524.200493)
assert_allclose(te,1e3)
assert_allclose(ti,1e3)
assert_allclose(v_Boltz,1.380658e-23)
assert_allclose(v_C_0,299792458)
assert_allclose(v_amu,1.6605402e-27)
assert_allclose(v_elemcharge,1.602177e-19)
assert_allclose(v_epsilon0,8.854188e-12)
assert_allclose(v_me,9.1093897e-31)

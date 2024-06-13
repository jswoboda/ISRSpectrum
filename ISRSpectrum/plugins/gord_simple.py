#!python

import numpy as np
import scipy.special as sp
import scipy.constants as spconst

import scipy.special as sp_spec


class GordPlug():

    def calcgordeyev(self,dataline, K, omeg,dFlag):
        (Ns, Ts, Vs, qs, ms, nus) = dataline[:6]

        C = np.sqrt(spconst.k * Ts / ms)
        omeg_s = omeg - K * Vs
        theta = omeg_s / (K * C * np.sqrt(2.0))

        num_g = np.sqrt(np.pi) * np.exp(-(theta**2)) - 1j * 2.0 * sp_spec.dawsn(theta)
        den_g = K * C * np.sqrt(2)
        gord = num_g / den_g
        if dFlag:
            print("\t No collisions No magnetic field,again")
        return (gord, Ts, Ns, omeg_s)
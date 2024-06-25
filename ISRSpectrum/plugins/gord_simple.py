#!python

import numpy as np
import scipy.special as sp
import scipy.constants as spconst
import scipy.special as sp_spec


class GordPlug:

    def calcgordeyev(self, dataline, K, omeg, dFlag, *args, **kwargs):
        """
        Direct calculation of Gordeyev integral with no magnetic field or collisions.

        Parameters
        ----------
        dataline : ndarray
            A numpy array of length that holds the plasma parameters needed to create the spectrum. Each row of the array will have the following set up.
                [Ns, Ts, Vs, qs, ms, nus]
                Ns - The density of the species in m^-3
                Ts - Temperature of the species in degrees K
                Vs - The Doppler velocity in m/s.
                qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
                ms - Mass of the species in AMU. (Value will be replaced for the electrons)
                nus - Collision frequency for species in s^-1.
        alphadeg : float
            The magnetic aspect angle in radians.
        K : float
            K value of the radar in rad/m
        omeg : ndarray
            Radian frequency vector the integral will be evaluated over.
        dFlag : float
            A debug flag, if set true will output debug text. Default is false.

        Returns
        -------
        gord : ndarray
            The result of the Gordeyev integral over Doppler corrected radian frequency
        hs : ndarray
            The Debye length in m.
        Ns : ndarray
            The density of the species in m^-3
        omeg_s : ndarray
            An array of the Doppler corrected radian frequency

        """
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

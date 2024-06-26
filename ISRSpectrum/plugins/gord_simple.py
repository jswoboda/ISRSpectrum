#!python

"""
This is an example of a plugin that performs the basic dawsons integral which is a call to a scipy function. This is simply an example of a plugin that a user can follow as an example.
"""

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
                [Ns, Ts, Vs, ms, nus]
                Ns - The density of the species in m^-3
                Ts - Temperature of the species in degrees K
                Vs - The Doppler velocity in m/s.
                qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
                ms - Mass of the species in kg.
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
        qs : float
            The charge of the species in elementary charges.
        omeg_s : ndarray
            An array of the Doppler corrected radian frequency

        """
        assert (
            len(dataline) >= 5
        ), "The dataline input needs to be length of at least 4 elements."
        (Ns, Ts, Vs, qs, ms) = dataline[:5]

        assert (
            ms > 0 and ms < 1e-20
        ), "Atomic masses should be very small and greater than 0."

        C = np.sqrt(spconst.k * Ts / ms)
        omeg_s = omeg - K * Vs
        theta = omeg_s / (K * C * np.sqrt(2.0))

        num_g = np.sqrt(np.pi) * np.exp(-(theta**2)) - 1j * 2.0 * sp_spec.dawsn(theta)
        den_g = K * C * np.sqrt(2)
        gord = num_g / den_g

        if dFlag:
            print("\t No collisions No magnetic field,again")
        return (gord, Ts, Ns, qs, omeg_s)

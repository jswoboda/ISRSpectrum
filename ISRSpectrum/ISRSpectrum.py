#!/usr/bin/env python
"""
Created on Tue Jul 15 16:12:05 2014

@author: John Swoboda
ISRSpectrum.py
This module will create ISR spectrums using plasma parameters. The code implements the
method of calculation found in the following article.

E. Kudeki and M. A. Milla, 2011.

The intent of the code is to be able to calculate an ISR spectrum in a number of
different conditions except for a very low magnetic aspect angles (<1deg).
"""
from __future__ import absolute_import
from pathlib import Path
from six import string_types
import numpy as np
import scipy as sp
import scipy.special as sp_spec
import pandas as pd
import scipy.constants as spconst
from .mathutils import sommerfelderfrep

INFODICT = {
    "O+": np.array([1, 16]),
    "NO+": np.array([1, 30]),
    "N2+": np.array([1, 28]),
    "O2+": np.array([1, 32]),
    "N+": np.array([1, 14]),
    "H+": np.array([1, 1]),
    "HE+": np.array([1, 2]),
    "e-": np.array([-1, 1]),
}


class Specinit(object):
    """ Class to create the spectrum.

    The instance of the class will hold infomation on the radar system such as
    sample frequency, center frequency and number of points for the spectrum.

    Attributes
    ----------
    bMag : float
        The magnetic field magnitude in Teslas.
    K : float
        The Bragg scattering vector magnitude corresponds to 1/2 radar wavelength.
    f : array_like
        Vector holding the frequency values in Hz.
    omeg : array_like
        Vector hold the frequency values in rad/s
    collfreqmin : float
        Minimum collision frequency before they are taken into account in the Gordeyev integral calculation
    alphamax : float
        The maximum aspect angle
    dFlag : bool
        A debug flag."""

    def __init__(
        self,
        centerFrequency=440.2 * 1e6,
        bMag=0.4e-4,
        nspec=64,
        sampfreq=50e3,
        collfreqmin=1e-2,
        alphamax=30.0,
        dFlag=False,
        f=None,
    ):
        """ Constructor for the class.
        Parameters
        ----------
        centerFrequency : float
            The radar center frequency in Hz.
        bMag : float
            The magnetic field magnitude in Teslas.
        nspec : int
            The number of points of the spectrum.
        sampfreq : float
            The sampling frequency of the A/Ds in Hz
        collfreqmin : float
            (Default 1e-2) The minimum collision frequency needed to incorporate it into Gordeyev
            integral calculations in units of K*sqrt(Kb*Ts/ms) for each ion species.
        alphamax : float
            (Default 30)  The maximum magnetic aspect angle in which the B-field will be taken into account.
        dFlag : float
            A debug flag, if set true will output debug text. Default is false.
        f : array_like
            Array of frequeny points, in Hz, the spectrum will be formed over. Default is
            None, at that point the frequency vector will be formed using the number of points for the spectrum
            and the sampling frequency to create a linearly sampled frequency vector. """
        self.bMag = bMag
        self.dFlag = dFlag
        self.collfreqmin = collfreqmin
        self.alphamax = alphamax
        curpath = Path(__file__).parent
        self.Bst = pd.read_csv(str(curpath / "ion2ion.csv"), index_col=0)

        self.Cin = pd.read_csv(str(curpath / "ion2neu.csv"), index_col=0)

        self.K = (
            2.0 * np.pi * 2 * centerFrequency / spconst.c
        )  # The Bragg scattering vector, corresponds to half the radar wavelength.
        if f is None:
            minfreq = -np.ceil((nspec - 1.0) / 2.0)
            maxfreq = np.floor((nspec - 1.0) / 2.0 + 1)
            self.f = np.arange(minfreq, maxfreq) * (
                sampfreq / (2 * np.ceil((nspec - 1.0) / 2.0))
            )
        else:
            self.f = f
        self.omeg = 2.0 * np.pi * self.f

    def getspec(self, datablock, alphadeg=90.0, rcsflag=False, seplines=False):
        """ Gives the spectrum and the frequency vectors.

        Get spectrum array given the block of data and the magnetic aspect angle.

        Parameters
        ----------
        datablock : array_like
            A numpy array of size 1+Nionsx6 that holds the plasma parameters needed
            to create the spectrum. The last row will hold the information for the electrons.
            Each row of the array will have the following set up.
                [Ns, Ts, Vs, qs, ms, nus]
                Ns - The density of the species in m^-3
                Ts - Tempretur of the species in degrees K
                Vs - The Doppler velocity in m/s.
                qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
                ms - Mass of the species in AMU. (Value will be replaced for the electrons)
                nus - Collision frequency for species in s^-1.
        alphadeg : float
            The magnetic aspect angle in degrees.
        rcsflag : float
            A bool that will determine if the reflected power is returned as well. (default is False)
        seplines :bool
            A bool that will change the output, spec to a list of numpy arrays that include
            return from the electron line and ion lines seperatly.

        Returns
        -------
        f : array_like
            A numpy array that holds the frequency vector associated with the spectrum, in Hz.
        spec : array_like
            The resulting ISR spectrum from the parameters with the max value set to 1.
        rcs : array_like
            The RCS from the parcle of plasma for the given parameters. The RCS is in m^2.
        """
        # perform a copy of the object to avoid it being written over with incorrect.
        datablock = datablock.copy()
        dFlag = self.dFlag
        alpha = alphadeg * np.pi / 180
        estuff = datablock[-1]
        Nions = datablock.shape[0] - 1
        estuff[3] = -spconst.e
        estuff[4] = spconst.m_e
        ionstuff = datablock[:-1]
        if dFlag:
            print("Calculating Gordeyev int for electons")
        (egord, Te, Ne, omeg_e) = self.__calcgordeyev__(estuff, alpha)
        h_e = np.sqrt(spconst.epsilon_0 * spconst.k * Te / (Ne * spconst.e * spconst.e))
        sig_e = (1j + omeg_e * egord) / (self.K ** 2 * h_e ** 2)
        nte = 2 * Ne * np.real(egord)

        # adjust ion stuff
        ionstuff[:, 3] = ionstuff[:, 3] * spconst.e
        ionstuff[:, 4] = ionstuff[:, 4] * spconst.m_p
        ionden = np.sum(ionstuff[:, 0])
        # normalize total ion density to be the same as electron density
        ionstuff[:, 0] = (estuff[0] / ionden) * ionstuff[:, 0]

        # ratio of charges between ion species and electrons
        qrotvec = ionstuff[:, 3] / estuff[3]
        firstion = True
        wevec = np.zeros(Nions)
        Tivec = np.zeros(Nions)
        for iion, iinfo in enumerate(ionstuff):
            if dFlag:
                print("Calculating Gordeyev int for ion species #{:d}".format(iion))

            (igord, Ti, Ni, omeg_i) = self.__calcgordeyev__(iinfo, alpha)

            wevec[iion] = Ni / Ne
            Tivec[iion] = Ti
            # sub out ion debye length because zero density of ion species
            # can cause a divid by zero error.
            # tempreture ratio
            mu = Ti / Te
            qrot = qrotvec[iion]
            sig_i = (
                (Ni / Ne)
                * (1j + omeg_i * igord)
                / (self.K ** 2 * mu * h_e ** 2 / qrot ** 2)
            )
            nti = 2 * Ni * np.real(igord)

            if firstion:
                sig_sum = sig_i
                nt_sum = nti
                firstion = False
            else:
                sig_sum = sig_i + sig_sum
                nt_sum = nti + nt_sum

        inum = np.absolute(sig_e) ** 2 * nt_sum
        enum = np.absolute(1j + sig_sum) ** 2 * nte
        den = np.absolute(1j + sig_sum + sig_e) ** 2
        iline = inum / den
        eline = enum / den
        if seplines:
            spec = [iline, eline]
        else:
            spec = iline + eline

        if rcsflag:
            Tr = Te / np.sum(wevec * Tivec)
            rcs = Ne / (
                (1 + self.K ** 2 * h_e ** 2) * (1 + self.K ** 2 * h_e ** 2 + Tr)
            )

            return (self.f, spec, rcs)
        else:
            return (self.f, spec)

    def getspecsep(
        self,
        datablock,
        species,
        vel=0.0,
        alphadeg=90.0,
        rcsflag=False,
        col_calc=False,
        n_datablock=None,
        n_species=None,
        seplines=False,
    ):
        """ This function is a different way of getting the spectrums.

        A datablock is still used, (Nsp x 2) numpy array, but it is filled in by using
        the species listed as a string in the list speces.

        Parameters
        ----------
        datablock : array_like
            A numpy array of size Nsp x2. The first element of the row is the
            density of species in m^-3 and the second element is the Tempreture in
            degrees K.
        species : array_like
            A Nsp list of strings that label each species.
        vel : array_like
            The line of site velocity of the ions in m/s, default = 0.0
        alphadeg : float
            The angle offset from the magnetic field in degrees, default is 90.0
        rcsflag : bool
            If this flag is True then a third output of the rcs will be made. Default value is False
        col_calc : bool
            If this flag is true then collisions will be calculated. (Default= False)
        n_datablock : array_like
            A numpy array of size Nnsp x2. The first element of the row is
            the density of the neutral species in m^-3 and the second element is the
            Tempreture in degrees K.  If set to None then (Default=None)
        n_species : array_like
            A Nsp list of strings that label each neutral species.(Default=None)
        seplines : bool
            A bool that will change the output, spec to a list of numpy arrays that include
            return from the electron line and ion lines seperatly.

        Returns
        -------
        f : array_like
            A numpy array that holds the frequency vector associated with the spectrum, in Hz.
        spec : array_like
            The resulting ISR spectrum from the parameters with the max value set to 1.
        rcs : array_like
            The RCS from the parcle of plasma for the given parameters. The RCS is in m^2.
        """
        assert Islistofstr(species), "Species needs to be a list of strings"
        assert allin(
            species, list(INFODICT.keys())
        ), "Have un named species in the list."
        nspec = datablock.shape[0]
        datablocknew = np.zeros((nspec, 7))

        (nuparr, nuperp) = get_collisionfreqs(
            datablock, species, self.Bst, self.Cin, n_datablock, n_species
        )
        for nspec, ispec in enumerate(species):
            datablocknew[nspec, :2] = datablock[nspec]
            datablocknew[nspec, 2] = vel
            datablocknew[nspec, 3:5] = INFODICT[ispec]
            if col_calc:
                datablocknew[nspec, 5] = nuparr[nspec]
                datablocknew[nspec, 6] = nuparr[nspec]

        return self.getspec(datablocknew, alphadeg, rcsflag, seplines=seplines)

    def __calcgordeyev__(self, dataline, alpha, alphdiff=10.0):
        """ Performs the Gordeyve integral calculation.

        Parameters
        ----------
        dataline : array_like
            A numpy array of length that holds the plasma parameters needed
            to create the spectrum.
            Each row of the array will have the following set up.
                [Ns, Ts, Vs, qs, ms, nus]
                Ns - The density of the species in m^-3
                Ts - Tempretur of the species in degrees K
                Vs - The Doppler velocity in m/s.
                qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
                ms - Mass of the species in AMU. (Value will be replaced for the electrons)
                nus - Collision frequency for species in s^-1.
        alphadeg : float
            The magnetic aspect angle in radians.

        Returns
        -------
        gord : array_like
            The result of the Gordeyev integral over Doppler corrected radian frequency
        hs : array_like
            The Debye length in m.
        Ns : array_like
            The density of the species in m^-3
        omeg_s : array_like
            An array of the Doppler corrected radian frequency
        """
        dFlag = self.dFlag
        K = self.K
        (Ns, Ts, Vs, qs, ms, nus) = dataline[:6]
        if len(dataline) == 7:
            nuperp = dataline[-1]
        else:
            nuperp = nus

        C = np.sqrt(spconst.k * Ts / ms)
        omeg_s = self.omeg - K * Vs
        theta = omeg_s / (K * C * np.sqrt(2.0))
        Om = qs * self.bMag / ms
        # determine what integral is used
        magbool = alpha * 180.0 / np.pi < self.alphamax
        collbool = self.collfreqmin * K * C < nus

        if not collbool and not magbool:
            # for case with no collisions or magnetic field just use analytic method
            num_g = np.sqrt(np.pi) * np.exp(-(theta ** 2)) - 1j * 2.0 * sp_spec.dawsn(
                theta
            )
            den_g = K * C * np.sqrt(2)
            gord = num_g / den_g
            if dFlag:
                print("\t No collisions No magnetic field")
            return (gord, Ts, Ns, omeg_s)

        elif collbool and not magbool:
            if dFlag:
                print("\t With collisions No magnetic field")
            gordfunc = collacf
            exparams = (K, C, nus)
        elif not collbool and magbool:
            if dFlag:
                print("\t No collisions with magnetic field")
            gordfunc = magacf
            exparams = (K, C, alpha, Om)
        else:
            if dFlag:
                print("\t With collisions with magnetic field")
            gordfunc = magncollacf
            exparams = (K, C, alpha, Om, nus)

        maxf = np.abs(self.f).max()
        T_s = 1.0 / (2.0 * maxf)

        #        N_somm = 2**15
        #        b1 = 10.0/(K*C*np.sqrt(2.0))
        # changed ot sample grid better
        N_somm = 2 ** 10
        b1 = T_s * N_somm
        #        b1 = interval/10.
        #        N_somm=np.minimum(2**10,np.ceil(b1/T_s))
        (gord, flag_c, outrep) = sommerfelderfrep(
            gordfunc, N_somm, omeg_s, b1, Lmax=500, errF=1e-7, exparams=exparams
        )
        if dFlag:
            yna = ["No", "Yes"]
            print(
                "\t Converged: {:s}, Number of iterations: {:d}".format(
                    yna[flag_c], outrep
                )
            )

        return (gord, Ts, Ns, omeg_s)


def magacf(tau, K, C, alpha, Om):
    """ Create a single particle acf for a species with magnetic field but no collisions.

    Parameters
    ----------
    tau : array_like
        The time vector for the acf.
    K : float
        Bragg scatter vector magnetude.
    C : float
        Thermal speed of the species.
    alpha : float
        Magnetic aspect angle in radians.
    Om : float
        The gyrofrequency of the particle.

    Returns
    -------
    acf : array_like
        The single particle acf.
    """
    Kpar = np.sin(alpha) * K
    Kperp = np.cos(alpha) * K
    return np.exp(
        -np.power(C * Kpar * tau, 2.0) / 2.0
        - 2.0 * np.power(Kperp * C * np.sin(Om * tau / 2.0) / Om, 2.0)
    )


def collacf(tau, K, C, nu):
    """ Create a single particle acf for a species with collisions, no magnetic field.

    Parameters
    ----------
    tau : array_like
        The time vector for the acf.
    K : float
        Bragg scatter vector magnetude.
    C : float
        Thermal speed of the species.
    nu : float
        The collision frequency in collisions/sec

    Returns
    -------
    acf : array_like
        The single particle acf.
    """
    return np.exp(-np.power(K * C / nu, 2.0) * (nu * tau - 1 + np.exp(-nu * tau)))


def magncollacf(tau, K, C, alpha, Om, nu):
    """ Create a single particle acf for a species with magnetic field and collisions.

    Parameters
    ----------
    tau : array_like
        The time vector for the acf.
    K : float
        Bragg scatter vector magnetude.
    C : float
        Thermal speed of the species.
    alpha : float
        Magnetic aspect angle in radians.
    Om : float
        The gyrofrequency of the particle.
    nu : float
        The collision frequency in collisions/sec

    Returns
    -------
    acf : array_like
        The single particle acf.
    """
    Kpar = np.sin(alpha) * K
    Kperp = np.cos(alpha) * K
    gam = np.arctan(nu / Om)

    deltl = np.exp(-np.power(Kpar * C / nu, 2.0) * (nu * tau - 1 + np.exp(-nu * tau)))
    deltp = np.exp(
        -np.power(C * Kperp, 2.0)
        / (Om * Om + nu * nu)
        * (
            np.cos(2 * gam)
            + nu * tau
            - np.exp(-nu * tau) * (np.cos(Om * tau - 2.0 * gam))
        )
    )
    return deltl * deltp


def Islistofstr(inlist):
    """ This function will determine if the input is a list of strings

    Parameters
    ----------
    inlist: list
        List to check for strings.

    Returns
    -------
    isstrs : bool
        If all the list is strings.
    """
    if not isinstance(inlist, list):
        return False
    for item in inlist:
        if not isinstance(item, string_types):
            return False
    return True


def allin(inp, reflist):
    """ This function will determine if all of the strings in the list input are in
    the list reflist.


    """
    for item in inp:
        if item not in reflist:
            return False
    return True


#%% Collision frequencies


def get_collisionfreqs(datablock, species, Bst, Cin, n_datablock=None, n_species=None):
    """Calculate collision frequencies
    Uses the methods shown in Schunk and Nagy (2009) chapter 4.

    Parameters
    ----------
    datablock : array_like
        A numpy array of size Nsp x2. The first element of the row is the
        density of species in m^-3 and the second element is the Tempreture in
        degrees K.
    Bst : array_like
        Ion ion collision constants read in from a file.
    Cin : array_like
        Ion neutral collision constants read in from a file.
    species : list
        A Nsp list of strings that label each species.
    col_calc : bool
        If this flag is true then collisions will be calculated. (Default= False)
    n_datablock : array_like
        A numpy array of size Nnsp x2. The first element of the row is
            the density of the neutral species in m^-3 and the second element is the
            Tempreture in degrees K.If set to None then (Default=None)
    n_species : array_like
        A Nnsp list of strings that label each neutral species.(Default=None)

    Returns
    -------
    nuparr : array_like
        A Nsp length numpy array that holds the parrallel collision frequencies in s^-1.
    nuperp : array_like
        A Nsp length numpy array that holds the perpendictular collision frequencies in s^-1.
    """

    nuperp = np.zeros((datablock.shape[0]))
    nuparr = np.zeros_like(nuperp)

    Ni = datablock[:-1, 0] * 1e-6
    Ti = datablock[:-1, 1]
    Ne = datablock[-1, 0] * 1e-6
    Te = datablock[-1, 1]

    # electron electron and electron ion collisions
    # Schunk and Nagy eq 4.144 and eq 4.145
    nuee = 54.5 / np.sqrt(2) * Ne / np.power(Te, 1.5)
    nuei = 54.5 * Ni / np.power(Te, 1.5)
    nuei = nuei.sum()

    nuperp[-1] = nuei + nuee
    nuparr[-1] = nuei
    # ionion collisions
    for si, s in enumerate(species[:-1]):
        for ti, t in enumerate(species[:-1]):

            nuparr[si] = nuparr[si] + Bst[s][t] * Ni[ti] / np.power(Ti[ti], 1.5)

    if (n_datablock is not None) and (n_species is not None):
        Nn = n_datablock[:, 0] * 1e-6
        Tn = n_datablock[:, 1]
        # ion neutral collisions
        for si, s in enumerate(species[:-1]):
            for ti, t in enumerate(n_species):

                curCin = Cin[s][t]

                if np.isnan(curCin):
                    neuts = r_ion_neutral(s, t, Ni[si], Nn[ti], Ti[si], Tn[ti])
                else:
                    neuts = curCin * Ni[si]
                nuparr[si] = nuparr[si] + neuts
        # electron neutral collision frequencies
        for ti, t in enumerate(n_species):
            nueneu = e_neutral(t, Nn[ti], Te)
            nuperp[-1] = nuperp[-1] + nueneu
            nuparr[-1] = nuparr[-1] + nueneu

    # according to milla and kudeki the collision for perpendicular and parrallel
    # directions only matter for the electrons
    nuperp[:-1] = nuparr[:-1]
    return (nuparr, nuperp)


def r_ion_neutral(s, t, Ni, Nn, Ti, Tn):
    """ This will calculate resonant ion - neutral reactions collision frequencies. See
    table 4.5 in Schunk and Nagy.
    Inputs
    s - Ion name string
    t - neutral name string
    Ni - Ion density cm^-3
    Nn - Neutral density cm^-3
    Ti - Ion tempreture K
    Tn - Neutral tempreture K
    Outputs
    nu_ineu - collision frequency s^-1
    """
    Tr = (Ti + Tn) * 0.5
    sp1 = (s, t)
    # from Schunk and Nagy table 4.5
    nudict = {
        ("H+", "H"): [2.65e-10, 0.083],
        ("He+", "He"): [8.73e-11, 0.093],
        ("N+", "N"): [3.84e-11, 0.063],
        ("O+", "O"): [3.67e-11, 0.064],
        ("N2+", "N"): [5.14e-11, 0.069],
        ("O2+", "O2"): [2.59e-11, 0.073],
        ("H+", "O"): [6.61e-11, 0.047],
        ("O+", "H"): [4.63e-12, 0.0],
        ("CO+", "CO"): [3.42e-11, 0.085],
        ("CO2+", "CO"): [2.85e-11, 0.083],
    }
    A = nudict[sp1][0]
    B = nudict[sp1][1]
    if sp1 == ("O+", "H"):

        nu_ineu = A * Nn * np.power(Ti / 16.0 + Tn, 0.5)
    elif sp1 == ("H+", "O"):
        nu_ineu = A * Nn * np.power(Ti, 0.5) * (1 - B * np.log10(Ti)) ** 2
    else:
        nu_ineu = A * Nn * np.power(Tr, 0.5) * (1 - B * np.log10(Tr)) ** 2
    return nu_ineu


def e_neutral(t, Nn, Te):
    """ This will calculate electron - neutral reactions collision frequencies. See
    table 4.6 in Schunk and Nagy.

    Inputs
    t - neutral name string
    Nn - Neutral density cm^-3
    Te - electron tempreture K
    Outputs
    nu_ineu - collision frequency s^-1"""
    if t == "N2":
        return 2.33e-11 * Nn * (1 - 1.21e-4 * Te) * Te
    elif t == "O2":
        return 1.82e-10 * Nn * (1 - 3.6e-2 * Te ** 0.5) * Te ** 0.5
    elif t == "O":
        return 8.9e-11 * Nn * (1 - 5.7e-4 * Te) * Te ** 0.5
    elif t == "He":
        return 4.6e-10 * Nn * Te ** 0.5
    elif t == "H":
        return 4.5e-9 * Nn * (1 - 1.35e-4 * Te) * Te ** 0.5
    elif t == "CO":
        return 2.34e-11 * Nn * (165 + Te)
    elif t == "CO2":
        return 3.68e-8 * Nn * (1 - 4.1e-11 * np.abs(4500.0 - Te) ** 2.93)
    else:
        return 0.0

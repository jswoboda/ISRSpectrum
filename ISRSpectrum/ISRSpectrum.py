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
from six import string_types
import numpy as np
import scipy.constants as spconst
from .plugins import gordplugs
from .collision_calc import get_collisionfreqs, getionmass, getionmass, INFODICT


class Specinit(object):
    """Class to create the spectrum.

    The instance of the class will hold infomation on the radar system such as
    sample frequency, center frequency and number of points for the spectrum.

    Attributes
    ----------
    bMag : float
        The magnetic field magnitude in Teslas.
    K : float
        The Bragg scattering vector magnitude corresponds to 1/2 radar wavelength.
    f : ndarray
        Vector holding the frequency values in Hz.
    omeg : ndarray
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
        """Constructor for the class.

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
            (Default 1e-2) The minimum collision frequency needed to incorporate it into Gordeyev integral calculations in units of K*sqrt(Kb*Ts/ms) for each ion species.
        alphamax : float
            (Default 30)  The maximum magnetic aspect angle in which the B-field will be taken into account.
        dFlag : float
            A debug flag, if set true will output debug text. Default is false.
        f : ndarray
            Array of frequeny points, in Hz, the spectrum will be formed over. Default is None, at that point the frequency vector will be formed using the number of points for the spectrum and the sampling frequency to create a linearly sampled frequency vector.
        plugin_folder : str
            A folder for possible plugins.
        """

        self._plugins = gordplugs

        self.bMag = bMag
        self.dFlag = dFlag
        self.collfreqmin = collfreqmin
        self.alphamax = alphamax

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

    def getspec(
        self,
        datablock,
        des_plug="default",
        alphadeg=90.0,
        rcsflag=False,
        seplines=False,
        *args,
        **kwargs
    ):
        """Gives the spectrum and the frequency vectors.

        Get spectrum array given the block of data and the magnetic aspect angle.

        Parameters
        ----------
        datablock : ndarray
            A numpy array of size 1+Nionsxa which is the size  that holds the plasma parameters needed to create the spectrum. The last row will hold the information for the electrons. The first two entries in each row must be the density and temperature.
        des_plug : str
            The desired pluggin for the Gordeyev integral
        alphadeg : float
            The magnetic aspect angle in degrees.
        rcsflag : float
            A bool that will determine if the reflected power is returned as well. (default is False)
        seplines :bool
            A bool that will change the output, spec to a list of numpy arrays that include return from the electron line and ion lines seperatly.
        args : list
            Any extra inputs needed for the Gordeyev integrals.
        kwargs : dict
            Extra inputs for the Gordeyev integrals.

        Returns
        -------
        f : ndarray
            A numpy array that holds the frequency vector associated with the spectrum, in Hz.
        spec : ndarray
            The resulting ISR spectrum from the parameters with the max value set to 1.
        rcs : ndarray
            The RCS from the parcle of plasma for the given parameters. The RCS is in m^2.
        """
        # perform a copy of the object to avoid it being written over with incorrect.
        datablock = datablock.copy()
        dFlag = self.dFlag
        alpha = alphadeg * np.pi / 180
        estuff = datablock[-1]
        Nions = datablock.shape[0] - 1

        ionstuff = datablock[:-1]
        if dFlag:
            print("Calculating Gordeyev int for electons")
        plug = self._plugins[des_plug]
        (egord, Te, Ne, qe, omeg_e) = plug.calcgordeyev(
            estuff,
            alpha=alpha,
            K=self.K,
            omeg=self.omeg,
            bMag=self.bMag,
            collfreqmin=self.collfreqmin,
            alphamax=self.alphamax,
            dFlag=self.dFlag,
        )
        h_e = np.sqrt(spconst.epsilon_0 * spconst.k * Te / (Ne * spconst.e * spconst.e))
        sig_e = (1j + omeg_e * egord) / (self.K**2 * h_e**2)
        nte = 2 * Ne * np.real(egord)

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

            (igord, Ti, Ni, qi, omeg_i) = plug.calcgordeyev(
                iinfo,
                alpha=alpha,
                K=self.K,
                omeg=self.omeg,
                bMag=self.bMag,
                collfreqmin=self.collfreqmin,
                alphamax=self.alphamax,
                dFlag=self.dFlag,
            )
            wevec[iion] = Ni / Ne
            Tivec[iion] = Ti
            # sub out ion debye length because zero density of ion species
            # can cause a divid by zero error.
            # temperature ratio
            mu = Ti / Te
            qrot = np.abs(qi / qe)
            sig_i = (
                (Ni / Ne) * (1j + omeg_i * igord) / (self.K**2 * mu * h_e**2 / qrot**2)
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
            rcs = Ne / ((1 + self.K**2 * h_e**2) * (1 + self.K**2 * h_e**2 + Tr))

            return (self.f, spec, rcs)
        else:
            return (self.f, spec)

    def getspecsimple(
        self,
        Ne,
        Te,
        Ti,
        ionspecies,
        ionfracs,
        vel=0.0,
        alphadeg=90.0,
        rcsflag=False,
        seplines=False,
    ):
        """Creates a spectrum using electron physical parameters not in a data block.

        Assuming only the following ionspecies 'O+', 'NO+', 'N2+', 'O2+', 'N+', 'H+', 'He+'

        Parameters
        ----------
        Ne : float
            electron density at the given altitude, m^-3
        Te : float
            electron temperature, K
        Ti : float
            ion temperature, K
        ionspecies : list
            Names of ionspecies used. Includes 'O+', 'NO+', 'N2+', 'O2+', 'N+', 'H+', 'He+'.
        ionfracs : list
            Fractions of each ionspecies. Will be normalized to sum to one.
        vel : ndarray
            The line of site velocity of the ions in m/s, default = 0.0
        alphadeg : float
            The angle offset from the magnetic field in degrees, default is 90.0
        rcsflag : bool
            If this flag is True then a third output of the rcs will be made. Default value is False
        seplines : bool
            A bool that will change the output, spec to a list of numpy arrays that include return from the electron line and ion lines seperatly.

        Returns
        -------
        f : ndarray
            A numpy array that holds the frequency vector associated with the spectrum, in Hz.
        spec : ndarray
            The resulting ISR spectrum from the parameters with the max value set to 1.
        rcs : ndarray
            The RCS from the parcle of plasma for the given parameters. The RCS is in m^2.
        """

        assert Islistofstr(ionspecies), "Species needs to be a list of strings"
        assert allin(
            ionspecies, list(INFODICT.keys())
        ), "Have unnamed species in the list."
        assert len(ionspecies) == len(
            ionfracs
        ), "Length of species names and fractions need to be same length"

        nspec = len(ionspecies) + 1
        datablock = np.zeros((nspec, 2))
        ionfracs = np.array(ionfracs) / sum(ionfracs)

        for inum, ifrac in enumerate(ionfracs):
            datablock[inum, 0] = ifrac * Ne
            datablock[inum, 1] = Ti

        datablock[-1, 0] = Ne
        datablock[-1, 1] = Te
        species = ionspecies + ["e-"]
        return self.getspecsep(datablock, species, vel, alphadeg, rcsflag, seplines)

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
        """This function is a different way of getting the spectrums.

        A datablock is still used, (Nsp x 2) numpy array, but it is filled in by using the species listed as a string in the list speces. This function will also get the collision freuqencies for calculating the spectrum. The colision frequency calculations by default use the CSV files included with the package.

        Parameters
        ----------
        datablock : ndarray
            A numpy array of size Nsp x2. The first element of the row is the density of species in m^-3 and the second element is the temperature in degrees K.
        species : ndarray
            A Nsp list of strings that label each species.
        vel : ndarray
            The line of site velocity of the ions in m/s, default = 0.0
        alphadeg : float
            The angle offset from the magnetic field in degrees, default is 90.0
        rcsflag : bool
            If this flag is True then a third output of the rcs will be made. Default value is False
        col_calc : bool
            If this flag is true then collisions will be calculated. (Default= False)
        n_datablock : ndarray
            A numpy array of size Nnsp x2. The first element of the row is the density of the neutral species in m^-3 and the second element is the temperature in degrees K.  If set to None then (Default=None)
        n_species : ndarray
            A Nsp list of strings that label each neutral species.(Default=None)
        seplines : bool
            A bool that will change the output, spec to a list of numpy arrays that include return from the electron line and ion lines seperatly.

        Returns
        -------
        f : ndarray
            A numpy array that holds the frequency vector associated with the spectrum, in Hz.
        spec : ndarray
            The resulting ISR spectrum from the parameters with the max value set to 1.
        rcs : ndarray
            The RCS from the parcle of plasma for the given parameters. The RCS is in m^2.
        """
        assert Islistofstr(species), "Species needs to be a list of strings"
        assert allin(
            species, list(INFODICT.keys())
        ), "Have unnamed species in the list."
        nspec = datablock.shape[0]
        datablocknew = np.zeros((nspec, 7))

        (nuparr, nuperp) = get_collisionfreqs(
            datablock, species, n_datablock=n_datablock, n_species=n_species
        )
        for nspec, ispec in enumerate(species):
            datablocknew[nspec, :2] = datablock[nspec]
            datablocknew[nspec, 2] = vel
            datablocknew[nspec, 3:5] = INFODICT[ispec]
            if col_calc:
                datablocknew[nspec, 5] = nuparr[nspec]
                datablocknew[nspec, 6] = nuparr[nspec]

        return self.getspec(
            datablocknew, alphadeg=alphadeg, rcsflag=rcsflag, seplines=seplines
        )


def Islistofstr(inlist):
    """This function will determine if the input is a list of strings

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
    """This function will determine if all of the strings in the list input are in the list reflist.

    Parameters
    ----------
    inp : list
        Input list of strings.
    reflist : list
        List of strings that will be compared.
    """
    for item in inp:
        if item not in reflist:
            return False
    return True

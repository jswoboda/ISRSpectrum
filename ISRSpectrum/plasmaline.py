#!/usr/bin/env python
"""plasmaline.py

This will hold material to create plasma line power spectral densities.


"""
import numpy as np
import scipy.signal as sig
import scipy.constants as sconst


class PLspecinit(object):
    """Class to create the spectrum.

    The instance of the class will hold infomation on the radar system such as
    sample frequency, center frequency and number of points for the spectrum.

    Attributes
    ----------
    bMag : float
        The magnetic field magnitude in Teslas.
    K : float
        The Bragg scattering vector magnitude corresponds to 1/2 radar wavelength.
    freq_vec : ndarray
        Center frequency of each channel.
    freq_ind : ndarray
        Index for each frequency channel.
    freq_l : ndarray
        Lower bound frequency of each channel.
    freq_h : ndarray
        Upper bound frequency of each channel
    cfreqvec : ndarray
        The sub-band frequency domain in Hz.
    alphamax : float
        The maximum aspect angle
    dFlag : bool
        A debug flag."""

    def __init__(
        self,
        centerFrequency=440.2 * 1e6,
        bMag=0.4e-4,
        dFlag=False,
        fs=25e6,
        nchans=250,
        nfreq_pfb=1024,
    ):
        """Initiates the container class for the plasma line spectrum functions.

        Parameters
        ----------
        centerFrequency : float
            The radar center frequency in Hz.
        bMag : float
            The magnetic field magnitude in Teslas.
        dFlag : bool
            A debug flag, if set true will output debug text. Default is false.
        fs : float
            Full sampling freuqency of the whole plasma channel.
        nchans : int
            Number of frequency channels used to represent the plasma line.
        nfreq_pfb : int
            Number of frequency samples for each frequency channel that will be used to create the plasma line spectrum.
        """
        self.K = (
            2.0 * np.pi * 2 * centerFrequency / sconst.c
        )  # The Bragg scattering vector, corresponds to half the radar wavelength.
        self.bMag = bMag
        nchans = 250
        c_bw = fs / nchans
        n_fpfb = 1024
        c_bw = fs / nchans

        self.freq_vec = np.fft.fftshift(np.fft.fftfreq(nchans, 1 / fs))
        self.freq_ind = np.arange(nchans)
        bw2 = fs / nchans / 2
        self.freq_l = self.freq_vec - bw2
        self.freq_h = self.freq_vec + bw2

        self.cfreqvec = np.fft.fftshift(np.fft.fftfreq(n_fpfb, 1 / c_bw))

    def get_ul_spec(self, data_vec, v_d, bMag, alphadeg=90.0, rcsflag=False, Tpe=1.0):

        Ne = data_vec[-1, 0]
        Te = data_vec[-1, -1]
        gam = 10000
        skirt = gam * 2

        lamb_d2 = sconst.epsilon_0 * sconst.k * Te / (sconst.e**2 * Ne)

        rcs_p = Tpe * Ne * self.K * 2 * lamb_d2 / 2

        frp, frm = get_freqs_default(Ne, Te, v_d, self.K, bMag, alphadeg)

        # Lower plasma line Spectrum
        f_0m = frm
        lb = self.freq_ind[f_0m - skirt > self.freq_l][-1]
        ub = self.freq_ind[f_0m + skirt < self.freq_h][0]
        cur_bin = self.freq_ind[lb : ub + 1]
        outlist = [frp, frm]
        f_lo = []
        spec_low = []
        for bin_i in cur_bin:
            cf_i = self.freq_vec[bin_i]
            f_lo.append(self.cfreqvec + cf_i)
            spec_i = make_pl_spec_default(self.cfreqvec + cf_i, f_0m, gam)
            spec_low.append(spec_i)

        # Upper plasma line spectrum.
        f_0p = frp
        lb = self.freq_ind[f_0p - skirt > self.freq_l][-1]
        ub = self.freq_ind[f_0p + skirt < self.freq_h][0]
        cur_bin = self.freq_ind[lb : ub + 1]
        outlist = [frp, frm]
        f_hi = []
        spec_hi = []
        for bin_i in cur_bin:
            cf_i = self.freq_vec[bin_i]
            f_hi.append(self.cfreqvec + cf_i)
            spec_i = make_pl_spec_default(self.cfreqvec + cf_i, f_0p, gam)
            spec_hi.append(spec_i)

        outlist = [f_lo, spec_low, f_hi, spec_hi]
        if rcsflag:
            outlist.append(rcs_p)

        return tuple(outlist)


def make_pl_spec_default(freq, f_0, gam):
    """Creates a Lorentzian function

    Parameters
    ----------
    freq : ndarray
        Frequency vector in Hz
    f_0 : float
        Center of Loerentzian in Hz.
    gam : float
        The HWHM of the function

    Returns
    -------
    lore : ndarray
        A lorentzian normalized to integrate to one.
    """

    lore = 1 / ((1 + (freq - f_0) ** 2 / gam**2) * gam * np.pi)
    return lore


def get_freqs_default(
    Ne,
    Te,
    v_d,
    k_num,
    bMag,
    alphadeg=90.0,
):
    """Calculate the the uper and lower plasma line frequences for Thompson scatter in the ionosphere.

    Parameters
    ----------
    Ne : float
        Electron density in m^-3
    Te : float
        Electron temperature in K
    v_d : float
        Doppler velocity in m/s
    k_num : float
        Bragg wave number in rad/m
    bMag : float
        Magnetic field in Teslas
    alphadeg : float
        The magnetic aspect angle in degrees.

    Returns
    -------
    frm : float
        Center frequency of lower plasma line in Hz.
    frp : float
        Center frequency of upper plasma line in Hz.
    """
    twopi = 2.0 * np.pi

    # debye length
    lamb_d2 = sconst.epsilon_0 * sconst.k * Te / (sconst.e**2 * Ne)
    # Electron cyclotron frequency
    f_c = sconst.e * bMag / (sconst.m_e * sconst.c) / twopi
    # plasma frequency squared
    f_p2 = sconst.e**2 * Ne / (sconst.epsilon_0 * sconst.m_e) / twopi**2
    # thermal speed squared
    f_th2 = lamb_d2 * f_p2

    alphrad = alphadeg * np.pi / 2
    # plasma frequency squared+3 *thermal frequency *K**2 + (fc*sinalph)**2
    fr_2 = f_p2 + 3 * k_num**2 * f_th2 + (f_c * np.sin(alphrad)) ** 2
    fr = np.sqrt(fr_2)
    f_o = k_num * sconst.c / (4 * np.pi)

    kp = (twopi / sconst.c) * (f_o + f_o + fr)
    km = (twopi / sconst.c) * (f_o + f_o - fr)

    # upper plasma line squared
    frp2 = f_p2 + 3 * kp**2 * f_th2 + (f_c * np.sin(alphrad)) ** 2
    # upper plasma line with doppler
    frp = np.sqrt(frp2) + (kp * v_d / (twopi))
    # lower plasma line squared
    frm2 = f_p2 + 3 * km**2 * f_th2 + (f_c * np.sin(alphrad)) ** 2
    # lower plasma line with doppler
    frm = -1 * np.sqrt(frm2) - (kp * v_d / (twopi))
    return frp, frm

from xmlrpc.client import boolean

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
        c_bw = fs / nchans
        n_fpfb = nfreq_pfb
        c_bw = fs / nchans

        self.freq_vec = np.fft.fftshift(np.fft.fftfreq(nchans, 1 / fs))
        self.freq_ind = np.arange(nchans)
        bw2 = fs / nchans / 2
        self.freq_l = self.freq_vec - bw2
        self.freq_h = self.freq_vec + bw2
        m_m = -((n_fpfb - 1) // -2)  # trick to get ceiling command with operators
        m_p = (n_fpfb - 1) // 2
        self.pos_dict = {
            ichan: (ichan * n_fpfb, (ichan + 1) * n_fpfb) for ichan in self.freq_ind
        }
        self.cfreqvec = np.fft.fftshift(np.fft.fftfreq(n_fpfb, 1 / c_bw))

    def get_ul_spec(
        self,
        data_vec,
        v_d,
        alphadeg=90.0,
        rcsflag=False,
        freqflag=False,
        Tpe=1.0,
        posflag=False,
        chanflag=False,
        heflag=False,
    ):
        """Gets the upper and lower plasma line spectra.

        Parameters
        ----------
        data_vec : ndarray
            A numpy array of Nsx2 array of densities and temperatures.
        v_d : float
            Doppler velocity in m/s
        bMag : float
            Magnetic field in Teslas
        alphadeg : float
            The magnetic aspect angle in degrees.
        rcsflag : bool
            A bool that will determine if the scaled density needed for RCS calculation is returned as well. (default is False)
        Tpe : float
            The ratio between plasma line temperature and electron temperature. Default is 1.0.
        chanflag : bool
            A bool that will determine if the channel entries are of the ouptuts. (default is False)
        posflag : bool
            A flag to give out the position within the full frequency space of the spectra.
        heflag : bool
            A bool that will determine if the debye length is one of the ouptuts. (default is False)
            
        Returns
        -------
        f_lo : list
            List of frequency vectors parameterizing the lower plasma line spectra in Hz.
        spec_low : list
            List of lower plasma line spectra.
        f_hi : list
            List of frequency vectors parameterizing the upper plasma line spectra in Hz.
        spec_hi : list
            List of upper plasma line spectra.
        rcs_p : float
            The scaled density needed to calculate the RCS in m^{-3}.
        h_e : float
             The debye length in m.
        freqlist : list
            Center frequencies in of lower and upper plasma lines in Hz.
        chan_nums_l : list
            List of channels that are assoicated with the output
        """

        Ne = data_vec[-1, 0]
        Te = data_vec[-1, -1]
        # HACK half with half maximum in Hz set to 10 kHz.
        gam = 10000
        # Skirt is FWHM
        skirt = gam * 2

        lamb_d2 = sconst.epsilon_0 * sconst.k * Te / (sconst.e**2 * Ne)

        rcs_p = Tpe * Ne * self.K * 2 * lamb_d2 / 2

        frp, frm = get_pl_freqs_default(Ne, Te, v_d, self.K, self.bMag, alphadeg)

        # Lower plasma line Spectrum
        f_0m = frm
        lb = self.freq_ind[f_0m - skirt > self.freq_l][-1]
        ub = self.freq_ind[f_0m + skirt < self.freq_h][0]
        # Find the bins that have data in them
        cur_bin = self.freq_ind[lb : ub + 1]
        outlist = [frp, frm]
        f_lo = []
        spec_low = []
        pos_low = []
        chan_nums_l = []
        for bin_i in cur_bin:
            cf_i = self.freq_vec[bin_i]
            f_lo.append(self.cfreqvec + cf_i)
            spec_i = make_pl_spec_default(self.cfreqvec + cf_i, f_0m, gam)
            spec_low.append(spec_i)
            pos_low.append(np.arange(*self.pos_dict[bin_i], dtype=int))
            chan_nums_l.append(bin_i)
            
        # Upper plasma line spectrum.
        f_0p = frp
        lb = self.freq_ind[f_0p - skirt > self.freq_l][-1]
        ub = self.freq_ind[f_0p + skirt < self.freq_h][0]
        cur_bin = self.freq_ind[lb : ub + 1]
        outlist = [frp, frm]
        f_hi = []
        spec_hi = []
        pos_hi = []
        chan_nums_h = []
        for bin_i in cur_bin:
            cf_i = self.freq_vec[bin_i]
            f_hi.append(self.cfreqvec + cf_i)
            spec_i = make_pl_spec_default(self.cfreqvec + cf_i, f_0p, gam)
            spec_hi.append(spec_i)
            pos_hi.append(np.arange(*self.pos_dict[bin_i], dtype=int))
            chan_nums_h.append(bin_i)
            
        outlist = [f_lo, spec_low, f_hi, spec_hi]
        if rcsflag:
            outlist.append(rcs_p)
        if heflag:
            outlist.append(np.sqrt(lamb_d2))
        if freqflag:
            outlist.append([frm, frp])
        if chanflag: 
            outlist.append(chan_nums_l)
            outlist.append(chan_nums_h)
        if posflag:
            outlist.append(pos_low)
            outlist.append(pos_hi)

        return tuple(outlist)


def make_pl_spec_default(freq, f_0, gam):
    """Creates a Lorentzian function using the center frequency and half width at half maximum (HWHM) parameters in Hz.

    Parameters
    ----------
    freq : ndarray
        Frequency vector in Hz
    f_0 : float
        Center of Loerentzian in Hz.
    gam : float
        The HWHM of the function in Hz.

    Returns
    -------
    lore : ndarray
        A lorentzian normalized to integrate to one.
    """

    lore = 1 / ((1 + (freq - f_0) ** 2 / gam**2) * gam * np.pi)
    return lore


def get_pl_freqs_default(
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

    phirad = np.abs(90 - alphadeg) * np.pi / 180
    # plasma frequency squared+3 *thermal frequency *K**2 + (fc*sinalph)**2
    fr_2 = f_p2 + 3 * k_num**2 * f_th2 + (f_c * np.sin(phirad)) ** 2
    fr = np.sqrt(fr_2)
    f_o = k_num * sconst.c / (4 * np.pi)

    kp = (twopi / sconst.c) * (f_o + f_o + fr)
    km = (twopi / sconst.c) * (f_o + f_o - fr)

    # upper plasma line squared
    frp2 = f_p2 + 3 * kp**2 * f_th2 + (f_c * np.sin(phirad)) ** 2
    # upper plasma line with doppler
    frp = np.sqrt(frp2) + (kp * v_d / (twopi))
    # lower plasma line squared
    frm2 = f_p2 + 3 * km**2 * f_th2 + (f_c * np.sin(phirad)) ** 2
    # lower plasma line with doppler
    frm = -1 * np.sqrt(frm2) - (kp * v_d / (twopi))
    return frp, frm

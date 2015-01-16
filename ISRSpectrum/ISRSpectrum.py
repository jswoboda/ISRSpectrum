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
import numpy as np
import scipy.special
import pdb
from const.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge, v_me, v_amu

class ISRSpectrum(object):
    """ Class to create the spectrum. The instance of the class will hold infomation on
    the radar system such as sample frequency, center frequency and number of points for
    the spectrum.
    Parameters
    bMag - The magnetic field magnitude in Teslas.
    K - The Bragg scattering vector magnitude corresponds to 1/2 radar wavelength.
    f - Vector holding the frequency values in Hz
    omeg - Vector hold the frequency values in rad/s"""
    def __init__(self,centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=64, sampfreq=50e3):
        """ Constructor for the class.
        Inputs :
        centerFrequency: The radar center frequency in Hz.
        bMag: The magnetic field magnitude in Teslas.
        nspec: the number of points of the spectrum.
        sampfreq: The sampling frequency of the A/Ds"""
        self.bMag = bMag

        self.K = 2.0*np.pi*2*centerFrequency/v_C_0 #The Bragg scattering vector, corresponds to half the radar wavelength.

        self.f = np.arange(-np.ceil((nspec-1.0)/2.0),np.floor((nspec-1.0)/2.0+1))*(sampfreq/(2*np.ceil((nspec-1.0)/2.0)))
        self.omeg = 2.0*np.pi*self.f

    def getspec(self,datablock,alphadeg=90.0):
        """ Gives the spectrum and the frequency vectors given the block of data and
        the magnetic aspect angle.
        Inputs
        datablock: A numpy array of size 1+Nionsx6 that holds the plasma parameters needed
        to create the spectrum. The first row will hold the information for the electrons.
        Each row of the array will have the following set up.
            [Ns, Ts, Vs, qs, ms, nus]
            Ns - The density of the species in m^-3
            Ts - Tempretur of the species in degrees K
            Vs - The Doppler velocity in m/s.
            qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
            ms - Mass of the species in AMU. (Value will be replaced for the electrons)
            nus - Collision frequency for species in s^-1.
        alphadeg: The magnetic aspect angle in degrees.
        Outputs
        freqvec: The frequency vector in Hz.
        spec: the spectrum in a numpy array.
        """
        alpha = alphadeg*np.pi/180
        estuff = datablock[0]
        estuff[3] = -v_elemcharge
        estuff[4] = v_me
        ionstuff = datablock[1:]
        (egord,h_e,Ne,omeg_e) = self.__calcgordeyev__(estuff,alpha)

        sig_e = (1-1j*omeg_e*egord)/(self.K**2*h_e**2)
        nte = 2*Ne*np.real(egord)

        #adjust ion stuff
        ionstuff[:,3] = ionstuff[:,3]*v_elemcharge
        ionstuff[:,4] = ionstuff[:,4]*v_amu
        ionden = np.sum(ionstuff[:,0])
        ionstuff[:,0] = (estuff[0]/ionden)*ionstuff[:,0]
        firstion = True
        for iinfo in ionstuff:
            (igord,h_i,Ni,omeg_i) = self.__calcgordeyev__(iinfo,alpha)

            sig_i = (1-1j*omeg_i*igord)/(self.K**2*h_i**2)
            nti = 2*Ni*np.real(igord)
            if firstion:
                sig_sum =sig_i
                nt_sum = nti
            else:
                sig_sum = sig_i+sig_sum
                nt_sum = nti+nt_sum

        inum = np.abs(sig_e)**2*nt_sum
        enum = np.abs(1+sig_sum)**2*nte
        den = np.abs(1 + sig_sum +sig_e)**2
        iline = inum/den
        eline = enum/den
        spec = iline+eline
        return (self.f,spec)

    def __calcgordeyev__(self,dataline,alpha):
        """ Performs the Gordeyve integral calculation.
        Inputs
        dataline: A numpy array of length that holds the plasma parameters needed
        to create the spectrum.
        Each row of the array will have the following set up.
            [Ns, Ts, Vs, qs, ms, nus]
            Ns - The density of the species in m^-3
            Ts - Tempretur of the species in degrees K
            Vs - The Doppler velocity in m/s.
            qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
            ms - Mass of the species in AMU. (Value will be replaced for the electrons)
            nus - Collision frequency for species in s^-1.
        alphadeg: The magnetic aspect angle in radians.
        Output
        gord: The result of the Gordeyev integral over Doppler corrected radian frequency
        hs: The Debye length in m.
        Ns: The density of the species in m^-3
        omeg_s: An array of the Doppler corrected radian frequency
        """
        aldiff = np.abs(alpha*180.0/np.pi-90.0)<1;
        K = self.K

        (Ns,Ts,Vs,qs,ms,nus) = dataline[:7]
        hs = np.sqrt(v_epsilon0*v_Boltz*Ts/(Ns*qs*qs))
        C = np.sqrt(v_Boltz*Ts/ms)
        omeg_s = self.omeg - K*Vs
        theta = omeg_s/(K*C*np.sqrt(2.0))


        if K*C>10.0*nus and aldiff:
            #for case with no collisions and magnetic field
            gord = (np.sqrt(np.pi)*np.exp(-theta**2)-1j*2.0*scipy.special.dawsn(theta))/(K*C*np.sqrt(2))


        return (gord,hs,Ns,omeg_s)
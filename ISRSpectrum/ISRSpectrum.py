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
from const.mathutils import sommerfelderfrep
class ISRSpectrum(object):
    """ Class to create the spectrum. The instance of the class will hold infomation on
    the radar system such as sample frequency, center frequency and number of points for
    the spectrum.
    Parameters
    bMag - The magnetic field magnitude in Teslas.
    K - The Bragg scattering vector magnitude corresponds to 1/2 radar wavelength.
    f - Vector holding the frequency values in Hz
    omeg - Vector hold the frequency values in rad/s
    self.collfreqmin - Minimum collision frequency before they are taken into account in the Gordeyev integral calculation
    self.alphamax - The maximum aspect angle
    self.dFlag - A debug flag."""
    def __init__(self,centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=64, sampfreq=50e3,collfreqmin=1e-2,alphamax=30.0,dFlag=False):
        """ Constructor for the class.
        Inputs :
        centerFrequency: The radar center frequency in Hz.
        bMag: The magnetic field magnitude in Teslas.
        nspec: the number of points of the spectrum.
        sampfreq: The sampling frequency of the A/Ds in Hz
        collfreqmin: (Default 1e-2) The minimum collision frequency needed to incorporate it into Gordeyev
            integral calculations in units of K*sqrt(Kb*Ts/ms) for each ion species.
        alphamax: (Default 30)  The maximum magnetic aspect angle in which the B-field will be taken into account.
        dFlag: A debug flag, if set true will output debug text. Default is false."""
        self.bMag = bMag
        self.dFlag = dFlag
        self.collfreqmin = collfreqmin
        self.alphamax = alphamax

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
        # perform a copy of the object to avoid it being written over with incorrect.
        datablock=datablock.copy()
        dFlag = self.dFlag
        alpha = alphadeg*np.pi/180
        estuff = datablock[0]
        estuff[3] = -v_elemcharge
        estuff[4] = v_me
        ionstuff = datablock[1:]
        if dFlag:
            print "Calculating Gordeyev int for electons"
        (egord,h_e,Ne,omeg_e) = self.__calcgordeyev__(estuff,alpha)

        sig_e = (1-1j*omeg_e*egord)/(self.K**2*h_e**2)
        nte = 2*Ne*np.real(egord)

        #adjust ion stuff
        ionstuff[:,3] = ionstuff[:,3]*v_elemcharge
        ionstuff[:,4] = ionstuff[:,4]*v_amu
        ionden = np.sum(ionstuff[:,0])
        ionstuff[:,0] = (estuff[0]/ionden)*ionstuff[:,0]
        firstion = True
        for iion,iinfo in enumerate(ionstuff):
            if dFlag:
                print "Calculating Gordeyev int for ion species #{:d}".format(iion)
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

    def __calcgordeyev__(self,dataline,alpha,alphdiff = 10.0):
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
        dFlag = self.dFlag
        K = self.K

        (Ns,Ts,Vs,qs,ms,nus) = dataline[:7]
        hs = np.sqrt(v_epsilon0*v_Boltz*Ts/(Ns*qs*qs))
        C = np.sqrt(v_Boltz*Ts/ms)
        omeg_s = self.omeg - K*Vs
        theta = omeg_s/(K*C*np.sqrt(2.0))
        Om = qs*self.bMag/ms
        # determine what integral is used
        magbool = alpha*180.0/np.pi < self.alphamax
        collbool = self.collfreqmin*K*C<nus

#        pdb.set_trace()
        if  not collbool and not magbool:
            #for case with no collisions or magnetic field just use analytic method
            gord = (np.sqrt(np.pi)*np.exp(-theta**2)-1j*2.0*scipy.special.dawsn(theta))/(K*C*np.sqrt(2))
            if dFlag:
                print '\t No collisions No magnetic field'
            return (gord,hs,Ns,omeg_s)
        elif collbool and not magbool:
            if dFlag:
                print '\t With collisions No magnetic field'
            gordfunc = collacf
            exparams = (K,C,nus)
        elif not collbool and magbool:
            if dFlag:
                print '\t No collisions with magnetic field'
            gordfunc = magacf
            exparams = (K,C,alpha,Om)
        else:
            if dFlag:
                print '\t With collisions with magnetic field'
            gordfunc = magncollacf
            exparams = (K,C,alpha,Om,nus)

        N_somm = 2**9
        b1 = 100.0/(K*C*np.sqrt(2.0))

        (gord,flag_c,outrep) = sommerfelderfrep(gordfunc,N_somm,omeg_s,b1,Lmax=100,errF=1e-2,exparams=exparams)
        if dFlag:
            yna = ['No','Yes']
            print '\t Converged: {:s}, Number of iterations: {:d}'.format(yna[flag_c],outrep)

        return (gord,hs,Ns,omeg_s)

def magacf(tau,K,C,alpha,Om):
    """ magacf(tau,K,C,alpha,Om)
        by John Swoboda
        This function will create a single particle acf for a particle species with magnetic
        field but no collisions.
        Inputs
        tau: The time vector for the acf.
        K: Bragg scatter vector magnetude.
        C: Thermal speed of the species.
        alpha: Magnetic aspect angle in radians.
        Om: The gyrofrequency of the particle.
        Output
        acf - The single particle acf.
        """
    Kpar = np.sin(alpha)*K
    Kperp = np.cos(alpha)*K
    return np.exp(-np.power(C*Kpar*tau,2.0)/2.0-2.0*np.power(Kperp*C*np.sin(Om*tau/2.0)/Om,2.0))
def collacf(tau,K,C,nu):
    """ collacf(tau,K,C,nu)
        by John Swoboda
        This function will create a single particle acf for a particle species with no magnetic
        field but with collisions.
        Inputs
        tau: The time vector for the acf.
        K: Bragg scatter vector magnetude.
        C: Thermal speed of the species.
        nu: The collision frequency in collisions/sec
        Output
        acf - The single particle acf.
        """
    return np.exp(-np.power(K*C/nu,2.0)*(nu*tau-1+np.exp(-nu*tau)))
def magncollacf(tau,K,C,alpha,Om,nu):
    """ magncollacf(tau,K,C,alpha,Om)
        by John Swoboda
        This function will create a single particle acf for a particle species with magnetic
        field and collisions.
        Inputs
        tau: The time vector for the acf.
        K: Bragg scatter vector magnetude.
        C: Thermal speed of the species.
        alpha: Magnetic aspect angle in radians.
        Om: The gyrofrequency of the particle.
        nu: The collision frequency in collisions/sec
        Output
        acf - The single particle acf.
        """
    Kpar = np.sin(alpha)*K
    Kperp = np.cos(alpha)*K
    gam = np.arctan(nu/Om)

    deltl = np.exp(-np.power(Kpar*C/nu,2.0)*(nu*tau-1+np.exp(-nu*tau)))
    deltp = np.exp(-np.power(C*Kperp,2.0)/(Om*Om+nu*nu)*(np.cos(2*gam)+nu*tau-np.exp(-nu*tau)*(np.cos(Om*tau-2.0*gam))))
    return deltl*deltp

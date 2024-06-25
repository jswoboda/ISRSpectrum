#!python

"""
common_examples.py

A set of examples that creates figures showing how the spectra change based off of different parameters varying.
"""

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import scipy.fft as scfft
sns.set_style("whitegrid")
sns.set_context("notebook")
#
from ISRSpectrum import Specinit


def temp_ratio():
    """Creates a figure showing the affect of temperature ratio on the spectra. The base spectrum is from a 100% O+ plasma, with an electron density of 1e11 and the ion temperature is 1000 K; the radar has a center frequency of 440MHz. The ratio varied between 1-3.    
    """
    nspec = 4097
    spfreq = 200e3
    ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=nspec, sampfreq=spfreq,dFlag=True)

    ratio_var = np.linspace(1,3,5)
    Ne = 1e11
    Ni = 1e11
    Ti  = 1000
    species=['O+','e-']
    n_lag = int((nspec-1)/2)
    fig,ax = plt.subplots(2,1,figsize=(8, 6))
    for ir in ratio_var:

        datablock = np.array([[Ni,Ti],[Ne,ir*Ti]])
        species = ['O+','e-']
        (f,spec1,rcs) = ISS2.getspecsep(datablock,species,vel = 0.0, alphadeg=90.0,rcsflag=True)

        acf=scfft.ifft(scfft.ifftshift(spec1)).real
        tau=scfft.ifftshift(np.arange(-np.ceil((float(nspec)-1)/2),np.floor((float(nspec)-1)/2)+1))/spfreq
        
        ax[0].plot(f*1e-3,spec1/np.nanmax(spec1), linestyle='-',linewidth=3,label=r'Tr={0}'.format(ir))
        l1=ax[1].plot(tau[:n_lag]*1e6,acf[:n_lag]/acf[0],'-',lw=3,label=r'Tr={0}'.format(ir))[0]
        sns.despine()
   
    ax[0].set_xlabel('f in kHz', fontsize=14)
    ax[0].set_ylabel('Amp',fontsize=14)
    ax[0].set_title('Spectra',fontsize=18)
    ax[0].set_xlim([-20,20])
    ax[0].set_ylim([0,1.1])
    
    ax[1].set_xlim([0,480])
    ax[1].spines['right'].set_visible(False)

    ax[1].set_xlabel('Lag in Microseconds ',fontsize=14)
    ax[1].set_title('ACF',fontsize=18)
    ax[1].set_ylabel('Normalized Magnitude',fontsize=14)
    ax[1].set_ylim([-.5,1.1])
    ax[1].legend()
    plt.tight_layout()
    plt.savefig('Tempratiovary.png',dpi=300)



def ion_temp_vary():
    """Creates a figure showing the affect of ion temperature on the spectra. The base spectrum is from a 100% O+ plasma, with electron density of 1e11, the Te/Ti ratio is fixed at 2; the radar has a center frequency of 440MHz. The ratio varies from 1-3.
    """
    nspec = 4097
    spfreq = 200e3
    ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=nspec, sampfreq=spfreq,dFlag=True)

    ti_list = np.linspace(500,2500,5)
    Ne = 1e11
    Ni = 1e11
    species=['O+','e-']
    n_lag = int((nspec-1)/2)
    fig,ax = plt.subplots(2,1,figsize=(8, 6))
    for ti in ti_list:

        datablock = np.array([[Ni,ti],[Ne,2*ti]])
        (f,spec1,rcs) = ISS2.getspecsep(datablock,species,vel = 0.0, alphadeg=90.0,rcsflag=True)

        acf=scfft.ifft(scfft.ifftshift(spec1)).real
        tau=scfft.ifftshift(np.arange(-np.ceil((float(nspec)-1)/2),np.floor((float(nspec)-1)/2)+1))/spfreq
        
        ax[0].plot(f*1e-3,spec1/np.nanmax(spec1), linestyle='-',linewidth=3,label="{0} K".format(int(ti)))
        l1=ax[1].plot(tau[:n_lag]*1e6,acf[:n_lag]/acf[0],'-',lw=3,label=r'{0:d} K'.format(int(ti)))[0]
        sns.despine()
   
    ax[0].set_xlabel('f in kHz', fontsize=14)
    ax[0].set_ylabel('Amp',fontsize=14)
    ax[0].set_title('Spectra',fontsize=18)
    ax[0].set_xlim([-20,20])
    ax[0].set_ylim([0,1.1])
    
    ax[1].set_xlim([0,480])
    ax[1].spines['right'].set_visible(False)

    ax[1].set_xlabel('Lag in Microseconds ',fontsize=14)
    ax[1].set_title('ACF',fontsize=18)
    ax[1].set_ylabel('Normalized Magnitude',fontsize=14)
    ax[1].set_ylim([-.5,1.1])
    ax[1].legend()
    plt.tight_layout()
    plt.savefig('Iontempvary.png',dpi=300)


def collision_frequency():
    """Creates a figure showing the affect of collision on the spectra. The base spectrum is from a 100% NO+ plasma, with electron density of 1e12, Ti=Te=500 K; the radar has a center frequency of 440MHz. The collision frequency varies from 10 Hz to 10 kHz
    """

    nspec = 4096
    spfreq = 100e3
    ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=nspec, sampfreq=spfreq,dFlag=True)

    col_freq = np.array([10,100,1000,5000,10000])
    Ne = 1e12
    Ti = 500
    Te = 500
    species=['NO+','e-']
    n_lag = int((nspec-1)/2)
    fig,ax = plt.subplots(2,1,figsize=(8, 6))
    # HACK the charge and amu of the electron will be replaced is the getspec method
    datablock = np.array([[Ne,Ti,0,1,30,0 ],[Ne,Te,0,-1,1,0]])
    for icol in col_freq:

        datablock[0,-1] =icol
        print("collision frequency {0} kHz".format(icol*1e-3))
        (f,spec1) = ISS2.getspec(datablock,des_plug="default")

        acf=scfft.ifft(scfft.ifftshift(spec1)).real
        tau=scfft.ifftshift(np.arange(-np.ceil((float(nspec)-1)/2),np.floor((float(nspec)-1)/2)+1))/spfreq
        
        ax[0].plot(f*1e-3,spec1/np.nanmax(spec1), linestyle='-',linewidth=3,label="col freq = {0} kHz".format(icol*1e-3))
        l1=ax[1].plot(tau[:n_lag]*1e6,acf[:n_lag]/acf[0],'-',lw=3,label=r'col freq = {0} kHz'.format(icol*1e-3))[0]
        sns.despine()
   
    ax[0].set_xlabel('f in kHz', fontsize=14)
    ax[0].set_ylabel('Amp',fontsize=14)
    ax[0].set_title('Spectra',fontsize=18)
    ax[0].set_xlim([-10,10])
    ax[0].set_ylim([0,1.1])
    ax[0].legend()

    ax[1].set_xlim([0,480])
    ax[1].spines['right'].set_visible(False)

    ax[1].set_xlabel('Lag in Microseconds ',fontsize=14)
    ax[1].set_title('ACF',fontsize=18)
    ax[1].set_ylabel('Normalized Magnitude',fontsize=14)
    ax[1].set_ylim([-.5,1.1])
    
    plt.tight_layout()
    plt.savefig('colfreqvary.png',dpi=300)
 

def species_vary():
    """Creates a figure showing the affect of different species of plasma spectra. The base spectrum is a mix of O+ and NO+ plasma, with electron density of 1e12, 2Ti=Te=3000 K; the radar has a center frequency of 440MHz. The collision frequency varies from species varies from all NO+ to all O+
    """

    nspec = 4096
    spfreq = 100e3
    ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=nspec, sampfreq=spfreq,dFlag=True)

    comp_frac = np.linspace(0.,1.,6)
    Ne = 1e12
    Ti = 1500
    Te = 3000
    species=['NO+',"O+",'e-']
    n_lag = int((nspec-1)/2)
    fig,ax = plt.subplots(2,1,figsize=(8, 6))
    # HACK the charge and amu of the electron will be replaced is the getspec method
    datablock = np.array([[Ne,Ti,0,1,30,0 ],[Ne,Ti,0,1,16,0 ],[Ne,Te,0,-1,1,0]])
    for icomp in comp_frac:
        no_c = 1-icomp
        o_c = icomp
        datablock[0,0] = Ne*no_c
        datablock[1,0] = Ne*o_c

        (f,spec1) = ISS2.getspec(datablock,des_plug="default")

        acf=scfft.ifft(scfft.ifftshift(spec1)).real
        tau=scfft.ifftshift(np.arange(-np.ceil((float(nspec)-1)/2),np.floor((float(nspec)-1)/2)+1))/spfreq
        
        ax[0].plot(f*1e-3,spec1/np.nanmax(spec1), linestyle='-',linewidth=3,label="O+/e-= {0:.1f}".format(icomp))
        l1=ax[1].plot(tau[:n_lag]*1e6,acf[:n_lag]/acf[0],'-',lw=3,label="O+/e-= {0:.1f}".format(icomp))[0]
        sns.despine()
   
    ax[0].set_xlabel('f in kHz', fontsize=14)
    ax[0].set_ylabel('Amp',fontsize=14)
    ax[0].set_title('Spectra',fontsize=18)
    ax[0].set_xlim([-20,20])
    ax[0].set_ylim([0,1.1])    
    ax[1].set_xlim([0,480])
    ax[0].legend()


    ax[1].spines['right'].set_visible(False)

    ax[1].set_xlabel('Lag in Microseconds ',fontsize=14)
    ax[1].set_title('ACF',fontsize=18)
    ax[1].set_ylabel('Normalized Magnitude',fontsize=14)
    ax[1].set_ylim([-.5,1.1])
    
    plt.tight_layout()
    plt.savefig('compvary.png',dpi=300)
 



 
if __name__== '__main__':
    temp_ratio()
    ion_temp_vary()
    collision_frequency()
    species_vary()
#!/usr/bin/env python
"""
elineexample.py
Created on Sun Dec 27 15:28:19 2015
This example shows everything up to electron line for magnitized and non-magnitized plasmas.
@author: John Swoboda
"""

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import scipy.fft as scfft
sns.set_style("whitegrid")
sns.set_context("notebook")
#
from ISRSpectrum import Specinit

if __name__== '__main__':

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
        species = ['O+','e-']
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

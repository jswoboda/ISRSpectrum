#!/usr/bin/env python
"""
Created on Fri Apr 22 20:24:52 2016

@author: John Swoboda
"""

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("notebook")
#
from ISRSpectrum.ISRSpectrum import ISRSpectrum


if __name__== '__main__':

   
    f = np.logspace(1,np.log10(2.5e6),2**10)
#    f = np.logspace(1,np.log10(8e3),2**10)
    ISpec = ISRSpectrum(nspec=2**16,bMag = 3.5e-5,sampfreq=15e6,alphamax=80,f=f,dFlag=True)

    species=['O+','e-']
#    databloc = np.array([[1.66e10,863.],[1.66e10,863.]])
    databloc = np.array([[1.66e10,1e3],[1.66e10,2.5e3]])
    #%% With B-Field
    aldeg = 30.
    f,[iline,eline] = ISpec.getspecsep(databloc,species,alphadeg=aldeg,seplines=True)

    plt.figure()
    ired = 100.
    spec1=iline/ired+eline
    maxy=spec1.max()
    maxi = iline.max()/ired
    f =f*1e-6
    l1=plt.plot(f,spec1,'-',label='Sum',lw=3)[0]
    l2=plt.plot(f,iline/ired,'g--',label='Ion Line',lw=3)[0]
    l3=plt.plot(f,eline,'r--',label='Elec. Line',lw=3)[0]

    plt.xlabel('Frequency in MHz')

    plt.ylabel('Amplitude')

    plt.xscale('linear')

    plt.yscale('linear')
#    plt.xlim(f.min(),8.)
    plt.ylim((1,maxi*1.2))
    plt.legend(loc='upper right',handles=[l1,l2,l3])

    plt.grid(False)

    plt.title('Spectrum with B-field at alpha at {0} deg'.format(int(aldeg)))
    plt.savefig('elinilinemagfinallinlin.png',dpi=300)

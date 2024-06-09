#!/usr/bin/env python
"""
elineexample.py
Created on Sun Dec 27 15:28:19 2015
This example shows everything up to electron line for magnitized and non-magnitized plasmas.
@author: John Swoboda
"""

import numpy as np

from ISRSpectrum import Specinit
import matplotlib.pylab as plt
import seaborn as sns


if __name__== '__main__':
    sns.set_style("whitegrid")
    sns.set_context("notebook")


#    f = np.lin(1,7,2**10)
    f = np.linspace(0.,5e6,2**12)

    ISpec = Specinit(centerFrequency = 440e6,nspec=2**16,sampfreq=15e6,alphamax=60,f=f,dFlag=True)

    species=['O+','e-']
    databloc = np.array([[1e11,1100.],[1e11,2500.]])
    #%% No B-Field
    f,[iline,eline] = ISpec.getspecsep(databloc,species,alphadeg=90,seplines=True)

    plt.figure()

    spec1=iline+eline
    maxy=spec1.max()

    l1=plt.plot(f,spec1,'-',label='Sum',lw=3)[0]
    l2=plt.plot(f,iline,':',label='Ion Line',lw=3)[0]
    l3=plt.plot(f,eline,'--',label='Elec. Line',lw=3)[0]

    plt.xlabel('Frequency in Hz')

    plt.ylabel('Amplitude')

    plt.xscale('log')

    plt.yscale('log')
    plt.ylim((maxy*1e-7,maxy))
    plt.legend(loc='lower left',handles=[l1,l2,l3])

    plt.grid(True)

    plt.title('Spectrum in Log Log No B-field')

    plt.savefig('elinilinenoB.png',dpi=300)
    #%% With B-Field
    aldeg = 15
    f,[iline,eline] = ISpec.getspecsep(databloc,species,alphadeg=aldeg,seplines=True)

    plt.figure()

    spec1=iline+eline
    maxy=spec1.max()

    l1=plt.plot(f,spec1,'-',label='Sum',lw=3)[0]
    l2=plt.plot(f,iline,':',label='Ion Line',lw=3)[0]
    l3=plt.plot(f,eline,'--',label='Elec. Line',lw=3)[0]

    plt.xlabel('Frequency in Hz')

    plt.ylabel('Amplitude')

    plt.xscale('log')

    plt.yscale('log')
    plt.ylim((maxy*1e-7,maxy))
    plt.legend(loc='lower left',handles=[l1,l2,l3])

    plt.grid(True)

    plt.title('Spectrum in Log Log With B-field at alpha at {0} deg'.format(int(aldeg)))
    plt.savefig('elinilinewB.png',dpi=300)

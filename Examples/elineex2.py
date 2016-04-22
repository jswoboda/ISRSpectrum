#!/usr/bin/env python
"""
elineexample2.py
Created on Sun Dec 27 15:28:19 2015
This example shows everything up to electron line for magnitized and non-magnitized plasmas.
@author: John Swoboda
"""

import numpy as np
import os,inspect
from ISRSpectrum.ISRSpectrum import ISRSpectrum
import matplotlib.pylab as plt
import seaborn as sns


if __name__== '__main__':
    sns.set_style("whitegrid")
    sns.set_context("notebook")
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    imagepath = os.path.join(os.path.split(curpath)[0],'Doc','Figs')

    f = np.linspace(0.,5e6,2**12)
    ISpec = ISRSpectrum(nspec=2**16,sampfreq=15e6,alphamax=60,f=f,dFlag=True)

    species=['O+','e-']
    databloc = np.array([[1e11,1100.],[1e11,2500.]])
    #%% No B-Field
    f,[iline,eline] = ISpec.getspecsep(databloc,species,alphadeg=90,seplines=True)

    plt.figure()

    spec1=iline+eline
    maxy=spec1.max()

    l1=plt.plot(f*1e-6,spec1,'-',label='Sum',lw=3)[0]
    l2=plt.plot(f*1e-6,iline,':',label='Ion Line',lw=3)[0]
    l3=plt.plot(f*1e-6,eline,'--',label='Elec. Line',lw=3)[0]

    plt.xlabel('Frequency in MHz')

    plt.ylabel('Amplitude')


    plt.yscale('log')
    plt.ylim((maxy*1e-7,maxy))
    plt.legend(loc='lower left',handles=[l1,l2,l3])

    plt.grid(True)

    plt.title('Spectrum in Log Y scale No B-field')

    plt.savefig(os.path.join(imagepath,'eliniline2noB.png'),dpi=300)
    #%% With B-Field
    aldeg = 20
    f,[iline,eline] = ISpec.getspecsep(databloc,species,alphadeg=aldeg,seplines=True)

    plt.figure()

    spec1=iline+eline
    maxy=spec1.max()

    l1=plt.plot(f*1e-6,spec1,'-',label='Sum',lw=3)[0]
    l2=plt.plot(f*1e-6,iline,':',label='Ion Line',lw=3)[0]
    l3=plt.plot(f*1e-6,eline,'--',label='Elec. Line',lw=3)[0]

    plt.xlabel('Frequency in MHz')

    plt.ylabel('Amplitude')


    plt.yscale('log')
    plt.ylim((maxy*1e-7,maxy))
    plt.legend(loc='lower left',handles=[l1,l2,l3])

    plt.grid(True)

    plt.title('Spectrum in Log Y scale With B-field at alpha at {0} deg'.format(int(aldeg)))
    plt.savefig(os.path.join(imagepath,'eliniline2wB.png'),dpi=300)

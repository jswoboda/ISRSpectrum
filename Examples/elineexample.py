#!/usr/bin/env python
"""
elineexample.py

This set of eamples show how to break up the electron line and ion line from the spectrum.
"""

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
#
from ISRSpectrum import Specinit



def elineplot():
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



    f = np.linspace(20e3,2.5e6,2**10)
    ISpec = Specinit(nspec=2**16,bMag = 3.5e-5,sampfreq=15e6,alphamax=80,f=f,dFlag=True)

    species=['O+','e-']
    databloc = np.array([[1.66e10,863.],[1.66e10,863.]])
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

    plt.savefig('eliniline2noB.png',dpi=300)
    #%% With B-Field
    aldeg = 30
    f,[iline,eline] = ISpec.getspecsep(databloc,species,alphadeg=aldeg,seplines=True)

    plt.figure()

    spec1=iline+eline
    maxy=spec1.max()

    l1=plt.plot(f*1e-6,spec1,'-',label='Sum',lw=3)[0]
    l2=plt.plot(f*1e-6,iline,':',label='Ion Line',lw=3)[0]
    l3=plt.plot(f*1e-6,eline,'--',label='Elec. Line',lw=3)[0]

    plt.xlabel('Frequency in MHz')

    plt.ylabel('Amplitude')


    plt.ylim((-maxy*1e-7,maxy*1.1))
    plt.legend(loc='lower left',handles=[l1,l2,l3])

    plt.grid(True)

    plt.title('Spectrum in Log Y scale With B-field at alpha at {0} deg'.format(int(aldeg)))
    plt.savefig('eliniline2wB.png',dpi=300)


    f = np.logspace(1,np.log10(2.5e6),2**10)
#    f = np.logspace(1,np.log10(8e3),2**10)
    ISpec = Specinit(nspec=2**16,bMag = 3.5e-5,sampfreq=15e6,alphamax=80,f=f,dFlag=True)

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

if __name__== '__main__':
    elineplot()

#!/usr/bin/env python
"""
Created on Tue Feb  3 09:39:46 2015

@author: John Swoboda
"""

import numpy as np
import matplotlib.pylab as plt
import RadarDataSim.ISSpectrum as ISSorig
import ISRSpectrum.ISRSpectrum as ISSnew
from ISRSpectrum.const.physConstants import v_Boltz, v_me, v_amu #v_C_0, v_epsilon0, v_elemcharge,



if __name__== '__main__':
    ISS1 = ISSorig.ISSpectrum(centerFrequency = 440.2, bMag = 0.4e-4, nspec=129, sampfreq=50e3)
    ISS2 = ISSnew.ISRSpectrum(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=129, sampfreq=50e3,dFlag=True)

    ti = 1e3
    te = 1e3
    Ne = 1e11
    mi = 16.0
    (omega,specorig) = ISS1.getSpectrum(ti,te/ti,np.log10(Ne),mi,1,0)
    plt.figure()
    plt.plot(omega,specorig,marker='o', linestyle='--')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Original Spectrum Program')
    plt.grid(True)

    plt.show(False)

    datablock = np.array([[Ne,te,0,1,1,0],[Ne,ti,0,1,mi,0]])

    (omega,specnew) = ISS2.getspec(datablock)
    plt.figure()
    plt.plot(omega,specnew,marker='o', linestyle='--')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Ion Line')
    plt.grid(True)

    plt.show(False)

    specnewred = specnew/np.nanmax(specnew)

    plt.figure()
    plt.plot(omega,specorig,marker='.', linestyle='--',color='b',label='Old')
    plt.hold(True)
    plt.plot(omega,specnewred,marker='o', linestyle='--',color='r',label='New')
    plt.hold(False)
    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Ion Line')
    plt.grid(True)
    #%% Vary temp keep ratio

    tevec = np.linspace(500.0,3000,4)
    figvarTe, axarr = plt.subplots(2, sharex=True,figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')

    for iten, ite in enumerate(tevec):

        (omega,specorig) = ISS1.getSpectrum(ti,ite/ite,np.log10(Ne),mi,1,0)

        axarr[0].plot(omega*1e-3,specorig/np.nanmax(specorig),linewidth=3,label=r'$T_e = {:.2f} ^\circ K$'.format(ite))

        if iten ==0:
            axarr[0].hold(True)

            axarr[0].set_ylabel('Amplitude')
            axarr[0].set_title(r'Old Spectrum, Vary $T_e$ keep $T_e/T_i=1$ ')
            axarr[0].grid(True)


        datablock = np.array([[Ne,ite,0,1,1,0],[Ne,ite,0,1,mi,0]])
        print("\n Evaluation with Te=Ti = {:.2f} K\n".format(iten))
        (omega,specoll) = ISS2.getspec(datablock)
        axarr[1].plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$T_e = {:.2f}  ^\circK$'.format(ite))
        if iten ==0:
            plt.hold(True)
            axarr[1].set_xlabel('Frequency (kHz)')

            axarr[1].set_ylabel('Amplitude')
            axarr[1].set_title(r'New Spectrums, Vary $T_e$ keep $T_e/T_i=1$ ')
            axarr[1].grid(True)


    axarr[0].legend(loc='right', shadow=True, fontsize='x-large')

#    axarr[1].legend(loc='upper right', shadow=True, fontsize='x-large')
#    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    plt.show(False)

     #%% Vary temp change ratio

    tevec = np.linspace(500.0,3000,4)
    figvarra, axarrf = plt.subplots(2, sharex=True,figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')

    for iten, ite in enumerate(tevec):

        (omega,specorig) = ISS1.getSpectrum(ti,ite/ti,np.log10(Ne),mi,1,0)

        axarrf[0].plot(omega*1e-3,specorig/np.nanmax(specorig),linewidth=3,label=r'$T_e = {:.2f} ^\circ K$'.format(ite))

        if iten ==0:
            axarrf[0].hold(True)

            axarrf[0].set_ylabel('Amplitude')
            axarrf[0].set_title(r'Old Spectrum, Vary $T_e$ keep $T_i = 1000^\circ K$ ')
            axarrf[0].grid(True)


        datablock = np.array([[Ne,ite,0,1,1,0],[Ne,ti,0,1,mi,0]])
        print("\n Evaluation with Te=Ti = {:.2f} K\n".format(iten))
        (omega,specoll) = ISS2.getspec(datablock)
        axarrf[1].plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$T_e = {:.2f}  ^\circK$'.format(ite))
        if iten ==0:
            plt.hold(True)
            axarrf[1].set_xlabel('Frequency (kHz)')

            axarrf[1].set_ylabel('Amplitude')
            axarrf[1].set_title(r'New Spectrums, Vary $T_e$ keep $T_i = 1000^\circ K$ ')
            axarrf[1].grid(True)


    axarrf[0].legend(loc='right', shadow=True, fontsize='x-large')

#    axarr[1].legend(loc='upper right', shadow=True, fontsize='x-large')
#    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    plt.show(False)


     #%% Vary temp change ratio

    tivec = np.linspace(500.0,3000,4)
    figvarra, axarrf = plt.subplots(2, sharex=True,figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')

    for itin, iti in enumerate(tevec):

        (omega,specorig) = ISS1.getSpectrum(iti,te/iti,np.log10(Ne),mi,1,0)

        axarrf[0].plot(omega*1e-3,specorig/np.nanmax(specorig),linewidth=3,label=r'$T_i = {:.2f} ^\circ K$'.format(iti))

        if itin ==0:
            axarrf[0].hold(True)

            axarrf[0].set_ylabel('Amplitude')
            axarrf[0].set_title(r'Old Spectrum, Vary $T_i$ keep $T_e = 1000^\circ K$ ')
            axarrf[0].grid(True)


        datablock = np.array([[Ne,te,0,1,1,0],[Ne,iti,0,1,mi,0]])
        print("\n Evaluation with Te=Ti = {:.2f} K\n".format(iten))
        (omega,specoll) = ISS2.getspec(datablock)
        axarrf[1].plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$T_i = {:.2f}  ^\circK$'.format(iti))
        if itin ==0:
            plt.hold(True)
            axarrf[1].set_xlabel('Frequency (kHz)')

            axarrf[1].set_ylabel('Amplitude')
            axarrf[1].set_title(r'New Spectrums, Vary $T_i$ keep $T_e = 1000^\circ K$ ')
            axarrf[1].grid(True)


    axarrf[0].legend(loc='right', shadow=True, fontsize='x-large')

#    axarr[1].legend(loc='upper right', shadow=True, fontsize='x-large')
#    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    plt.show(False)

#%% Vary percent

    ratiovec = np.linspace(0,1,5)
    figvarra, axarrf = plt.subplots(2, sharex=True,figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')

    for irotn, irot in enumerate(ratiovec):

        (omega,specorig) = ISS1.getSpectrum(ti,te/ti,np.log10(Ne),mi,1,irot)

        axarrf[0].plot(omega*1e-3,specorig/np.nanmax(specorig),linewidth=3,label=r'$H+ = {:.0f} %$'.format(irot*100))

        if irotn ==0:
            axarrf[0].hold(True)

            axarrf[0].set_ylabel('Amplitude')
            axarrf[0].set_title(r'Old Spectrum, Vary $ H^+$ ions ')
            axarrf[0].grid(True)

        datablock = np.array([[Ne,te,0,1,1,0],[(1.0-irot),ti,0,1,mi,0],[irot,ti,0,1,1.0,0]])
      # print "\n Evaluation with Te=Ti = {:.2f} K\n".format(iten)
        (omega,specoll) = ISS2.getspec(datablock)
        axarrf[1].plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$H^+ = {:.0f} %$'.format(irot*100))
        if irotn ==0:
            plt.hold(True)
            axarrf[1].set_xlabel('Frequency (kHz)')

            axarrf[1].set_ylabel('Amplitude')
            axarrf[1].set_title(r'New Spectrums, Vary $ H^+$ ions ')
            axarrf[1].grid(True)


    axarrf[0].legend(loc='right', shadow=True, fontsize='x-large')

#    axarr[1].legend(loc='upper right', shadow=True, fontsize='x-large')
#    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    plt.show(False)
    #%% set up for the rest
    Ce = np.sqrt(v_Boltz*te/v_me)
    Ci = np.sqrt(v_Boltz*ti/(v_amu*mi))

    datablock = np.array([[Ne,te,0,1,1,0],[Ne,ti,0,1,mi,0]])
    datablockn = datablock.copy()
    print("\n Evaluation using Dawson's function\n")
    (omega,specorig) = ISS2.getspec(datablock.copy())
    #%% Vary collision freq
    nuvec = np.logspace(-2.0,2.0,10)
    plt.figure()
    plt.plot(omega*1e-3,specorig/np.nanmax(specorig),marker='o', linestyle='--',linewidth=3,label="Original")
    plt.hold(True)
    plt.xlabel('Frequency (kHz)')
    plt.ylabel('Amplitude')
    plt.title('Varying Collision Frequency for O+ plasma')
    plt.grid(True)
    for inun, inu in enumerate(nuvec):
        datablockn[0,5] = inu*ISS2.K*Ce
        datablockn[1,5] = inu*ISS2.K*Ci
        print("\n Evaluation with nu = {:.2f} KC\n".format(inu))
        (omega,specoll) = ISS2.getspec(datablockn)
        plt.plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$\nu = {:.2f} KC$'.format(inu))

    plt.legend()
    plt.show(False)

    #%% vary magnetic field

    alpha = np.linspace(19,1,10)
    plt.figure()
    plt.plot(omega*1e-3,specorig/np.nanmax(specorig),marker='o', linestyle='--',linewidth=3,label="Original")
    plt.hold(True)
    plt.xlabel('Frequency (kHz)')
    plt.ylabel('Amplitude')
    plt.title('Varying Magnetic aspect angle for O+ plasma')
    plt.grid(True)
    for ialn, ial in enumerate(alpha):

        print("\n Evaluation with alpha = {:.2f} deg \n".format(ial))
        (omega,specoll) = ISS2.getspec(datablock,ial)
        plt.plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$\alpha = {:.0f}^\circ$'.format(ial))
    plt.legend()
    plt.show(False)
    #%% Adjust doppler

    dopmax = np.nanmax(omega)
    # multiply by 2 pi because the omega term is frequency in Hz not rad/s
    maxv = 2*np.pi*dopmax/ISS2.K
    vels = np.linspace(-maxv,maxv,5)
    plt.figure()
    datablockv = datablock.copy()
    for iveln, ivel in enumerate(vels):

        print("\n Evaluation with Vs = {:.2f} m/s\n".format(ivel))
        datablockv[:,2] = ivel
        (omega,specoll) = ISS2.getspec(datablockv)
        plt.plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$V_s = {:.0f}^\circ$'.format(ivel))
        if iveln==0:
            plt.hold(True)
    plt.xlabel('Frequency (kHz)')
    plt.ylabel('Amplitude')
    plt.title('Varying Doppler velocity for O+ plasma')
    plt.grid(True)
    plt.legend()
    plt.show(False)
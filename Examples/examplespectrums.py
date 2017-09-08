#!/usr/bin/env python
"""
basictest2
Created on Thu Jan 29 13:10:33 2015

@author: John Swoboda
"""


import numpy as np
import time
import ISRSpectrum.ISRSpectrum as ISSnew
from isrutilities.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge, v_me, v_amu
import matplotlib.pylab as plt

if __name__== '__main__':

    ISS2 = ISSnew.ISRSpectrum(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=129, sampfreq=100e3,dFlag=True)

    ti = 2.8e3
    te = 2.8e3
    Ne = 1e11
    mi = 16
    Ce = np.sqrt(v_Boltz*te/v_me)
    Ci = np.sqrt(v_Boltz*ti/(v_amu*mi))

    datablock = np.array([[Ne,ti,0,1,mi,0],[Ne,te,0,1,1,0]])
    datablockn = datablock.copy()
    print("\n Evaluation using Dawson's function\n")
    (omega,specorig) = ISS2.getspec(datablock.copy())
    #%% Vary collision freq
    nuvec = np.logspace(-2.0,2.0,10)
    plt.figure(facecolor='w', edgecolor='k')
    plt.plot(omega*1e-3,specorig/np.nanmax(specorig),marker='o', linestyle='--',linewidth=3,label="Original")
    plt.hold(True)
    plt.xlabel('Frequency (kHz)')
    plt.ylabel('Amplitude')
    plt.title('Varying Collision Frequency for O+ plasma')
    plt.grid(True)
    for inun, inu in enumerate(nuvec):
        datablockn[0,5] = inu*ISS2.K*Ci
        datablockn[1,5] = inu*ISS2.K*Ce
        print("\n Evaluation with nu = {:.2f} KC\n".format(inu))
        (omega,specoll) = ISS2.getspec(datablockn)
        plt.plot(omega*1e-3,specoll/np.nanmax(specoll),linewidth=3,label=r'$\nu = {:.2f} KC$'.format(inu))

    plt.legend()


    #%% vary magnetic field

    alpha = np.linspace(19,1,10)
    plt.figure(facecolor='w', edgecolor='k')
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

    #%% Adjust doppler

    dopmax = np.nanmax(omega)
    # multiply by 2 pi because the omega term is frequency in Hz not rad/s
    maxv = 2*np.pi*dopmax/ISS2.K
    vels = np.linspace(-maxv,maxv,5)
    plt.figure(facecolor='w', edgecolor='k')
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

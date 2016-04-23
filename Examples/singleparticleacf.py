#!/usr/bin/env python
"""
Created on Fri Apr 22 14:42:49 2016

@author: John Swoboda
"""

import os
import numpy as np
import scipy.fftpack as fftsy
import scipy.special
import matplotlib.pylab as plt
from matplotlib import rc
from ISRSpectrum.ISRSpectrum import magacf, collacf,magncollacf
from ISRSpectrum.const.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge, v_me, v_amu
from ISRSpectrum.const.mathutils import chirpz, sommerfeldchirpz, sommerfelderf, sommerfelderfrep


if __name__== '__main__':
    # directory stuff
    curdir = os.path.dirname(os.path.realpath(__file__))
    ISRdir = os.path.split(curdir)[0]
    figsdir = os.path.join(ISRdir,'Doc','Figs')
    #%% Sim set up
    centerFrequency = 440.2*1e6
    nspec=129
    sampfreq=50e3
    bMag = 0.4e-4
    Ts = 1e3
    Partdict = {0:('Electron',v_me,v_elemcharge),1:('H+ Ion',v_amu,v_elemcharge),2:('O+ Ion',16*v_amu,v_elemcharge)}
    particle = Partdict[0]
    pname = particle[0]

    ms = particle[1]
    q_ch = particle[2]
    K = 2.0*np.pi*2*centerFrequency/v_C_0
    f = np.arange(-np.ceil((nspec-1.0)/2.0),np.floor((nspec-1.0)/2.0+1))*(sampfreq/(2*np.ceil((nspec-1.0)/2.0)))
    C = np.sqrt(v_Boltz*Ts/ms)
    Om = q_ch*bMag/ms

    omeg_s = 2.0*np.pi*f
    theta = omeg_s/(K*C*np.sqrt(2.0))

    dtau = 2e-2*2.0/(K*C*np.sqrt(2.0))
    N = 2**12
    tau = np.arange(N)*dtau
    d2r = np.pi/180.0

#%% No collisions or magnetic field
    gordnn = np.exp(-np.power(C*K*tau,2.0)/2.0)

    plt.figure()
    plt.plot(tau,gordnn,linewidth=3)
    plt.title(r'Single Particle ACF for ' +pname)
    plt.grid(True)
    plt.show(False)
    plt.savefig(os.path.join(figsdir,'ACF'+pname.replace(" ", "")+'.png'))
#%% With collisions
    nuvec = np.logspace(-2.0,2.0,10)*K*C

    taumat = np.tile(tau[np.newaxis,:],(len(nuvec),1))
    numat = np.tile(nuvec[:,np.newaxis],(1,len(tau)))

    gordnun = collacf(taumat,K,C,numat)
    #np.exp(-np.power(K*C/numat,2.0)*(numat*taumat-1+np.exp(-numat*taumat)))

    plt.figure()
    plt.plot(tau,gordnn,linestyle='--',color='b',linewidth=4,label=r'No Collisions')
    plt.hold(True)
    
    for inun, inu in enumerate(nuvec):
        numult = inu/(K*C)
        plt.plot(tau,gordnun[inun].real,linewidth=3,label=r'$\nu = {:.2f} KC$'.format(numult))

    plt.grid(True)
    plt.title(r'Single Particle ACF W/ Collisions for '+pname)
    plt.legend()
    plt.show(False)
    plt.savefig(os.path.join(figsdir,'ACFwcolls'+pname.replace(" ", "")+'.png'))

#%% With magnetic field
    alpha = np.linspace(19,1,10)
    taumat = np.tile(tau[np.newaxis,:],(len(alpha),1))
    almat = np.tile(alpha[:,np.newaxis],(1,len(tau)))
    Kpar = np.sin(d2r*almat)*K
    Kperp = np.cos(d2r*almat)*K

#    gordmag = np.exp(-np.power(C*Kpar*taumat,2.0)/2.0-2.0*np.power(Kperp*C*np.sin(Om*taumat/2.0)/Om,2.0))
    gordmag = magacf(taumat,K,C,d2r*almat,Om)
    plt.figure()
    plt.plot(tau,gordnn,linestyle='--',color='b',linewidth=4,label='No B-field')
    plt.hold(True)
    
    for ialn, ial in enumerate(alpha):
        plt.plot(tau,gordmag[ialn].real,linewidth=3,label=r'$\alpha = {:.0f}^\circ$'.format(ial))

    plt.grid(True)
    plt.title('Single Particle ACF W/ Mag for ' +pname)
    plt.legend()
    plt.show(False)
    plt.savefig(os.path.join(figsdir,'ACFwmag'+pname.replace(" ", "")+'.png'))
    
    gordfft = np.fft()

#%% Error surface with both
    almat3d = np.tile(alpha[:,np.newaxis,np.newaxis],(1,len(nuvec),len(tau)))
    numat3d = np.tile(nuvec[np.newaxis,:,np.newaxis],(len(alpha),1,len(tau)))
    taumat3d = np.tile(tau[np.newaxis,np.newaxis,:],(len(alpha),len(nuvec),1))
    Kpar3d = np.sin(d2r*almat3d)*K
    Kperp3d = np.cos(d2r*almat3d)*K
    gam = np.arctan(numat3d/Om)

#    deltl = np.exp(-np.power(Kpar3d*C/numat3d,2.0)*(numat3d*taumat3d-1+np.exp(-numat3d*taumat3d)))
#    deltp = np.exp(-np.power(C*Kperp3d,2.0)/(Om*Om+numat3d*numat3d)*(np.cos(2*gam)+numat3d*taumat3d-np.exp(-numat3d*taumat3d)*(np.cos(Om*taumat3d-2.0*gam))))

#    gordall = deltl*deltp
    gordall=magncollacf(taumat3d,K,C,d2r*almat3d,Om,numat3d)
    gordnnmat = np.tile(gordnn[np.newaxis,np.newaxis,:],(len(alpha),len(nuvec),1))

    gorddiff = np.abs(gordall-gordnnmat)**2
    err = np.sqrt(gorddiff.sum(2))/np.sqrt(np.power(gordnn,2.0).sum())
    extent = [np.log10(nuvec[0]/(K*C)),np.log10(nuvec[-1]/(K*C)),alpha[0],alpha[-1]]
    plt.figure()
    myim = plt.imshow(err*100,extent = extent,origin='lower',aspect='auto')
    myim.set_clim(0.0,5)
    plt.xlabel(r'$\log_{10}(\nu /KC)$')
    plt.ylabel(r'$^\circ\alpha$')
    cbar = plt.colorbar()
    cbar.set_label('% Error', rotation=270)
    cbar.ax.get_yaxis().labelpad = 15
    plt.title('Error between ideal ACF and with Collisions and B-field for '+pname)
    plt.savefig(os.path.join(figsdir,'ACFerr'+pname.replace(" ", "")+'.png'))
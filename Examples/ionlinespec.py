#!/usr/bin/env python
"""
Created on Mon Aug 29 16:29:42 2016

@author: John Swoboda
"""
import numpy as np
import scipy.fftpack as scfft
import os,inspect
from ISRSpectrum.ISRSpectrum import ISRSpectrum
import matplotlib.pylab as plt
import seaborn as sns
import pdb
if __name__== '__main__':
    sns.set_style("white")
    sns.set_context("notebook")
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    imagepath = os.path.join(os.path.split(curpath)[0],'Doc','Figs')
    databloc = np.array([[1e11,1e3],[1e11,2.5e3]])
    nspec=256
    spfreq=50e3
    ISpec_ion = ISRSpectrum(centerFrequency = 449e6, nspec=nspec, sampfreq=spfreq,dFlag=True)
    species=['O+','e-']
#    databloc = np.array([[1.66e10,863.],[1.66e10,863.]])
    ylim=[0,1.4]
    ylim2=[-.5,1.4]
    flims=np.array([3.03e6,3.034e6])
    #%% With B-Field
    
    
    fion,ionline= ISpec_ion.getspecsep(databloc,species)
    
    acf=scfft.ifft(scfft.ifftshift(ionline)).real
    tau=scfft.ifftshift(np.arange(-np.ceil((float(nspec)-1)/2),np.floor((float(nspec)-1)/2)+1))/spfreq
    
    
    
    
    fig,ax = plt.subplots(1,1,sharey=True, figsize=(6,4),facecolor='w')
    
    l1=ax.plot(fion*1e-3,ionline/ionline.max(),'-',lw=3)[0]
    sns.despine()

    ax.set_xlim([-15,15])
    ax.spines['right'].set_visible(False)

    ax.set_xlabel('Frequency (kHz)',fontsize=14)
    ax.set_title(r'$\langle|n_e(|\mathbf{k}|=18.5,\omega)|^2\rangle$',fontsize=18)
    ax.set_ylabel(r'Normalized Magnitude',fontsize=14)
    plt.tight_layout()
    ax.set_ylim(ylim)
    plt.savefig(os.path.join(imagepath,'Specion.png'),dpi=300)
    
    
    fig,ax = plt.subplots(1,1,sharey=True, figsize=(6,4),facecolor='w')
    
    l1=ax.plot(tau[:64]*1e6,acf[:64]/acf[0],'-',lw=3)[0]
    sns.despine()

    ax.set_xlim([0,280])
    ax.spines['right'].set_visible(False)

    ax.set_xlabel(r'$\tau$ in $\mu$s ',fontsize=14)
    ax.set_title(r'$\langle|n_e(|\mathbf{k}|=18.5,\tau)|^2\rangle$',fontsize=18)
    ax.set_ylabel(r'Normalized Magnitude',fontsize=14)
    plt.tight_layout()
    ax.set_ylim(ylim2)
    plt.savefig(os.path.join(imagepath,'acfion.png'),dpi=300)
    #%% With random values
    n_pulse=np.array([1,10,50,200])
    lab_strs=['J={0}'.format(int(i))for i in n_pulse]
    lab_strs.insert(0,'Original')
    np_1=n_pulse.max()
    
    x=(np.random.randn(np_1,nspec)+1j*np.random.randn(np_1,nspec))/np.sqrt(2.)
    filt=np.tile(np.sqrt(ionline[np.newaxis]),(np_1,1)).astype(x.dtype)
    y=x*filt
    ysqrt=y.real**2+y.imag**2
    ysqrtfft=scfft.ifft(scfft.ifftshift(ysqrt,axes=-1),axis=-1)
    ystats=np.array([ysqrt[:i].mean(axis=0) for i in n_pulse])
    #ystatsacf=scfft.ifft(scfft.ifftshift(ystats,axes=-1),axis=-1).real
    ystatsacf=np.array([ysqrtfft[:i].mean(axis=0).real for i in n_pulse])
    fig,ax = plt.subplots(1,1,sharey=True, figsize=(6,4),facecolor='w')
    
    l1=ax.plot(fion*1e-3,ionline/ionline.max(),'-',lw=3,zorder=len(n_pulse))[0]
    sns.despine()

    ax.set_xlim([-15,15])
    ax.spines['right'].set_visible(False)
    hand=[l1]    
    for iyn,iy in enumerate(ystats):
        l1=ax.plot(fion*1e-3,iy/ionline.max(),'-',lw=3,zorder=iyn)[0]
        hand.append(l1)
    ax.set_xlabel('Frequency (kHz)',fontsize=14)
    ax.set_title(r'Ion Line Averaging',fontsize=18)
    ax.set_ylabel(r'Normalized Magnitude',fontsize=14)
    ax.set_ylim(ylim)
    plt.tight_layout()
    ax.legend(hand,lab_strs,)
    plt.savefig(os.path.join(imagepath,'Specionave.png'),dpi=300)
    
    
    
    fig,ax = plt.subplots(1,1,sharey=True, figsize=(6,4),facecolor='w')
    l1=ax.plot(tau[:64]*1e6,acf[:64]/acf[0],'-',lw=3,zorder=len(n_pulse))[0]
    hand=[l1]
    for iyn,iy in enumerate(ystatsacf):
        l1=ax.plot(tau[:64]*1e6,iy[:64]/acf[0],'-',lw=3,zorder=iyn)[0]
        hand.append(l1)
    sns.despine()

    ax.set_xlim([0,280])
    ax.spines['right'].set_visible(False)

    ax.set_xlabel(r'$\tau$ in $\mu$s ',fontsize=14)
    ax.set_title(r'ACF Averaging',fontsize=18)
    ax.set_ylabel(r'Normalized Magnitude',fontsize=14)
    plt.tight_layout()
    ax.set_ylim(ylim2)
    ax.legend(hand,lab_strs,)
    plt.savefig(os.path.join(imagepath,'acfionave.png'),dpi=300)
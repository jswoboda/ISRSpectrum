#!/usr/bin/env python
"""
Created on Mon Aug 29 16:29:42 2016

@author: John Swoboda
"""
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style("white")
sns.set_context("notebook")
#
from ISRSpectrum.ISRSpectrum import ISRSpectrum


if __name__== '__main__':

    databloc = np.array([[1e11,1e3],[1e11,2.5e3]])
    nspec=256
    ISpec_ion = ISRSpectrum(centerFrequency = 449e6, nspec=nspec, sampfreq=50e3,dFlag=True)
    species=['O+','e-']
#    databloc = np.array([[1.66e10,863.],[1.66e10,863.]])
    
    
    flims=np.array([3.03e6,3.034e6])
    #%% With B-Field
    
    
    fion,ionline= ISpec_ion.getspecsep(databloc,species)
    
    
    
    
    fig,ax = plt.subplots(1,1,sharey=True, figsize=(6,4),facecolor='w')
    
    l1=ax.plot(fion*1e-3,ionline/ionline.max(),'-',lw=3)[0]
    sns.despine()

    ax.set_xlim([-15,15])
    ax.spines['right'].set_visible(False)

    ax.set_xlabel('Frequency (kHz)',fontsize=14)
    ax.set_title(r'$\langle|n_e(|\mathbf{k}|=18.5,\omega)|^2\rangle$',fontsize=18)
    ax.set_ylabel(r'Normalized Magnitude',fontsize=14)
    plt.tight_layout()

    plt.savefig('Specion.png',dpi=300)
    
    #%% With random values
    n_pulse=np.array([50,200,500,1000])
    lab_strs=['J={0}'.format(int(i))for i in n_pulse]
    lab_strs.insert(0,'Original')
    np_1=n_pulse.max()
    
    x=(np.random.randn(np_1,nspec)+1j*np.random.randn(np_1,nspec))/np.sqrt(2.)
    filt=np.tile(np.sqrt(ionline[np.newaxis]),(np_1,1)).astype(x.dtype)
    y=x*filt
    ysqrt=y.real**2+y.imag**2
    
    ystats=np.array([ysqrt[:i].mean(axis=0) for i in n_pulse])
    
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
    plt.tight_layout()
    ax.legend(hand,lab_strs,)
    plt.savefig('Specionave.png',dpi=300)

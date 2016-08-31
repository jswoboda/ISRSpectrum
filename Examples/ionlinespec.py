#!/usr/bin/env python
"""
Created on Mon Aug 29 16:29:42 2016

@author: John Swoboda
"""
import numpy as np
import os,inspect
from ISRSpectrum.ISRSpectrum import ISRSpectrum
import matplotlib.pylab as plt
import seaborn as sns
if __name__== '__main__':
    sns.set_style("white")
    sns.set_context("notebook")
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    imagepath = os.path.join(os.path.split(curpath)[0],'Doc','Figs')
    databloc = np.array([[1e11,1e3],[1e11,2.5e3]])
   
    ISpec_ion = ISRSpectrum(centerFrequency = 449e6, nspec=256, sampfreq=50e3,dFlag=True)
    species=['O+','e-']
#    databloc = np.array([[1.66e10,863.],[1.66e10,863.]])
    
    
    flims=np.array([3.03e6,3.034e6])
    #%% With B-Field
    
    
    fion,ionline= ISpec_ion.getspecsep(databloc,species)
    
    
    
    
    fig,ax = plt.subplots(1,1,sharey=True, figsize=(4,4),facecolor='w')
    
    l1=ax.plot(fion*1e-3,ionline/ionline.max(),'-',lw=3)[0]
    sns.despine()

    ax.set_xlim([-25,25])
    ax.spines['right'].set_visible(False)

    ax.set_xlabel('Frequency (kHz)',fontsize=14)
    ax.set_title(r'$\langle|n_e(|\mathbf{k}|=18.5,\omega)|^2\rangle$',fontsize=18)
    ax.set_ylabel(r'Normalized Magnitude',fontsize=14)
    plt.tight_layout()

    plt.savefig(os.path.join(imagepath,'Specion.png'),dpi=300)
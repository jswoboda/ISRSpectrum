#!/usr/bin/env python
"""
This example creates a 2-d mesh plot varying the radar wavelength to show how this impacts the IS spectrum. Every spectrum is normalized to it largest entry. Plots two figures one in radar frequency and the other in k.
"""

import numpy as np
import scipy.constants as spconst
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style("white")
sns.set_context("notebook")
#
from ISRSpectrum import Specinit

def spec2dplot():


    databloc = np.array([[1e11,1e3],[1e11,2.5e3]])
    species=['O+','e-']
    newcentfreq=    449e6 +np.linspace(-200,550,250)*1e6
    k_all= 2*2*np.pi*newcentfreq/spconst.c
    k_lims=np.array([k_all.min(),k_all.max()])
    freq_lims=np.array([newcentfreq.min(),newcentfreq.max()])
    oution = []
    for i,if_0 in enumerate(newcentfreq):
        ISpec_ion = Specinit(centerFrequency = if_0, nspec=256, sampfreq=50e3,dFlag=False)
        fion,ionline= ISpec_ion.getspecsep(databloc,species)
        # Normalize each spectrum individually
        oution.append(ionline/ionline.max())
    oution=np.array(oution)


    F,K_b=np.meshgrid(fion,k_all)

    fig,ax = plt.subplots(1,1,sharey=True, figsize=(4,4),facecolor='w')

    l1=ax.pcolor(F*1e-3,K_b,oution,cmap='viridis')

    cb1 = plt.colorbar(l1, ax=ax)
    ax.set_xlim([-25,25])
    ax.set_ylim(k_lims)
    ax.set_xlabel('Frequency (kHz)',fontsize=14)
    ax.set_title(r'$\langle|n_e(\mathbf{k},\omega)|^2\rangle$',fontsize=18)
    ax.set_ylabel(r'$|\mathbf{k}|$ (rad/m)',fontsize=14)
    plt.tight_layout()



    plt.savefig('Spec2dion.png',dpi=300)

    fig,ax = plt.subplots(1,1,sharey=True, figsize=(4,4),facecolor='w')
    l1=ax.pcolor(F*1e-3,1e-9*K_b*spconst.c/(2*2*np.pi),oution,cmap='viridis')

    cb1 = plt.colorbar(l1, ax=ax)
    ax.set_xlim([-25,25])
    ax.set_ylim(freq_lims*1e-9)
    ax.set_xlabel('Frequency (kHz)',fontsize=14)
    ax.set_title('Ion Line',fontsize=18)
    ax.set_ylabel(r'Radar Center Frequency (GHz)',fontsize=14)
    plt.tight_layout()
    plt.savefig('Spec2dionf0.png',dpi=300)

if __name__ == "__main__":
    spec2dplot()
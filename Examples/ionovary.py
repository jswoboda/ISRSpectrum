#!/usr/bin/env python
"""
elineexample.py
Created on Sun Dec 27 15:28:19 2015
This example shows everything up to electron line for magnitized and non-magnitized plasmas.
@author: John Swoboda
"""

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("notebook")
#
from ISRSpectrum import Specinit

if __name__== '__main__':


    ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=129, sampfreq=50e3,dFlag=True)

    ti_list = np.linspace(1000,5000,5)
    te = 2e3
    Ne = 1e11
    Ni = 1e11
    species=['O+','e-']
    databloc = np.array([[1e11,1100.],[1e11,2500.]])

    fig,ax = plt.subplots()
    for ti in ti_list:

        datablock = np.array([[Ni,ti],[Ne,te]])
        species = ['O+','e-']
        (f,spec1,rcs) = ISS2.getspecsep(datablock,species,vel = 0.0, alphadeg=90.0,rcsflag=True)
        ax.plot(f*1e-3,spec1/np.nanmax(spec1),marker='o', linestyle='--',linewidth=3,label="{0} K".format(int(ti)))
    ax.legend()
    plt.savefig('Iontempvary.png')

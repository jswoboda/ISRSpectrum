#!/usr/bin/env python
"""
Created on Thu Nov 12 10:51:26 2015

@author: John Swoboda
"""
from pathlib import Path
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
#
from ISRSpectrum import Specinit


def main():
    sns.set_style("whitegrid")
    sns.set_context("notebook")
    n_species=['H', 'He','N','N2', 'O', 'O2']
    Nn = np.array([22340974.0, 1.170203e+08, 269365.12, 8.7995607e+12,5.6790339e+11,1.9389578e+12])*1e6
    Tn = 179.29138
    n_datablock = np.zeros((len(Nn),2))
    n_datablock[:,0] = Nn
    n_datablock[:,1] = Tn


    species = ['NO+','O2+','e-']
    Ni = np.array([1998.2988800000001, 57.700355999999999])*1e6
    Ti = 175.1
    Ne= 2055.9992320000001 *1e6
    Te = 175.1
    datablock = np.zeros((len(Ni)+1,2))
    datablock[:-1,0] = Ni
    datablock[:-1,1] = Ti
    datablock[-1,0] = Ne
    datablock[-1,1] = Te

    #
    ISS1 = Specinit(centerFrequency = 449e6, bMag = 0.4e-4, nspec=512, sampfreq=25e3,dFlag=True)

    (omeg,spec_1)=ISS1.getspecsep(datablock,species,vel = 0.0,rcsflag=False)
    (omeg,spec_2)=ISS1.getspecsep(datablock,species,vel = 0.0,rcsflag=False,col_calc = True)
    (omeg,spec_3)=ISS1.getspecsep(datablock,species,vel = 0.0,rcsflag=False,col_calc = True,n_datablock=n_datablock,n_species=n_species)


    (figmplf, axmat) = plt.subplots(1, 1,figsize=(20, 15), facecolor='w')
    curax = axmat
    lines = [None]*3
    labels = ['Spectrum','Spectrum with Coulomb Collisions','Spectrum with Columb and Neutral Collisions']
    lines[0] = curax.plot(omeg*1e-3,spec_1/np.max(spec_1),label='Output',linewidth=5)[0]
    lines[1] = curax.plot(omeg*1e-3,spec_2/np.max(spec_2),label='Output',linewidth=5)[0]
    lines[2] = curax.plot(omeg*1e-3,spec_3/np.max(spec_3),label='Output',linewidth=5)[0]

    figmplf.suptitle('Spectrums with and without collisions')
    plt.figlegend( lines, labels, loc = 'lower center', ncol=5, labelspacing=0. )

    plt.savefig('Collisions.png')

if __name__== '__main__':

    main()

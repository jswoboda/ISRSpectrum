#!/usr/bin/env python
"""
This script will create plots of ISR spectra for different radar systems

@author: John Swoboda
"""

import scipy as sp
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("notebook")
#
import ISRSpectrum.ISRSpectrum as ISSnew
from isrutilities.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge, v_me, v_amu


def main():

    Ne = 1e11
    Ti = 3e3
    Te = 3e3
    databloc =  sp.array([[Ne,Ti],[Ne,Te]])
    species = ['O+','e-']

    Cia = sp.sqrt(v_Boltz*(Te+3.*Ti)/16./v_amu)

    #make list of dictionaries
    biglist = [{'name':'AMISR','Fo':449e6,'Fs':50e3,'alpha':70.},
        {'name':'Sondrestrom','Fo':1290e6,'Fs':100e3,'alpha':80.},
        {'name':'Haystack','Fo':440e6,'Fs':50e3,'alpha':65.},
        {'name':'Arecibo','Fo':430e6,'Fs':50e3,'alpha':45.},
        {'name':'Jicamarca','Fo':50e6,'Fs':10e3,'alpha':1.}]

    (figmplf, axmat) = plt.subplots(3, 2,figsize=(20, 15), facecolor='w')
    axvec = axmat.flatten()
    lines = [None]*2
    labels = ['Spectrum','Ion Acoustic Frequency']
    for ima,idict in enumerate(biglist):
        curax = axvec[ima]
        k = 2.0*sp.pi*2*idict['Fo']/v_C_0
        xloc = sp.array([-k*Cia,k*Cia])/2/sp.pi/2
        ISS1 = ISSnew.ISRSpectrum(centerFrequency = idict['Fo'], bMag = 0.4e-4, nspec=256, sampfreq=idict['Fs'],dFlag=True)
        (omeg,spec)=ISS1.getspecsep(databloc,species,vel = 0.0, alphadeg=idict['alpha'],rcsflag=False)

        lines[0] = curax.plot(omeg*1e-3,spec,label='Output',linewidth=5)[0]
        lines[1] = curax.stem(xloc*1e-3, sp.ones(2)*sp.amax(spec), linefmt='g--', markerfmt='go', basefmt=' ')[0]
        curax.set_xlabel('f in kHz')
        curax.set_ylabel('Amp')
        curax.set_title(idict['name']+ ' Spectrtum')
    figmplf.suptitle(r'Spectrums $N_e$ = {0:.1e}m$^{{-3}}$, $T_e$ = {1:.0f}$^o$K, $T_i$ = {2:.0f}$^o$K'.format(Ne,Te,Ti), fontsize=20)
    plt.figlegend( lines, labels, loc = 'lower center', ncol=5, labelspacing=0. )

    plt.savefig('DifferentSystems.png')
if __name__== '__main__':

    main()

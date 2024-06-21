#!/usr/bin/env python
"""
Created on Thu Aug 25 21:43:45 2016

@author: John Swoboda
"""
import numpy as np
import scipy as sp
import scipy.constants as spconst
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("notebook")

from ISRSpectrum import Specinit

def main():
    Ne = 1e11
    Ti = 1.1e3
    Te = 3e3
    databloc =  np.array([[Ne,Ti],[Ne,Te]])
    species = ['O+','e-']

    mult=np.arange(1,4)
    mult[-1]=50
    Cia = np.sqrt(spconst.k*(2.*Ti)/16./spconst.m_p)

    #make list of dictionaries
    dict1={'name':'AMISR','Fo':449e6,'Fs':50e3,'alpha':70.}

    ISS1 = Specinit(centerFrequency = dict1['Fo'], bMag = 0.4e-4, nspec=256, sampfreq=dict1['Fs'],dFlag=True)
    (figmplf, ax) = plt.subplots(2, 1,figsize=(8, 6), facecolor='w')
    lines = []
    labels = ['Spectrum','Ion Acoustic Frequency']
    for ima,imult in enumerate(mult):
        k = 2*dict1['Fo']/spconst.c
        databloc[1,1]=Ti*imult
        Cia = np.sqrt(spconst.k*(imult*Ti+Ti)/(16.*spconst.m_p))
        xloc = np.array([-k*Cia,k*Cia])

        (omeg,spec)=ISS1.getspecsep(databloc,species,vel = 0.0, alphadeg=dict1['alpha'],rcsflag=False)

        if ima==0:

            ax[0].set_xlabel('f in kHz')
            ax[0].set_ylabel('Amp')
            ax[0].set_title('Spectra')
            ax[1].set_xlim([0,480])
            ax[1].spines['right'].set_visible(False)

            ax[1].set_xlabel('Lag in Microseconds ',fontsize=14)
            ax[1].set_title('ACF K={:0.1f}'.format(kvec),fontsize=18)
            ax[1].set_ylabel('Normalized Magnitude',fontsize=14)
        
        l1=ax[1].plot(tau[:64]*1e6,acf[:64]/acf[0],'-',lw=3)[0]
        sns.despine()


        lines.append( ax[0].stem(xloc*1e-3, np.ones(2)*np.amax(spec), linefmt='g--', markerfmt='go', basefmt=' ')[0])
        lines.append( ax[0].plot(omeg*1e-3,spec,label='Output',linewidth=5)[0])

    plt.savefig('DifferentTemps.png')
if __name__== '__main__':

    main()

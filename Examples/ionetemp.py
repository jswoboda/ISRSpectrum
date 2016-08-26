#!/usr/bin/env python
"""
Created on Thu Aug 25 21:43:45 2016

@author: John Swoboda
"""

import scipy as sp
import inspect
import os
import ISRSpectrum.ISRSpectrum as ISSnew
from ISRSpectrum.const.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge, v_me, v_amu
import matplotlib.pylab as plt
import seaborn as sns

def main():
    sns.set_style("whitegrid")
    sns.set_context("notebook")
    Ne = 1e11
    Ti = 1.1e3
    Te = 3e3
    databloc =  sp.array([[Ne,Ti],[Ne,Te]])
    species = ['O+','e-']
    
    mult=sp.arange(1,4)
    mult[-1]=50
    Cia = sp.sqrt(v_Boltz*(2.*Ti)/16./v_amu)

    #make list of dictionaries
    dict1={'name':'AMISR','Fo':449e6,'Fs':50e3,'alpha':70.}
    
    ISS1 = ISSnew.ISRSpectrum(centerFrequency = dict1['Fo'], bMag = 0.4e-4, nspec=256, sampfreq=dict1['Fs'],dFlag=True)
    (figmplf, axmat) = plt.subplots(1, 1,figsize=(8, 6), facecolor='w')
    lines = []
    labels = ['Spectrum','Ion Acoustic Frequency']
    for ima,imult in enumerate(mult):
        k = 2*dict1['Fo']/v_C_0
        databloc[1,1]=Ti*imult
        Cia = sp.sqrt(v_Boltz*(imult*Ti+Ti)/(16.*v_amu))
        xloc = sp.array([-k*Cia,k*Cia])
        
        (omeg,spec)=ISS1.getspecsep(databloc,species,vel = 0.0, alphadeg=dict1['alpha'],rcsflag=False)

        if ima==0:
            
            axmat.set_xlabel('f in kHz')
            axmat.set_ylabel('Amp')
            axmat.set_title('Spectrtum')
        lines.append( axmat.stem(xloc*1e-3, sp.ones(2)*sp.amax(spec), linefmt='g--', markerfmt='go', basefmt=' ')[0])
        lines.append( axmat.plot(omeg*1e-3,spec,label='Output',linewidth=5)[0])
    

    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    imagepath = os.path.join(os.path.split(curpath)[0],'Doc','Figs')
    plt.savefig(os.path.join(imagepath,'DifferentTemps.png'))
if __name__== '__main__':

    main()
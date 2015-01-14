#!/usr/bin/env python
"""
Created on Tue Jul 15 16:12:05 2014

@author: John Swoboda
"""
import numpy as np
import scipy as sp
import scipy.special
import time
import pdb
from const.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge

class ISRSpectrum(object):
    def __init__(self,centerFrequency = 440.2, bMag = 0.4e-4, nspec=64, sampfreq=50e3):
        self.bMag = bMag
        self.K = 2*np.pi*v_C_0/centerFrequency
        self.f = np.arange(-np.ceil((nspec-1.0)/2.0),np.floor((nspec-1.0)/2.0+1))*(sampfreq/(2*np.ceil((nspec-1.0)/2.0)))
        self.omeg = 2.0*np.pi*self.f

    def getspec(self,datablock,alphadeg=90.0):

        alpha = alphadeg*np.pi/180
        estuff = datablock[0]
        ionstuff = datablock[1:]
        (egord,h_e,Ne,omeg_e) = self.__calcoordeyev__(estuff,alpha)

        sig_e = 1j*v_epsilon0*(1-omeg_e*egord)/(self.K**2*h_e**2)
        nte = 2*N_e*np.real(egord)

        firstion = True
        for iinfo in ionostuff:
            (igord,h_i,Ni,omeg_i) = self.__calcoordeyev__(iinfo,alpha)

            sig_i = 1j*v_epsilon0*(1-omeg_i*igord)/(self.K**2*h_i**2)
            nti = 2*N_i*np.real(igord)
            if firstion:
                sig_sum =sig_i
                nt_sum = np.abs(nti)**2
            else:
                sig_sum = sig_i+sig_sum
                nt_sum = np.abs(nti)**2+nt_sum

        num = np.abs(sig_e)**2*nt_sum +np.abs(1j*self.omeg*v_epsilon0+sig_sum)**2*np.abs(nte)**2
        den = np.abs(1j*self.omeg)
        spec = num/den
        return (self.f,spec)

    def __calcgordeyev__(self,dataline,alpha):
        aldiff = np.abs(alpha*180.0/np.pi-90.0)<1;
        K = self.k

        (Ns,Ts,Vs,qs,ms,nus) = dataline[:7]
        omeg_s = self.omeg - K*Vs
        hs = np.sqrt(v_epsilon0*v_Boltz*Ts/(Ns*qs*v_elemcharge))
        C = np.sqrt(v_Boltz*Ts/ms)


        theta = np.sqrt(2.0)*v_Boltz;
        jomeg = spsf.dawsn(omeg)

        if K*C>10.0*nus and aldiff:
            #for case with no collisions and magnetic field
            gord = scipy.special.dawsn(omeg_s/(K*C*np.sqrt(2)))/(K*C*np.sqrt(2))


        return (gord,hs,Ns,omeg_s)
if __name__== '__main__':


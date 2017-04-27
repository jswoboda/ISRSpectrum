#!/usr/bin/env python
"""
Created on Tue Jul 15 16:12:05 2014

@author: John Swoboda
ISRSpectrum.py
This module will create ISR spectrums using plasma parameters. The code implements the
method of calculation found in the following article.

E. Kudeki and M. A. Milla, 2011.

The intent of the code is to be able to calculate an ISR spectrum in a number of
different conditions except for a very low magnetic aspect angles (<1deg).
"""
from __future__ import absolute_import
from ISRSpectrum import Path
from six import string_types
import scipy as sp
import scipy.special
import pandas as pd
#
from isrutilities.physConstants import v_Boltz, v_C_0, v_epsilon0, v_elemcharge, v_me, v_amu
from isrutilities.mathutils import sommerfelderfrep

INFODICT = {'O+':sp.array([1,16]),'NO+':sp.array([1,30]),
                'N2+':sp.array([1,28]),'O2+':sp.array([1,32]),
                'N+':sp.array([1,14]), 'H+':sp.array([1,1]),
                'e-':sp.array([-1,1])}
class ISRSpectrum(object):
    """ Class to create the spectrum. The instance of the class will hold infomation on
    the radar system such as sample frequency, center frequency and number of points for
    the spectrum.
    Parameters
    bMag - The magnetic field magnitude in Teslas.
    K - The Bragg scattering vector magnitude corresponds to 1/2 radar wavelength.
    f - Vector holding the frequency values in Hz
    omeg - Vector hold the frequency values in rad/s
    collfreqmin - Minimum collision frequency before they are taken into account in the Gordeyev integral calculation
    alphamax - The maximum aspect angle
    dFlag - A debug flag."""
    def __init__(self,centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=64, sampfreq=50e3,collfreqmin=1e-2,alphamax=30.0,dFlag=False,f=None):
        """ Constructor for the class.
        Inputs :
        centerFrequency: The radar center frequency in Hz.
        bMag: The magnetic field magnitude in Teslas.
        nspec: the number of points of the spectrum.
        sampfreq: The sampling frequency of the A/Ds in Hz
        collfreqmin: (Default 1e-2) The minimum collision frequency needed to incorporate it into Gordeyev
            integral calculations in units of K*sqrt(Kb*Ts/ms) for each ion species.
        alphamax: (Default 30)  The maximum magnetic aspect angle in which the B-field will be taken into account.
        dFlag: A debug flag, if set true will output debug text. Default is false.
        f: A numpy array of frequeny points, in Hz, the spectrum will be formed over. Default is
            None, at that point the frequency vector will be formed using the number of points for the spectrum
            and the sampling frequency to create a linearly sampled frequency vector. """
        self.bMag = bMag
        self.dFlag = dFlag
        self.collfreqmin = collfreqmin
        self.alphamax = alphamax

        self.K = 2.0*sp.pi*2*centerFrequency/v_C_0 #The Bragg scattering vector, corresponds to half the radar wavelength.
        if f is None:
            minfreq = -sp.ceil((nspec-1.0)/2.0)
            maxfreq = sp.floor((nspec-1.0)/2.0+1)
            self.f = sp.arange(minfreq,maxfreq)*(sampfreq/(2*sp.ceil((nspec-1.0)/2.0)))
        else:
            self.f=f
        self.omeg = 2.0*sp.pi*self.f

    def getspec(self,datablock,alphadeg=90.0,rcsflag=False,seplines=False):
        """ Gives the spectrum and the frequency vectors given the block of data and
        the magnetic aspect angle.
        Inputs
        datablock: A numpy array of size 1+Nionsx6 that holds the plasma parameters needed
        to create the spectrum. The last row will hold the information for the electrons.
        Each row of the array will have the following set up.
            [Ns, Ts, Vs, qs, ms, nus]
            Ns - The density of the species in m^-3
            Ts - Tempretur of the species in degrees K
            Vs - The Doppler velocity in m/s.
            qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
            ms - Mass of the species in AMU. (Value will be replaced for the electrons)
            nus - Collision frequency for species in s^-1.
        alphadeg: The magnetic aspect angle in degrees.
        rcsflag: A bool that will determine if the reflected power is returned as well. (default is False)
        seplines: A bool that will change the output, spec to a list of numpy arrays that include
            return from the electron line and ion lines seperatly.
        Outputs
        f - A numpy array that holds the frequency vector associated with the spectrum,
            in Hz.
        spec - The resulting ISR spectrum from the parameters with the max value set to 1.
        rcs - The RCS from the parcle of plasma for the given parameters. The RCS is
            in m.
        """
        # perform a copy of the object to avoid it being written over with incorrect.
        datablock=datablock.copy()
        dFlag = self.dFlag
        alpha = alphadeg*sp.pi/180
        estuff = datablock[-1]
        Nions = datablock.shape[0]-1
        estuff[3] = -v_elemcharge
        estuff[4] = v_me
        ionstuff = datablock[:-1]
        if dFlag:
            print("Calculating Gordeyev int for electons")
        (egord,Te,Ne,omeg_e) = self.__calcgordeyev__(estuff,alpha)
        h_e = sp.sqrt(v_epsilon0*v_Boltz*Te/(Ne*v_elemcharge*v_elemcharge))
        sig_e = (1j+omeg_e*egord)/(self.K**2*h_e**2)
        nte = 2*Ne*sp.real(egord)

        #adjust ion stuff
        ionstuff[:,3] = ionstuff[:,3]*v_elemcharge
        ionstuff[:,4] = ionstuff[:,4]*v_amu
        ionden = sp.sum(ionstuff[:,0])
        # normalize total ion density to be the same as electron density
        ionstuff[:,0] = (estuff[0]/ionden)*ionstuff[:,0]

        # ratio of charges between ion species and electrons
        qrotvec = ionstuff[:,3]/estuff[3]
        firstion = True
        wevec = sp.zeros(Nions)
        Tivec = sp.zeros(Nions)
        for iion,iinfo in enumerate(ionstuff):
            if dFlag:
                print("Calculating Gordeyev int for ion species #{:d}".format(iion))

            (igord,Ti,Ni,omeg_i) = self.__calcgordeyev__(iinfo,alpha)


            too_big = self.f>1e5
#            igord[too_big] = 0.

            wevec[iion] = Ni/Ne
            Tivec[iion] = Ti
            # sub out ion debye length because zero density of ion species can cause a divid by zero error.
            # mu=Ti/Te, tempreture ratio
            mu = Ti/Te
            qrot = qrotvec[iion]
            sig_i = (Ni/Ne)*(1j+omeg_i*igord)/(self.K**2*mu*h_e**2/qrot**2)
            nti = 2*Ni*sp.real(igord)
#            nti[too_big] = 0.
            if firstion:
                sig_sum =sig_i
                nt_sum = nti
                firstion=False
            else:
                sig_sum = sig_i+sig_sum
                nt_sum = nti+nt_sum

        inum = sp.absolute(sig_e)**2*nt_sum
        enum = sp.absolute(1j+sig_sum)**2*nte
        den = sp.absolute(1j + sig_sum +sig_e)**2
        iline = inum/den
        eline = enum/den
        if seplines:
            spec=[iline,eline]
        else:
            spec = iline+eline
        if rcsflag:
            Tr = Te/sp.sum(wevec*Tivec)
            rcs = Ne/((1+self.K**2*h_e**2)*(1+self.K**2*h_e**2+Tr))

            return (self.f,spec,rcs)
        else:
            return (self.f,spec)
    def getspecsep(self,datablock,species,vel = 0.0, alphadeg=90.0,rcsflag=False,col_calc = False,n_datablock=None,n_species=None,seplines=False):
        """ This function is a different way of getting the spectrums. A datablock is
        still used, (Nsp x 2) numpy array, but it is filled in by using the species
        listed as a string in the list speces.
        inputs
        datablock - A numpy array of size Nsp x2. The first element of the row is the
            density of species in m^-3 and the second element is the Tempreture in
            degrees K.
        species - A Nsp list of strings that label each species.
        vel - The line of site velocity of the ions in m/s, default = 0.0
        alphadeg - The angle offset from the magnetic field in degrees, default is 90.0
        rcsflag - If this flag is True then a third output of the rcs will be made.
            Default value is False
        col_calc - If this flag is true then collisions will be calculated. (Default= False)
        n_datablock - A numpy array of size Nnsp x2. The first element of the row is
            the density of the neutral species in m^-3 and the second element is the
            Tempreture in degrees K.If set to None then (Default=None)
        n_species - A Nnsp list of strings that label each neutral species.(Default=None)
        seplines: A bool that will change the output, spec to a list of numpy arrays that include
            return from the electron line and ion lines seperatly.
        Outputs
        f - A numpy array that holds the frequency vector associated with the spectrum,
            in Hz.
        spec - The resulting ISR spectrum from the parameters with the max value set to 1.
        rcs - The RCS from the parcle of plasma for the given parameters. The RCS is
            in m.

        """
        assert Islistofstr(species),"Species needs to be a list of strings"
        assert allin(species,list(INFODICT.keys())), "Have un named species in the list."
        nspec = datablock.shape[0]
        datablocknew = sp.zeros((nspec,7))

        (nuparr,nuperp) = get_collisionfreqs(datablock,species,n_datablock, n_species)
        for nspec, ispec in enumerate(species):
            datablocknew[nspec,:2] = datablock[nspec]
            datablocknew[nspec,2] = vel
            datablocknew[nspec,3:5] = INFODICT[ispec]
            if col_calc:
                datablocknew[nspec,5] =nuparr[nspec]
                datablocknew[nspec,6] =nuparr[nspec]

        return self.getspec(datablocknew,alphadeg,rcsflag,seplines=seplines)

    def __calcgordeyev__(self,dataline,alpha,alphdiff = 10.0):
        """ Performs the Gordeyve integral calculation.
        Inputs
        dataline: A numpy array of length that holds the plasma parameters needed
        to create the spectrum.
        Each row of the array will have the following set up.
            [Ns, Ts, Vs, qs, ms, nus]
            Ns - The density of the species in m^-3
            Ts - Tempretur of the species in degrees K
            Vs - The Doppler velocity in m/s.
            qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
            ms - Mass of the species in AMU. (Value will be replaced for the electrons)
            nus - Collision frequency for species in s^-1.
        alphadeg: The magnetic aspect angle in radians.
        Output
        gord: The result of the Gordeyev integral over Doppler corrected radian frequency
        hs: The Debye length in m.
        Ns: The density of the species in m^-3
        omeg_s: An array of the Doppler corrected radian frequency
        """
        dFlag = self.dFlag
        K = self.K
        (Ns,Ts,Vs,qs,ms,nus) = dataline[:6]
        if len(dataline)==7:
            nuperp = dataline[-1]
        else:
            nuperp =nus

        C = sp.sqrt(v_Boltz*Ts/ms)
        omeg_s = self.omeg - K*Vs
        theta = omeg_s/(K*C*sp.sqrt(2.0))
        Om = qs*self.bMag/ms
        # determine what integral is used
        magbool = alpha*180.0/sp.pi < self.alphamax
        collbool = self.collfreqmin*K*C<nus

        if  not collbool and not magbool:
            #for case with no collisions or magnetic field just use analytic method
            gord = (sp.sqrt(sp.pi)*sp.exp(-theta**2)-1j*2.0*scipy.special.dawsn(theta))/(K*C*sp.sqrt(2))
            if dFlag:
                print('\t No collisions No magnetic field')
            return (gord,Ts,Ns,omeg_s)
        elif collbool and not magbool:
            if dFlag:
                print('\t With collisions No magnetic field')
            gordfunc = collacf
            exparams = (K,C,nus)
        elif not collbool and magbool:
            if dFlag:
                print('\t No collisions with magnetic field')
            gordfunc = magacf
            exparams = (K,C,alpha,Om)
        else:
            if dFlag:
                print('\t With collisions with magnetic field')
            gordfunc = magncollacf
            exparams = (K,C,alpha,Om,nus)


        maxf = sp.absolute(self.f).max()
        T_s = 1./(2.*maxf)

#        N_somm = 2**15
#        b1 = 10.0/(K*C*sp.sqrt(2.0))
        # changed ot sample grid better
        N_somm=2**10
        b1 = T_s*N_somm
#        b1 = interval/10.
#        N_somm=sp.minimum(2**10,sp.ceil(b1/T_s))
        (gord,flag_c,outrep) = sommerfelderfrep(gordfunc,N_somm,omeg_s,b1,Lmax=500,errF=1e-7,exparams=exparams)
        if dFlag:
            yna = ['No','Yes']
            print('\t Converged: {:s}, Number of iterations: {:d}'.format(yna[flag_c],outrep))

        return (gord,Ts,Ns,omeg_s)

def magacf(tau,K,C,alpha,Om):
    """ magacf(tau,K,C,alpha,Om)
        by John Swoboda
        This function will create a single particle acf for a particle species with magnetic
        field but no collisions.
        Inputs
        tau: The time vector for the acf.
        K: Bragg scatter vector magnetude.
        C: Thermal speed of the species.
        alpha: Magnetic aspect angle in radians.
        Om: The gyrofrequency of the particle.
        Output
        acf - The single particle acf.
        """
    Kpar = sp.sin(alpha)*K
    Kperp = sp.cos(alpha)*K
    return sp.exp(-sp.power(C*Kpar*tau,2.0)/2.0  -2.0*sp.power(Kperp*C*sp.sin(Om*tau/2.0)/Om,2.0))
def collacf(tau,K,C,nu):
    """ collacf(tau,K,C,nu)
        by John Swoboda
        This function will create a single particle acf for a particle species with no magnetic
        field but with collisions.
        Inputs
        tau: The time vector for the acf.
        K: Bragg scatter vector magnetude.
        C: Thermal speed of the species.
        nu: The collision frequency in collisions/sec
        Output
        acf - The single particle acf.
        """
    return sp.exp(-sp.power(K*C/nu,2.0)*(nu*tau-1+sp.exp(-nu*tau)))
def magncollacf(tau,K,C,alpha,Om,nu):
    """ magncollacf(tau,K,C,alpha,Om)
        by John Swoboda
        This function will create a single particle acf for a particle species with magnetic
        field and collisions.
        Inputs
        tau: The time vector for the acf.
        K: Bragg scatter vector magnetude.
        C: Thermal speed of the species.
        alpha: Magnetic aspect angle in radians.
        Om: The gyrofrequency of the particle.
        nu: The collision frequency in collisions/sec
        Output
        acf - The single particle acf.
        """
    Kpar = sp.sin(alpha)*K
    Kperp = sp.cos(alpha)*K
    gam = sp.arctan(nu/Om)

    deltl = sp.exp(-sp.power(Kpar*C/nu,2.0)*(nu*tau-1+sp.exp(-nu*tau)))
    deltp = sp.exp(-sp.power(C*Kperp,2.0)/(Om*Om+nu*nu)*(sp.cos(2*gam)+nu*tau-sp.exp(-nu*tau)*(sp.cos(Om*tau-2.0*gam))))
    return deltl*deltp

def Islistofstr(inlist):
    """ This function will determine if the input is a list of strings"""
    if not isinstance(inlist,list):
        return False
    for item in inlist:
        if not isinstance(item,string_types):
            return False
    return True

def allin(inp,reflist):
    """ This function will determine if all of the strings in the list input are in
    the list reflist."""
    for item in inp:
        if item not in reflist:
            return False
    return True

#%% Collision frequencies

def get_collisionfreqs(datablock,species,n_datablock=None, n_species=None):
    """ This function will calculate collision frequencies based off of the methods
    shown in Schunk and Nagy (2009) chapter 4.
        inputs
        datablock - A numpy array of size Nsp x2. The first element of the row is the
            density of species in m^-3 and the second element is the Tempreture in
            degrees K.
        species - A Nsp list of strings that label each species.
        col_calc - If this flag is true then collisions will be calculated. (Default= False)
        n_datablock - A numpy array of size Nnsp x2. The first element of the row is
            the density of the neutral species in m^-3 and the second element is the
            Tempreture in degrees K.If set to None then (Default=None)
        n_species - A Nnsp list of strings that label each neutral species.(Default=None)
        Outputs
        nuparr - A Nsp length numpy array that holds the parrallel collision frequencies in s^-1.
        nuperp - A Nsp length numpy array that holds the perpendictular collision frequencies in s^-1.
    """
    curpath = Path(__file__).parent

    nuperp = sp.zeros((datablock.shape[0]))
    nuparr = sp.zeros_like(nuperp)


    Ni = datablock[:-1,0]* 1e-6
    Ti = datablock[:-1,1]
    Ne = datablock[-1,0] *1e-6
    Te = datablock[-1,1]
    Bst = pd.read_csv(str(curpath/'ion2ion.csv'),index_col=0)

    Cin = pd.read_csv(str(curpath/'ion2neu.csv'),index_col=0)
    # electron electron and electron ion collisions
    #Schunk and Nagy eq 4.144 and eq 4.145
    nuee = 54.5/sp.sqrt(2) * Ne/sp.power(Te,1.5)
    nuei = 54.5*Ni/sp.power(Te,1.5)
    nuei = nuei.sum()

    nuperp[-1] = nuei +nuee
    nuparr[-1] = nuei
    # ionion collisions
    for si, s in enumerate(species[:-1]):
        for ti, t in enumerate(species[:-1]):
            nuparr[si]=nuparr[si]+Bst[s][t]*Ni[ti]/sp.power(Ti[ti],1.5)


    if (n_datablock is not None) and (n_species is not None):
        Nn = n_datablock[:,0]* 1e-6
        Tn = n_datablock[:,1]
        # ion neutral collisions
        for si, s in enumerate(species[:-1]):
            for ti, t in enumerate(n_species):

                curCin = Cin[s][t]

                if sp.isnan(curCin):
                    neuts = r_ion_neutral(s,t,Ni[si],Nn[ti],Ti[si],Tn[ti])
                else:
                    neuts = curCin* Ni[si]
                nuparr[si ]=nuparr[si]+neuts
        # electron neutral collision frequencies
        for ti, t in enumerate(n_species):
            nueneu =e_neutral(t,Nn[ti],Te)
            nuperp[-1] = nuperp[-1]+nueneu
            nuparr[-1]  = nuparr[-1]+ nueneu


    # according to milla and kudeki the collision for perpendicular and parrallel
    #directions only matter for the electrons
    nuperp[:-1]=nuparr[:-1]
    return(nuparr,nuperp)

def r_ion_neutral(s,t,Ni,Nn,Ti,Tn):
    """ This will calculate resonant ion - neutral reactions collision frequencies. See
    table 4.5 in Schunk and Nagy.
    Inputs
    s - Ion name string
    t - neutral name string
    Ni - Ion density cm^-3
    Nn - Neutral density cm^-3
    Ti - Ion tempreture K
    Tn - Neutral tempreture K
    Outputs
    nu_ineu - collision frequency s^-1
    """
    Tr = (Ti+Tn)*0.5
    sp1 = (s,t)
    # from Schunk and Nagy table 4.5
    nudict={('H+','H'):[2.65e-10,0.083],('He+','He'):[8.73e-11,0.093], ('N+','N'):[3.84e-11,0.063],
            ('O+','O'):[3.67e-11,0.064], ('N2+','N'):[5.14e-11,0.069], ('O2+','O2'):[2.59e-11,0.073],
            ('H+','O'):[6.61e-11,0.047],('O+','H'):[4.63e-12,0.],('CO+','CO'):[3.42e-11,0.085],
            ('CO2+','CO'):[2.85e-11,0.083]}
    A = nudict[sp1][0]
    B = nudict[sp1][1]
    if sp1==('O+','H'):

        nu_ineu = A*Nn*sp.power(Ti/16.+Tn,.5)
    elif sp1==('H+','O'):
        nu_ineu = A*Nn*sp.power(Ti,.5)*(1-B*sp.log10(Ti))**2
    else:
        nu_ineu = A*Nn*sp.power(Tr,.5)*(1-B*sp.log10(Tr))**2
    return nu_ineu

def e_neutral(t,Nn,Te):
    """ This will calculate electron - neutral reactions collision frequencies. See
    table 4.6 in Schunk and Nagy.

    Inputs
    t - neutral name string
    Nn - Neutral density cm^-3
    Te - electron tempreture K
    Outputs
    nu_ineu - collision frequency s^-1"""
    if t == 'N2':
        return 2.33e-11 * Nn * (1-1.21e-4 * Te) * Te
    elif t == 'O2':
        return 1.82e-10 * Nn * (1 - 3.6e-2 * Te**.5) * Te**.5
    elif t == 'O':
        return 8.9e-11 * Nn * (1 - 5.7e-4 * Te) * Te**.5
    elif t == 'He':
        return 4.6e-10 * Nn * Te**.5
    elif t == 'H':
        return 4.5e-9 * Nn * (1 - 1.35e-4 * Te) * Te**.5
    elif t == 'CO':
        return 2.34e-11 * Nn * (165 + Te)
    elif t == 'CO2':
        return 3.68e-8 *Nn*(1-4.1e-11 * sp.absolute(4500 - Te)**2.93)
    else:
        return 0.

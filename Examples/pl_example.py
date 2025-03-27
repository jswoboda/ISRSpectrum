
import numpy as np
import matplotlib.pyplot as plt
from ISRSpectrum import PLspecinit

#%% Test functions
def chapman_func(z, H_0, Z_0, N_0):
    """This function will return the Chapman function for a given altitude vector z.  All of the height values are assumed km.

    Parameters
    ----------
    z : ndarray
        An array of z values in km.
    H_0 : float
        A single float of the scale height in km.
    Z_0 : float
        The peak density location.
    N_0 : float
        The peak electron density.

    Returns
    -------
    Ne : ndarray
        Electron density as a function of z in m^{-3}
    """
    z1 = (z-Z_0)/H_0
    Ne = N_0*np.exp(0.5*(1-z1-np.exp(-z1)))
    return Ne

def temp_profile(z, T0=1000., z0=100.):
    """This function creates a tempreture profile using arc tan functions for test purposes.

    Parameters
    ----------
    z : ndarray
        An array of z values in km.
    T0 : ndarray
        The value of the lowest tempretures in K.
    z0 : ndarray
        The middle value of the atan functions along alitutude. In km.

    Returns
    -------
    Te : ndarray
        The electron density profile in K. 1700*(atan((z-z0)2*exp(1)/400-exp(1))+1)/2 +T0
    Ti : ndarray
        The ion density profile in K. 500*(atan((z-z0)2*exp(1)/400-exp(1))+1)/2 +T0
    """
    zall = (z-z0)*2.*np.exp(1)/400. -np.exp(1)
    atanshp = (np.tanh(zall)+1.)/2
    Te = 1700*atanshp+T0
    Ti = 500*atanshp+T0

    return (Te,Ti)

def example_profile():
    """ Creates profiles for ion temperature, electron density and electron temperature.
    Returns
    -------
    z : ndarray
        Height in km.
    Ne : ndarray
        Electron density as a function of z in m^{-3}
    Te : ndarray
        The electron density profile in K. 1700*(atan((z-z0)2*exp(1)/400-exp(1))+1)/2 +T0
    Ti : ndarray
        The ion density profile in K. 500*(atan((z-z0)2*exp(1)/400-exp(1))+1)/2 +T0
    """
    z = np.linspace(100,1000,100)
    N_0=1e12
    z_0=250.0
    H_0=50.0
    Ne_profile = chapman_func(z, H_0, z_0, N_0)
    (Te_prof, Ti_prof) = temp_profile(z)

    return (z,Ne_profile, Te_prof,Ti_prof)

def pl_demo():
    z, Ne, Te, Ti = example_profile()
    pls = PLspecinit(centerFrequency=440.2 * 1e6, bMag=0.4e-4, dFlag=False, fs=25e6, nchans=250, nfreq_pfb=1024)

    pl_cf = np.zeros((len(z),2))
    for inum,(irng,ine,ite,iti) in enumerate(zip(z,Ne,Te,Ti)):
        data_vec = np.array([[ine,iti],[ine,ite]])


        pld = pls.get_ul_spec(data_vec, 0., alphadeg=90.0, rcsflag=True, freqflag =True, Tpe=1.0)
        (f_lo,f_hi,sp_lo,sp_hi,rcs,freqs) = pld
        pl_cf[inum] = freqs



    fig = plt.figure()
    plt.plot(pl_cf*1e-6,z)
    plt.grid(True)
    plt.xlabel("Frequency in MHz")
    plt.ylabel("Altitude in km")
    plt.title('Plasma Frequency by Altitude.')
    fig.savefig("plfreqs.png",dpi=300)

if __name__== '__main__':
    pl_demo()

#!/usr/bin/env python

"""
This module will calcualte the Gordeyve integral for various cases, with collisions, with a magnetic field up to .1 degrees off of perp to B. This cut is controlled by an assert statement.
"""


import numpy as np
import scipy.special as sp_spec
import scipy.constants as spconst
import numpy as np
import scipy.fftpack as fftsy


class GordPlug:

    def calcgordeyev(
        self, dataline, alpha, K, omeg, bMag, collfreqmin=1e-2, alphamax=30, dFlag=False
    ):
        """
        Performs the Gordeyve integral calculation for main cases from the first Kudeki/Milla paper.

        Parameters
        ----------
        dataline : ndarray
            A numpy array of length that holds the plasma parameters needed to create the spectrum. Each row of the array will have the following set up.
                [Ns, Ts, Vs, qs, ms, nus]
                Ns - The density of the species in m^-3
                Ts - Temperature of the species in degrees K
                Vs - The Doppler velocity in m/s.
                qs - The charge of the species in elementary charges. (Value will be replaced for the electrons)
                ms - Mass of the species in AMU. (Value will be replaced for the electrons)
                nus - Collision frequency for species in s^-1.
        alphadeg : float
            The magnetic aspect angle in radians.
        K : float
            K value of the radar in rad/m
        omeg : ndarray
            Radian frequency vector the integral will be evaluated over.
        bMag : float
            Magnetic field in Teslas
        collfreqmin : float
            (Default 1e-2) The minimum collision frequency needed to incorporate it into Gordeyev integral calculations in units of K*sqrt(Kb*Ts/ms) for each ion species.
        alphamax : float
            (Default 30)  The maximum magnetic aspect angle in which the B-field will be taken into account.
        dFlag : float
            A debug flag, if set true will output debug text. Default is false.

        Returns
        -------
        gord : ndarray
            The result of the Gordeyev integral over Doppler corrected radian frequency
        hs : ndarray
            The Debye length in m.
        Ns : ndarray
            The density of the species in m^-3
        qs : float
            The charge of the species in elementary charges.
        omeg_s : ndarray
            An array of the Doppler corrected radian frequency
        """

        assert (
            alpha > 0.1
        ), "Angle off of perp to B must be greater than .1 degrees otherwise integral will not converge."
        assert (
            len(dataline) >= 6
        ), "The dataline input needs to be length of at least 6 elements."
        (Ns, Ts, Vs, qs, ms, nus) = dataline[:6]

        if qs < 0:
            qs = -spconst.e
            ms = spconst.m_e
        else:
            qs = qs * spconst.e
            ms = ms * spconst.m_p

        nur = np.pi * 2 * nus
        C = np.sqrt(spconst.k * Ts / ms)
        omeg_s = omeg - K * Vs
        theta = omeg_s / (K * C * np.sqrt(2.0))
        Om = qs * bMag / ms
        # determine what integral is used
        magbool = alpha * 180.0 / np.pi < alphamax
        collbool = collfreqmin * C * K < nur

        if not collbool and not magbool:
            # for case with no collisions or magnetic field just use analytic method
            num_g = np.sqrt(np.pi) * np.exp(-(theta**2)) - 1j * 2.0 * sp_spec.dawsn(
                theta
            )
            den_g = K * C * np.sqrt(2)
            gord = num_g / den_g

            if dFlag:
                print("\t No collisions No magnetic field,again")
            return (gord, Ts, Ns, qs, omeg_s)

        elif collbool and not magbool:
            if dFlag:
                print("\t With collisions No magnetic field")
            gordfunc = collacf
            exparams = (K, C, nur)
        elif not collbool and magbool:
            if dFlag:
                print("\t No collisions with magnetic field")
            gordfunc = magacf
            exparams = (K, C, alpha, Om)
        else:
            if dFlag:
                print("\t With collisions with magnetic field")
            gordfunc = magncollacf
            exparams = (K, C, alpha, Om, nur)

        maxf = np.abs(omeg / (2 * np.pi)).max()
        T_s = 1.0 / (2.0 * maxf)

        #        N_somm = 2**15
        #        b1 = 10.0/(K*C*np.sqrt(2.0))
        # changed ot sample grid better
        N_somm = 2**10
        b1 = T_s * N_somm
        #        b1 = interval/10.
        #        N_somm=np.minimum(2**10,np.ceil(b1/T_s))
        (gord, flag_c, outrep) = sommerfelderfrep(
            gordfunc, N_somm, omeg_s, b1, Lmax=500, errF=1e-7, exparams=exparams
        )
        if dFlag:
            yna = ["No", "Yes"]
            print(
                "\t Converged: {:s}, Number of iterations: {:d}".format(
                    yna[flag_c], outrep
                )
            )

        return (gord, Ts, Ns, qs, omeg_s)


def magacf(tau, K, C, alpha, Om):
    """Create a single particle acf for a species with magnetic field but no collisions.

    Parameters
    ----------
    tau : ndarray
        The time vector for the acf.
    K : float
        Bragg scatter vector magnetude.
    C : float
        Thermal speed of the species.
    alpha : float
        Magnetic aspect angle in radians.
    Om : float
        The gyrofrequency of the particle.

    Returns
    -------
    spacf : ndarray
        The single particle acf.
    """
    Kpar = np.sin(alpha) * K
    Kperp = np.cos(alpha) * K

    par_ex = -np.power(C * Kpar * tau, 2.0) / 2.0
    perp_ex = -2.0 * np.power((Kperp * C * np.sin(Om * tau / 2.0)) / Om, 2.0)
    spacf = np.exp(par_ex + perp_ex)
    return spacf


def collacf(tau, K, C, nu):
    """Create a single particle acf for a species with collisions, no magnetic field. See equation 48 in Kudeki and milla 2011.

    Parameters
    ----------
    tau : ndarray
        The time vector for the acf.
    K : float
        Bragg scatter vector magnetude.
    C : float
        Thermal speed of the species.
    nu : float
        The collision frequency in collisions*radians/sec

    Returns
    -------
    spacf : ndarray
        The single particle acf.
    """
    spacf = np.exp(-np.power(K * C / nu, 2.0) * (nu * tau - 1 + np.exp(-nu * tau)))
    return spacf


def magncollacf(tau, K, C, alpha, Om, nu):
    """Create a single particle acf for a species with magnetic field and collisions.

    Parameters
    ----------
    tau : ndarray
        The time vector for the acf.
    K : float
        Bragg scatter vector magnetude.
    C : float
        Thermal speed of the species.
    alpha : float
        Magnetic aspect angle in radians.
    Om : float
        The gyrofrequency of the particle.
    nu : float
        The collision frequency in collisions*radians/sec

    Returns
    -------
    acf : ndarray
        The single particle acf.
    """
    Kpar = np.sin(alpha) * K
    Kperp = np.cos(alpha) * K
    gam = np.arctan(nu / Om)

    deltl = np.exp(-np.power(Kpar * C / nu, 2.0) * (nu * tau - 1 + np.exp(-nu * tau)))
    deltp = np.exp(
        -np.power(C * Kperp, 2.0)
        / (Om * Om + nu * nu)
        * (
            np.cos(2 * gam)
            + nu * tau
            - np.exp(-nu * tau) * (np.cos(Om * tau - 2.0 * gam))
        )
    )
    return deltl * deltp


def chirpz(Xn, A, W, M):
    """Chirpz calculation for a single array.

    This function calculates the chirpz transfrom for the numpy array Xn given the complex constants A and W along with the length of the final array M.

    Parameters
    -----------
    Xn : ndarray
        The signal that the Chirp z transform will be calculated for.
    A : float
        A complex constant used to determine the direction of the integration in Z.
    W : float
        Another complex constant that will determine the direction of the integration in Z.
    M : int
        The length of the final chirpz transfrom.

    Returns
    -------
    yk : ndarray
        The M length chirp z tranfrom given Xn and the complex constants.
    """
    N = Xn.shape[0]
    # Make an L length output so the circular convolution does not wrap. Added 1 extra sample to make coding easier.
    L = np.power(2, np.ceil(np.log2(N + M - 1)))
    # chirpz numbering
    k = np.arange(M)
    # numbering for time axis
    n = np.arange(N)
    # Make complex arrays
    xn = np.zeros(L) + 1j * np.zeros(L)
    yn = np.zeros(L) + 1j * np.zeros(L)
    # complex constants raised to power
    An = A ** (-n)
    Wn = W ** (n**2 / 2.0)

    xn[:N] = Xn * An * Wn
    # Make the chirp kernal
    yn[:M] = W ** (-(k**2) / 2.0)
    nn = np.arange(L - N + 1, L)
    yn[L - N + 1 :] = W ** (-((L - nn) ** 2) / 2.0)
    # perform the circular convolution and multiply by chirp
    gk = fftsy.ifft(fftsy.fft(xn) * fftsy.fft(yn))[:M]
    yk = gk * W ** (k**2 / 2.0)

    return yk


def sommerfeldchirpz(
    func, N, M, dk, Lmax=1, errF=0.1, a=-1.0j, p=1.0, x0=None, exparams=()
):
    """Numerically integrate Sommerfeld like integral using chirpz.

    This function will numerically integrate a Sommerfeld like integral, int(exp(awk)f(k),t=0..inf) using at most Lmax number of N length chirpz transforms to make an M length array. If the normalized difference between the previous estimate of the output Xk is less then the parameter errF then the loop stops and a flag that represents convergence is set to true. A number of repeats are also output as well. This technique is based off the article Li et. al Adaptive evaluation of the Sommerfeld-type integral using the chirp z-transform, 1991. This function also uses a modified Simpsons rule for integration found in Milla PhD Thesis (2010).

    Parameters
    -----------
    func : func
        A function that is used to create f(k).
    N : int
        The length of the chirp z transform used.
    M : int
        The length of the output array.
    dk : float
        The sample period of k.
    Lmax : int
        The maximum number of repeats of the integration before the loop finishes, default 1.
    errF : float
        The threshold of the normalized difference between the new iteration and the old to stop to stop the iteration, default .1.
    a : float
        A complex number that determines the trejectory of the integration on the z plane, default -1.0*1j.
    p : float
        A real number that helps to define the spacing on the omega plane.
    x0 : float
        The starting point on the omega plane, default p*np.pi/dk,
    exparams :tuple
        Any extra params other then k to create f(k), default ().

    Returns
    -------
    Xk : ndarray
        The integrated data that is of length M.
    flag_c : bool
        A convergence flag.
    outrep : int
        The number of repetitions until convergence.
    """

    k = np.arange(N) * dk
    if x0 is None:
        x0 = p * np.pi / dk

    # set up simpson rule
    wk = np.ones(N)
    wk[np.mod(k, 2) == 0] = 2.0 / 3.0
    wk[np.mod(k, 2) == 1] = 4.0 / 3.0
    wk[0] = 1.5 / 3.0
    wk[-1] = 1.5 / 3.0
    # Make A_0 and W_0
    A_0 = np.exp(a * dk * x0)
    M2 = M - np.mod(M, 2)
    W_0 = np.exp(a * 2.0 * p * np.pi / M2)
    # Make frequency array for phase offset
    freqm = np.arange(-np.ceil((M - 1) / 2.0), np.floor((M - 1) / 2.0) + 1)
    Xk = np.zeros(M) * 1j
    flag_c = False

    for irep in range(Lmax):

        fk = func(k + N * dk * irep, *exparams)
        Xkold = Xk
        Xk = chirpz(fk * wk, A_0, W_0, M) * np.power(W_0, N * dk * irep * freqm) + Xk

        Xkdiff = np.sqrt(np.sum(np.power(np.abs(Xk - Xkold), 2.0)))
        Xkpow = np.sqrt(np.sum(np.power(np.abs(Xk), 2.0)))
        outrep = irep + 1
        # check for convergence
        if Xkdiff / Xkpow < errF:
            flag_c = True

            break

    return (Xk, flag_c, outrep)


def sommerfelderfrep(func, N, omega, b1, Lmax=1, errF=0.1, exparams=()):
    """Numerically integrate Sommerfeld like integral using erf transform function loop.

    This function will numerically integrate a Sommerfeld like integral, int(exp(-jwk)f(k),k=0..inf) using the ERF transform and 2N+1 samples and at most Lmax loops. If the normalized difference between the previous estimate of the output Xk is less then the parameter errF then the loop stops and a flag that represents convergence is set to true. A number of loops is also output as well.

    This function uses sommerfelderf to do the integration

    Parameters
    -----------
    func : func
        A function that is used to create f(k).
    N : int
        The integration uses 2N+1 samples to do the integration.
    omega : ndarray
        The Frequency array in radians per second.
    b1 : int
        The inital bounds of the first try and then step size for each subsiquent integral.
    Lmax : int
        The maximum number of repeats of the integration before the loop finishes, default 1.
    errF : float
        The threshold of the normalized difference between the new iteration and the old to stop to stop the iteration, default .1.
    exparams :tuple
        Any extra params other then k to create f(k), default ().

    Returns
    -------
    Xk : ndarray
        The integrated data that is of length M.
    flag_c : bool
        A convergence flag.
    outrep : int
        The number of repetitions until convergence.
    """
    Xk = np.zeros_like(omega) * (1 + 1j)
    flag_c = False
    for irep in range(Lmax):

        Xktemp = sommerfelderf(func, N, omega, b1 * irep, b1 * (irep + 1), exparams)
        Xkdiff = Xktemp.real**2 + Xktemp.imag**2
        Xk = Xk + Xktemp
        Xkpow = Xk.real**2 + Xk.imag**2
        #        Xkdiff = np.sqrt(np.sum(np.power(np.abs(Xktemp),2.0)))
        #        Xkpow = np.sqrt(np.sum(np.power(np.abs(Xk+Xktemp),2.0)))
        #        Xk = Xk+Xktemp
        outrep = irep + 1
        # check for convergence
        if np.sum(Xkdiff / Xkpow) < errF:
            flag_c = True
            outrep = irep
            break
    return (Xk, flag_c, outrep)


def sommerfelderf(func, N, omega, a, b, exparams=()):
    """Integrate somerfeld integral using ERF transform for single portion.

    This function will numerically integrate a Sommerfeld like integral, int(exp(-jwk)f(k),k=a..b) using the ERF transform and 2N+1 samples. This technique is from the paper B. L. Ooi 2007.

    Parameters
    -----------
    func : func
        A function that is used to create f(k).
    N : int
        The integration uses 2N+1 samples to do the integration.
    omega : ndarray
        The Frequency array in radians per second.
    a : float
        Lower bound of the integral.
    b : float
        Upper bound of teh integral.
    exparams :tuple
        Any extra params other then k to create f(k), default ().

    Returns
    -------
    Xk : ndarray
        The integrated data that is of length of omega.
    """

    nvec = np.arange(-N, N + 1)

    h = np.log(1.05 * np.sqrt(2 * N)) / N
    kn = 0.5 * (b + a) + 0.5 * (b - a) * sp_spec.erf(np.sinh(nvec * h))

    An = np.cosh(nvec * h) * np.exp(-np.power(np.sinh(nvec * h), 2))

    fk = func(kn, *exparams)
    kmat = np.tile(kn[:, np.newaxis], (1, len(omega)))
    omegamat = np.tile(omega[np.newaxis, :], (len(kn), 1))
    Xk3 = (An * fk).dot(np.exp(-1j * kmat * omegamat))

    return Xk3 * h * (b - a) / np.sqrt(np.pi)

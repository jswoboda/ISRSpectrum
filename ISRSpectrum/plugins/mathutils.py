#!python

"""
A number of mathematical tools needed to calculate the spectra.
"""
import numpy as np
import scipy.fftpack as fftsy
import scipy.special


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
    kn = 0.5 * (b + a) + 0.5 * (b - a) * scipy.special.erf(np.sinh(nvec * h))

    An = np.cosh(nvec * h) * np.exp(-np.power(np.sinh(nvec * h), 2))

    fk = func(kn, *exparams)
    kmat = np.tile(kn[:, np.newaxis], (1, len(omega)))
    omegamat = np.tile(omega[np.newaxis, :], (len(kn), 1))
    Xk3 = (An * fk).dot(np.exp(-1j * kmat * omegamat))

    return Xk3 * h * (b - a) / np.sqrt(np.pi)

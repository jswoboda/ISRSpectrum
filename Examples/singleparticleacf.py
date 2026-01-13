#!/usr/bin/env python
"""
This will create images of the signal particle ACFs using the default gordyev plugin
"""

from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
import scipy.constants as spconst
import scipy.fftpack as fftsy
import scipy.special as sp_spec
from matplotlib import rc
from matplotlib.markers import MarkerStyle

#
from ISRSpectrum.plugins.gord_default import (
    collacf,
    magacf,
    magncollacf,
    simpacf,
    sommerfelderfrep,
)


def singleparticalacfs():
    # %% Sim set up
    centerFrequency = 440.2 * 1e6
    nspec = 1024
    sampfreq = 50e3
    bMag = 0.4e-4
    Ts = 1e3
    Partdict = {
        0: ("Electron", spconst.m_e, spconst.e),
        1: ("H+ Ion", spconst.m_p, spconst.e),
        2: ("O+ Ion", 16 * spconst.m_p, spconst.e),
    }
    curpart = 1
    particle = Partdict[curpart]
    pname = particle[0]

    ms = particle[1]
    q_ch = particle[2]
    K = 2.0 * np.pi * 2 * centerFrequency / spconst.c

    f = np.arange(-np.ceil((nspec - 1.0) / 2.0), np.floor((nspec - 1.0) / 2.0 + 1)) * (
        sampfreq / (2 * np.ceil((nspec - 1.0) / 2.0))
    )
    C = np.sqrt(spconst.k * Ts / ms)
    Om = q_ch * bMag / ms / (K * C * np.sqrt(2))

    omeg_s = 2.0 * np.pi * f
    theta = omeg_s / (K * C * np.sqrt(2.0))

    dtau = 2e-3
    N = 2**12
    N_somm = 2**7
    b1 = N_somm * dtau
    tau = np.arange(N) * dtau
    d2r = np.pi / 180.0

    # %% No collisions or magnetic field
    gordnn = np.exp(-np.power(tau, 2.0) / 4.0)

    fig1 = plt.figure()
    plt.plot(tau, gordnn, linewidth=3)
    plt.title(r"Single Particle ACF for " + pname)
    plt.grid(True)
    plt.savefig("ACF" + pname.replace(" ", "") + ".png")
    plt.close(fig1)
    num_g = np.sqrt(np.pi) * np.exp(-(theta**2)) - 1j * 2.0 * sp_spec.dawsn(theta)
    den_g = K * C * np.sqrt(2)
    gord = num_g / den_g

    (Gn, flag_c, outrep) = sommerfelderfrep(
        simpacf, N_somm, theta, b1, Lmax=500, errF=1e-7
    )
    gordn = Gn / den_g

    f_khz = f * 1e-3
    fig2, ax = plt.subplots(1, 2)
    ax[0].plot(f_khz, gord.real, label="analytic Real")
    ax[1].plot(f_khz, gord.imag, label="analytic imag")
    ax[0].plot(f_khz, gordn.real, label="numerical Real", linestyle="dashed")
    ax[1].plot(f_khz, gordn.imag, label="numerical imag", linestyle="dashed")
    ax[0].set_xlabel("Frequency kHz")
    ax[1].set_xlabel("Frequency kHz")
    ax[0].set_ylabel("Ampltiude")

    ax[0].set_title("Gordyev Real Part")
    ax[1].set_title("Gordyev Imag Part")
    ax[0].grid(True)
    ax[1].grid(True)
    fig2.tight_layout()
    plt.savefig("gord" + pname.replace(" ", "") + ".png")
    plt.close(fig2)

    # %% With collisions
    orig_vec = np.logspace(-5.0, -1, 5)
    zvec = np.zeros_like(orig_vec)
    nucoll_vec = np.concat([orig_vec, zvec, orig_vec], axis=0)
    nuneu_vec = np.concat([zvec, orig_vec, orig_vec])

    sp_acf = np.zeros((len(nuneu_vec), len(tau)))
    alpha = np.pi / 2
    for inum, (inuc, inuneu) in enumerate(zip(nucoll_vec, nuneu_vec)):
        sp_acf[inum, :] = magncollacf(tau, K, C, alpha, Om, inuc / 2, inuc / 2, inuneu)
    # np.exp(-np.power(K*C/numat,2.0)*(numat*taumat-1+np.exp(-numat*taumat)))

    plt.figure()
    plt.plot(
        tau, sp_acf, linestyle="--", color="b", linewidth=4, label=r"No Collisions"
    )

    for inun, (inuc, inuneu) in enumerate(zip(nucoll_vec, nuneu_vec)):
        curlabel = r"$\nu_c$ = {:.2f} KC, $\nu_n$ = {:.2f} KC".format(inuc, inuneu)
        plt.plot(
            tau,
            sp_acf[inun].real,
            linewidth=3,
            label=curlabel,
        )

    plt.grid(True)
    plt.title(r"Single Particle ACF W/ Collisions for " + pname)
    plt.legend()
    plt.savefig("ACFwcolls" + pname.replace(" ", "") + ".png")

    # %% With magnetic field
    alpha = np.linspace(19, 1, 10)
    taumat = np.tile(tau[np.newaxis, :], (len(alpha), 1))
    almat = np.tile(alpha[:, np.newaxis], (1, len(tau)))
    Kpar = np.sin(d2r * almat) * K
    Kperp = np.cos(d2r * almat) * K

    #    gordmag = np.exp(-np.power(C*Kpar*taumat,2.0)/2.0-2.0*np.power(Kperp*C*np.sin(Om*taumat/2.0)/Om,2.0))
    gordmag = magacf(taumat, K, C, d2r * almat, Om)
    plt.figure()
    plt.plot(tau, gordnn, linestyle="--", color="b", linewidth=4, label="No B-field")

    for ialn, ial in enumerate(alpha):
        plt.plot(
            tau,
            gordmag[ialn].real,
            linewidth=3,
            label=r"$\alpha = {:.0f}^\circ$".format(ial),
        )

    plt.grid(True)
    plt.title("Single Particle ACF W/ Mag for " + pname)
    plt.legend()
    plt.savefig("ACFwmag" + pname.replace(" ", "") + ".png")

    # # %% Error surface with both
    # almat3d = np.tile(alpha[:, np.newaxis, np.newaxis], (1, len(nuvec), len(tau)))
    # numat3d = np.tile(nuvec[np.newaxis, :, np.newaxis], (len(alpha), 1, len(tau)))
    # taumat3d = np.tile(tau[np.newaxis, np.newaxis, :], (len(alpha), len(nuvec), 1))
    # Kpar3d = np.sin(d2r * almat3d) * K
    # Kperp3d = np.cos(d2r * almat3d) * K
    # gam = np.arctan(numat3d / Om)

    #    deltl = np.exp(-np.power(Kpar3d*C/numat3d,2.0)*(numat3d*taumat3d-1+np.exp(-numat3d*taumat3d)))
    #    deltp = np.exp(-np.power(C*Kperp3d,2.0)/(Om*Om+numat3d*numat3d)*(np.cos(2*gam)+numat3d*taumat3d-np.exp(-numat3d*taumat3d)*(np.cos(Om*taumat3d-2.0*gam))))

    # #    gordall = deltl*deltp
    # gordall = magncollacf(taumat3d, K, C, d2r * almat3d, Om, numat3d)
    # gordnnmat = np.tile(gordnn[np.newaxis, np.newaxis, :], (len(alpha), len(nuvec), 1))

    # gorddiff = np.abs(gordall - gordnnmat) ** 2
    # err = np.sqrt(gorddiff.mean(2)) / np.sqrt(np.power(gordnn, 2.0).sum())
    # extent = [
    #     np.log10(nuvec[0] / (K * C)),
    #     np.log10(nuvec[-1] / (K * C)),
    #     alpha[0],
    #     alpha[-1],
    # ]

    # plt.figure()
    # myim = plt.imshow(err * 100, extent=extent, origin="lower", aspect="auto")
    # myim.set_clim(0.0, 5.0)
    # plt.xlabel(r"$\log_{10}(\nu /KC)$")
    # plt.ylabel(r"$^\circ\alpha$")
    # cbar = plt.colorbar()
    # cbar.set_label("% Error", rotation=270)
    # cbar.ax.get_yaxis().labelpad = 15
    # plt.title("Error between ideal ACF and with Collisions and B-field for " + pname)
    # plt.savefig("ACFerr" + pname.replace(" ", "") + ".png")


if __name__ == "__main__":
    singleparticalacfs()

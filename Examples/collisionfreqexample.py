#!/usr/bin/env python
"""
Created on Thu Nov 12 10:51:26 2015

@author: John Swoboda
"""

from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
import scipy.constants as spconst
import seaborn as sns
from pyparsing.common import OneOrMore
from seaborn._core.typing import DataSource

#
from ISRSpectrum import INFODICT, Specinit, get_collisionfreqs
from ISRSpectrum.plugins.gord_default import collacf, magacf, magncollacf


def main():
    sns.set_style("whitegrid")
    sns.set_context("notebook")
    n_species = ["H", "He", "N", "N2", "O", "O2"]
    Nn = (
        np.array(
            [
                22340974.0,
                1.170203e08,
                269365.12,
                8.7995607e12,
                5.6790339e11,
                1.9389578e12,
            ]
        )
        * 1e6
    )
    Tn = 179.29138
    n_datablock = np.zeros((len(Nn), 2))
    n_datablock[:, 0] = Nn
    n_datablock[:, 1] = Tn

    species = ["NO+", "O2+", "e-"]
    Ni = np.array([1998.2988800000001, 57.700355999999999]) * 1e6
    Ti = 175.1
    Ne = 2055.9992320000001 * 1e6
    Te = 175.1
    datablock = np.zeros((len(Ni) + 1, 2))
    datablock[:-1, 0] = Ni
    datablock[:-1, 1] = Ti
    datablock[-1, 0] = Ne
    datablock[-1, 1] = Te

    nuparr, nuperp, nuneu = get_collisionfreqs(
        datablock, species, n_datablock=n_datablock, n_species=n_species
    )
    cent_freq = 449e6
    bMag = 4e-4
    K = 2.0 * np.pi * 2 * cent_freq / spconst.c

    print(f"Parallel collisions: {nuparr}")
    print(f"Perpendicular collisions: {nuperp}")
    print(f"Neutral collisions: {nuneu}")
    m_list = np.array([INFODICT[isp][-1] for isp in species]).astype(float)
    m_ar = np.empty_like(m_list, dtype=float)
    m_ar[:-1] = m_list[:-1] * spconst.m_p
    m_ar[-1] = spconst.m_e
    C_vec = np.sqrt(spconst.k * datablock[:, -1] / m_ar)
    onorm = K * C_vec * np.sqrt(2)

    Om_vec = spconst.e * bMag / m_ar / onorm
    nuparr_n = nuparr / onorm
    nuperp_n = nuperp / onorm
    nuneu_n = nuneu / onorm
    n_pnts = int(2**12)
    dtau = 2e-3
    tau = np.arange(n_pnts) * dtau
    gordnun = np.empty((3, n_pnts), dtype=float)
    for inum in range(len(species)):
        gordnun[inum, :] = magncollacf(
            tau,
            np.pi / 2,
            Om_vec[inum],
            nuparr_n[inum],
            nuperp_n[inum],
            nuneu_n[inum],
        )
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    for inum, igord in enumerate(gordnun):
        ax[0].plot(tau, igord, label=species[inum])
        ax[1].plot(tau, igord, label=species[inum])
    ax[0].plot(tau, np.exp(-(tau**2)), label="simp")
    ax[1].plot(tau, np.exp(-(tau**2)), label="simp")
    ax[1].set_yscale("log")
    ISS1 = Specinit(
        centerFrequency=cent_freq, bMag=0.4e-4, nspec=512, sampfreq=25e3, dFlag=True
    )

    (omeg, spec_1) = ISS1.getspecsep(datablock, species, vel=0.0, rcsflag=False)
    (omeg, spec_2) = ISS1.getspecsep(
        datablock, species, vel=0.0, rcsflag=False, col_calc=True
    )
    (omeg, spec_3) = ISS1.getspecsep(
        datablock,
        species,
        vel=0.0,
        rcsflag=False,
        col_calc=True,
        n_datablock=n_datablock,
        n_species=n_species,
    )

    (figmplf, axmat) = plt.subplots(1, 1, figsize=(20, 15), facecolor="w")
    curax = axmat
    lines = [None] * 3
    labels = [
        "Spectrum",
        "Spectrum with Coulomb Collisions",
        "Spectrum with Columb and Neutral Collisions",
    ]
    lines[0] = curax.plot(
        omeg * 1e-3, spec_1 / np.max(spec_1), label="Output", linewidth=5
    )[0]
    lines[1] = curax.plot(
        omeg * 1e-3, spec_2 / np.max(spec_2), label="Output", linewidth=5
    )[0]
    lines[2] = curax.plot(
        omeg * 1e-3, spec_3 / np.max(spec_3), label="Output", linewidth=5
    )[0]

    figmplf.suptitle("Spectrums with and without collisions")
    plt.figlegend(lines, labels, loc="lower center", ncol=5, labelspacing=0.0)

    plt.savefig("Collisions.png")


if __name__ == "__main__":
    main()

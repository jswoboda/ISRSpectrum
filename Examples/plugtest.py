#!python

import numpy as np
import scipy.constants as spconst
from ISRSpectrum import Specinit
import matplotlib.pyplot as plt


def main_plugtest():
    """This is an example of using the pluggin capability."""
    spfreq = 50e3
    nspec = 512
    ISS2 = Specinit(
        centerFrequency=440.2 * 1e6,
        bMag=0.4e-4,
        nspec=nspec,
        sampfreq=spfreq,
        dFlag=True,
    )

    ti = 1e3
    te = 1e3
    Ne = 1e11
    Ni = 1e11
    vi = 0
    mi = 16
    species = ["O+", "e-"]

    datablock = np.array([[Ni, ti], [Ne, te]])
    Nui = 0  # Collision frequency of ions
    Nue = 0  # Collision frequency of electrons
    i_list_s = [Ni, ti, vi, 1, mi, Nui]
    e_list_s = [Ne, te, vi, -1, 1, Nue]
    datablock_s = np.array([i_list_s, e_list_s])
    (omega, specorig, rcs) = ISS2.getspec(datablock_s, des_plug="default", rcsflag=True)

    i_list_simp = [Ni, ti, vi, 1, mi * spconst.m_p]
    e_list_simp = [Ne, te, vi, -1, spconst.m_e]
    datablocksimp = np.array([i_list_simp, e_list_simp])

    (omega, specorig_simp, rcs) = ISS2.getspec(
        datablocksimp, des_plug="simple", rcsflag=True
    )

    Nui = 10e3
    i_list_c = [Ni, ti, vi, 1, mi, Nui]
    e_list_c = [Ne, te, vi, -1, 1, Nue]
    datablockdefault = np.array([i_list_c, e_list_c])
    (omega, specorig_default, rcs) = ISS2.getspec(
        datablockdefault, des_plug="default", rcsflag=True
    )

    plt.plot(
        1e-3 * omega,
        specorig_default / specorig_default.max(),
        linewidth=3,
        label="Default plugin",
    )
    plt.plot(
        1e-3 * omega,
        specorig_simp / specorig_simp.max(),
        linewidth=3,
        label="Simple plugin",
    )

    plt.ylabel("Normalized Amplitude")
    plt.xlabel("Frequency in kHz")
    plt.title("Spectra from different plugins")
    plt.grid(True)
    plt.legend()
    plt.savefig("plugintest.png")


if __name__ == "__main__":

    main_plugtest()

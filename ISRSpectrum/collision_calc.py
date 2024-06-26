#!python
"""
collision_calc.py

This module is used to calculate collision frequencies.
"""
import numpy as np
from pathlib import Path
import pandas as pd

# %% Collision frequencies


# AMU for important molicules.
INFODICT = {
    "O+": np.array([1, 16]),
    "NO+": np.array([1, 30]),
    "N2+": np.array([1, 28]),
    "O2+": np.array([1, 32]),
    "N+": np.array([1, 14]),
    "H+": np.array([1, 1]),
    "He+": np.array([1, 4]),
    "e-": np.array([-1, 1]),
}


def ioncheck(ionspecies):
    """Used for assert statments to check if species is valid.

    Parameters
    ----------
    ionspecies : list
        Names of ionspecies used. Includes 'O+', 'NO+', 'N2+', 'O2+', 'N+', 'H+', 'He+'.

    Returns
    -------
    allgood : boolean
        True if all list members are valid species.
    """

    maslist = INFODICT.keys()

    allgood = True
    for ion in ionspecies:
        if not ion in maslist:
            print("{} not a recognized ionspecies".format(ion))
            allgood = False
    return allgood


def getionmass(ion):
    """Gets ion mass of a species.

    Parameters
    ----------
    ion : str
        Name of ionspecies used. Includes 'O+', 'NO+', 'N2+', 'O2+', 'N+', 'H+', 'He+'.

    Returns
    -------
    int
        Atomic mass of species.
    """

    assert ioncheck([ion])
    return INFODICT[ion][-1]


def get_collisionfreqs(
    datablock, species, Bstfile=None, Cinfile=None, n_datablock=None, n_species=None
):
    """Calculate collision frequencies

    Uses the methods shown in Schunk and Nagy (2009) chapter 4.

    Parameters
    ----------
    datablock : ndarray
        A numpy array of size Nsp x2. The first element of the row is the density of species in m^-3 and the second element is the temperature in degrees K.
    species : list
        A Nsp list of strings that label each species.
    Bstfile : str
        CSV file holding Ion ion collision constants
    Cinfile : ndarray
        CSV file holding Ion neutral collision
    n_datablock : ndarray
        A numpy array of size Nnsp x2. The first element of the row is the density of the neutral species in m^-3 and the second element is the temperature in degrees K. If set to None then (Default=None)
    n_species : ndarray
        A Nnsp list of strings that label each neutral species.(Default=None)

    Returns
    -------
    nuparr : ndarray
        A Nsp length numpy array that holds the parrallel collision frequencies in s^-1.
    nuperp : ndarray
        A Nsp length numpy array that holds the perpendictular collision frequencies in s^-1.
    """

    nuperp = np.zeros((datablock.shape[0]))
    nuparr = np.zeros_like(nuperp)

    Ni = datablock[:-1, 0] * 1e-6
    Ti = datablock[:-1, 1]
    Ne = datablock[-1, 0] * 1e-6
    Te = datablock[-1, 1]
    curpath = Path(__file__).parent
    if Bstfile is None:
        Bstfile = str(curpath.joinpath("ion2ion.csv"))

    Bst = pd.read_csv(Bstfile, index_col=0)

    if Cinfile is None:
        Cinfile = str(curpath.joinpath("ion2neu.csv"))
    Cin = pd.read_csv(Cinfile, index_col=0)
    # electron electron and electron ion collisions
    # Schunk and Nagy eq 4.144 and eq 4.145
    nuee = 54.5 / np.sqrt(2) * Ne / np.power(Te, 1.5)
    nuei = 54.5 * Ni / np.power(Te, 1.5)
    nuei = nuei.sum()

    nuperp[-1] = nuei + nuee
    nuparr[-1] = nuei
    # ionion collisions
    for si, s in enumerate(species[:-1]):
        for ti, t in enumerate(species[:-1]):

            nuparr[si] = nuparr[si] + Bst[s][t] * Ni[ti] / np.power(Ti[ti], 1.5)

    if (n_datablock is not None) and (n_species is not None):
        Nn = n_datablock[:, 0] * 1e-6
        Tn = n_datablock[:, 1]
        # ion neutral collisions
        for si, s in enumerate(species[:-1]):
            for ti, t in enumerate(n_species):

                curCin = Cin[s][t]

                if np.isnan(curCin):
                    neuts = r_ion_neutral(s, t, Ni[si], Nn[ti], Ti[si], Tn[ti])
                else:
                    neuts = curCin * Ni[si]
                nuparr[si] = nuparr[si] + neuts
        # electron neutral collision frequencies
        for ti, t in enumerate(n_species):
            nueneu = e_neutral(t, Nn[ti], Te)
            nuperp[-1] = nuperp[-1] + nueneu
            nuparr[-1] = nuparr[-1] + nueneu

    # according to milla and kudeki the collision for perpendicular and parrallel
    # directions only matter for the electrons
    nuperp[:-1] = nuparr[:-1]
    return (nuparr, nuperp)


def r_ion_neutral(s, t, Ni, Nn, Ti, Tn):
    """This will calculate resonant ion - neutral reactions collision frequencies. See table 4.5 in Schunk and Nagy.

    Parameters
    ----------
    s : str
        Ion species name
    t : str
        Neutral species name
    Ni : float
        Ion density cm^-3
    Nn : float
        Neutral density cm^-3
    Ti : float
        Ion temperature K
    Tn float
        Neutral temperature K

    Returns
    -------
    nu_ineu : float
        collision frequency s^-1
    """
    Tr = (Ti + Tn) * 0.5
    sp1 = (s, t)
    # from Schunk and Nagy table 4.5
    nudict = {
        ("H+", "H"): [2.65e-10, 0.083],
        ("He+", "He"): [8.73e-11, 0.093],
        ("N+", "N"): [3.84e-11, 0.063],
        ("O+", "O"): [3.67e-11, 0.064],
        ("N2+", "N"): [5.14e-11, 0.069],
        ("O2+", "O2"): [2.59e-11, 0.073],
        ("H+", "O"): [6.61e-11, 0.047],
        ("O+", "H"): [4.63e-12, 0.0],
        ("CO+", "CO"): [3.42e-11, 0.085],
        ("CO2+", "CO"): [2.85e-11, 0.083],
    }
    A = nudict[sp1][0]
    B = nudict[sp1][1]
    if sp1 == ("O+", "H"):

        nu_ineu = A * Nn * np.power(Ti / 16.0 + Tn, 0.5)
    elif sp1 == ("H+", "O"):
        nu_ineu = A * Nn * np.power(Ti, 0.5) * (1 - B * np.log10(Ti)) ** 2
    else:
        nu_ineu = A * Nn * np.power(Tr, 0.5) * (1 - B * np.log10(Tr)) ** 2
    return nu_ineu


def e_neutral(t, Nn, Te):
    """This will calculate electron - neutral reactions collision frequencies. See
    table 4.6 in Schunk and Nagy.

    Parameters
    ----------
    t : str
        Neutral species name
    Nn : float
        Neutral density cm^-3
    Te : float
        electron temperature K

    Returns
    -------
    nu_ineu : float
        Collision frequency s^-1
    """
    if t == "N2":
        return 2.33e-11 * Nn * (1 - 1.21e-4 * Te) * Te
    elif t == "O2":
        return 1.82e-10 * Nn * (1 - 3.6e-2 * Te**0.5) * Te**0.5
    elif t == "O":
        return 8.9e-11 * Nn * (1 - 5.7e-4 * Te) * Te**0.5
    elif t == "He":
        return 4.6e-10 * Nn * Te**0.5
    elif t == "H":
        return 4.5e-9 * Nn * (1 - 1.35e-4 * Te) * Te**0.5
    elif t == "CO":
        return 2.34e-11 * Nn * (165 + Te)
    elif t == "CO2":
        return 3.68e-8 * Nn * (1 - 4.1e-11 * np.abs(4500.0 - Te) ** 2.93)
    else:
        return 0.0

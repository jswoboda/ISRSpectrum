Introduction
============

This is a Python module to calculate an incoherent scatter spectrum based off of Kudeki and Milla's 2011 IEEE Geophysics paper.

    Kudeki, E.; Milla, M.A., "Incoherent Scatter Spectral Theories—Part I: A General Framework and Results for Small Magnetic Aspect Angles," Geoscience and Remote Sensing, IEEE Transactions on , vol.49, no.1, pp.315,328, Jan. 2011

 Like the model covered in the paper the software can calculate a spectra given, any magnetic aspect angle not perpendicular to B, any number of ion species, and the collision frequencies associated with those ion species. As the magnetic aspect angles get closer to perpendicular to B, usually &lt; 1 degree perp to B, more and more calculations are needed for the Gordeyev to converge.


Install
**********

This module is made from pure Python, it can be installed via pip after downloading it from GitHub

Usage
**********

There are a number of examples in the Examples directory that show how to use the software. The basic way to use the code is to import the Specinit object. Then the object can then call one of three different methods depending on how your parameters are defined. For example the getspecsimple method takes physical electron density, temperature, ion temperature, velocity and the species along with their concenterations. The getspecsep and getspec methods uses a numpy array to hold physical information about the spectra. See the documentation for more details.


Example
**********

This example shows how to use the Specinit class and then create a spectrum using the getspecsimple,  getspecsep and getspec. The two functions getspecsimple and getspecsep take the inputs different but calculate collision frequencies based off of the physical inputs, i.e. density, temperture, ion species. The method getspec takes as an input, in the datablock array, the collision frequency so the user can have more control over defining the frequency

.. code-block:: python

    import numpy as np
    import scipy.constants as spconsts
    from ISRSpectrum import Specinit
    

    spfreq = 50e3
    nspec = 512
    ISS2 = Specinit(centerFrequency = 440.2*1e6, bMag = 0.4e-4, nspec=nspec, sampfreq=spfreq,dFlag=True)

    ti = 1e3
    te = 1e3
    Ne = 1e11
    Ni = 1e11
    vi = 0
    mi = 16
    species = ['O+','e-']

    omeg,spec = ISS2.getspecsimple(Ne,te,ti,['O+'],[1])
    datablock = np.array([[Ni,ti],[Ne,te]])
    (omega,specorig,rcs) = ISS2.getspecsep(datablock,species,vi,alphadeg=90,rcsflag = True)
    Nui = 0 # Collision frequency of ions
    Nue = 0 # Collision frequency of electrons
    datablock90 = np.array([[Ni,ti,vi,1,mi,Nui],[Ne,te,vi,-1,1,Nue]])
    (omega,specorig,rcs) = ISS2.getspec(datablock90, rcsflag = True)

.. figure:: imgs/introexamplespec.png
   :class: with-border

   Example output spectrum
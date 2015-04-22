# ISRSpectrum
![alt text](https://raw.github.com/jswoboda/ISRSpectrum/master/logofig.png "ISR Spectrum")
# Overview
This is a Python module to calculate an ISR spectrum  based off of Kudeki and Milla's 2011 IEEE Geophysics paper. 

	Kudeki, E.; Milla, M.A., "Incoherent Scatter Spectral Theories—Part I: A General Framework and Results for Small Magnetic Aspect Angles," Geoscience and Remote Sensing, IEEE Transactions on , vol.49, no.1, pp.315,328, Jan. 2011

The code has been written to be able to produce a spectrum with an arbitrary number of ion species with an arbitrary collsion frequency and magnetic aspect angle. The only issue is that I will now promise that the code will run quickly at magnetic aspect angles <1 degree perp to B. 
 
#Installation

To install first clone repository:

	$ git clone https://github.com/jswoboda/ISRSpectrum.git

Then move to the const directory. This is a directory that was turned into a submodule because it was being used by multiple repositories.

	$ cd ISRSpectrum/ISRSpectrum/const
	$ git pull origin
	
Run the Python setup script. I would suggest running it in development mode.
	



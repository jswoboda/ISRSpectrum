Introduction
============

This is a Python module to calculate an incoherent scatter spectrum based off of Kudeki and Milla's 2011 IEEE Geophysics paper.

    Kudeki, E.; Milla, M.A., "Incoherent Scatter Spectral Theories—Part I: A General Framework and Results for Small Magnetic Aspect Angles," Geoscience and Remote Sensing, IEEE Transactions on , vol.49, no.1, pp.315,328, Jan. 2011

 Like the model covered in the paper the software can calculate a spectra given, any magnetic aspect angle not perpendicular to B, any number of ion species, and the collision frequencies associated with those ion species. As the magnetic aspect angles get closer to perpendicular to B, usually &lt; 1 degree perp to B, more and more calculations are needed for the Gordeyev to converge.
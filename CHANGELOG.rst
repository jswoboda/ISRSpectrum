======================
ISRSpectrum Change Log
======================

.. current developments

v4.2.0
====================

**Added:**

* Created plasmaline.py to for plasma line spectra.
* Example plasma line code in both notebook and script.
* Added notebooks to docs

**Changed:**

* The standard ISR spectra code will now output the Debye length with each spectrum.
* Added the frequency bin that the plasma line spectra will be in if you have a n channels times m channel resolution FFT.
* The slider notebook now uses temperature ratio



v4.1.0
====================

**Added:**

* Workflows for github that will create deployments automatically.

**Changed:**

* Testing method so we can test against numpys fft frequency vector.
* Clarified docs and how to use plugins.
* Allow for the find_gord_plugs.py to be used as a script.



v4.0.9
====================



v4.0.8
====================



v4.0.7
====================



v4.0.6
====================



v4.0.5
====================



v4.0.4
====================



v4.0.2
====================

**Added:**

* Workflow for anaconda publishing



v4.0.1
====================



v4.0.0
====================

**Added:**

* Plugin framework

* The find_gord_plugs.py bin file to find the plugins

* A script call common_examples.py shows a number of cases where various paramerters are varied and the affects on the spectrum

**Changed:**

* Gordeyev integral calculation moved to plugin directory.

**Removed:**

* mathutils.py removed and materials added to gord_default.py

**Fixed:**

* Collision frequency has been changed to radian's/second in the gord_default plugin.



v3.2.4
====================

**Added:**

* Note to README about default branch change

**Changed:**

* Change default branch to main



v3.2.3
====================

**Added:**

* Applied black tool to main code.
* Added versioneer.
* Added Sphinx.



v3.2.2
====================



v3.2.1
====================

**Fixed:**

* Now importing ioncheck and the amu lookup into main package.



v3.2.0
====================

**Added:**

* Added function to get atomic mass units for ion species.
* Added function to check if valid ion species.



v3.1.0
====================

**Added:**

* Added getspecsimple function so users can avoid using data blocks.

**Changed:**

* Applied black to mathutils.py

**Fixed:**

* String naming convention for ions. Second letters are always lower in the names.



v3.0.4
====================



v3.0.3
====================



v3.0.2
====================



v3.0.1
====================

**Fixed:**

* Fixed test code.



v3.0.0
====================

**Added:**

* Black formatting for main code.

**Changed:**

* API now changed to change name of spectrum generating function.
* Actually properly using the __init__ file now.
* Versioning now coming from __init__ file.

**Deprecated:**

* Python 2 no longer supported.

**Removed:**

* ionlinespec.py



v2.0.2
====================



v2.0.1
====================

**Added:**

* Added rever to update versions.

* <news item>



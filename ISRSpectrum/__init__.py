from .ISRSpectrum import Specinit, ioncheck, getionmass
from .plugins import dirpath as __pluginpath__
from .plugins import gordplugs

from . import _version

__version__ = _version.get_versions()["version"]

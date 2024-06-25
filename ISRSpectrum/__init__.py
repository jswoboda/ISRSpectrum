from .ISRSpectrum import Specinit
from .collision_calc import ioncheck, getionmass, get_collisionfreqs
from .plugins import dirpath as __pluginpath__
from .plugins import gordplugs

from . import _version

__version__ = _version.get_versions()["version"]

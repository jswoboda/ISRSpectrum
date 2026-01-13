from . import _version
from .collision_calc import INFODICT, get_collisionfreqs, getionmass, ioncheck
from .ISRSpectrum import Specinit
from .plasmaline import PLspecinit, get_pl_freqs_default, make_pl_spec_default
from .plugins import dirpath as __pluginpath__
from .plugins import gordplugs

__version__ = _version.get_versions()["version"]

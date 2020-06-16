try:
    from pathlib import Path
    Path().expanduser()
except (ImportError,AttributeError):
    raise ValueError('Please use python 3.')

from ._version import get_versions

from .ISRSpectrum import Specinit

from .mathutils import *

__version__ = get_versions()["version"]

del get_versions

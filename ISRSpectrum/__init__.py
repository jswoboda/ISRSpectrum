try:
    from pathlib import Path
    Path().expanduser()
except (ImportError,AttributeError):
    from pathlib2 import Path

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

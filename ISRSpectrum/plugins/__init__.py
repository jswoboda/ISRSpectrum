# Get current path
from importlib import util
import importlib
from pathlib import Path
from .mathutils import sommerfelderfrep


def load_module(p2):
    """Load the plugin modules

    Parameters
    ----------
    p2 : Path
        File that contains the plugin.

    Returns
    -------
    mod_name : str
        Name of the module

    """
    name = p2.stem
    mod_list = name.split("_")[1:]
    mod_name = "_".join(mod_list)
    spec = util.spec_from_file_location("GordPlug", str(p2))
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return mod_name, module


p1 = Path(__file__).absolute()
dirpath = p1.parent


gordplugs = {}
for fname in dirpath.glob("gord_*.py"):
    # Load load the modles

    mod_name, mod_call = load_module(fname)
    gordplugs[mod_name] = mod_call.GordPlug()

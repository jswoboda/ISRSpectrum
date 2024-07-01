#!python
"""

"""
import argparse
import shutil
import sys
from pathlib import Path
import ISRSpectrum


def parse_command_line(str_input=None):
    """This will parse through the command line arguments

    Function to go through the command line and if given a list of strings all
    also output a namespace object.

    Parameters
    ----------
    str_input : list
        A list of strings or the input from the command line.

    Returns
    -------
    input_args : Namespace
        An object holding the input arguments wrt the variables.
    """
    scriptpath = Path(sys.argv[0])
    scriptname = scriptpath.name

    formatter = argparse.RawDescriptionHelpFormatter(scriptname)
    width = formatter._width
    title = "Gordeyve Plugin Helper"
    shortdesc = "Helps with moving and understanding plugins."
    desc = "\n".join(
        (
            "*" * width,
            "*{0:^{1}}*".format(title, width - 2),
            "*{0:^{1}}*".format("", width - 2),
            "*{0:^{1}}*".format(shortdesc, width - 2),
            "*" * width,
        )
    )
    # desc = "This is the run script for SimVSR."
    # if str_input is None:
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-f",
        "--files",
        dest="files",
        nargs="+",
        default=None,
        help="Plugin files to add to your plugin area.",
        required=False,
    )

    if str_input is None:
        return parser.parse_args()
    return parser.parse_args(str_input)


def findthegord(files=None):
    """Finds the plugins

    Parameters
    ----------
    files : list
        Files that will be moved into the plugin folder.
    """
    pathstr = f"Plugin directory: {str(ISRSpectrum.__pluginpath__)}\n"

    str_list = [pathstr]

    for ikey, imod in ISRSpectrum.gordplugs.items():
        mod_path = ISRSpectrum.__pluginpath__.joinpath("gord_" + ikey + ".py")
        istr = f"Plugin name: {ikey}\nLocation: {str(mod_path)}\n"
        str_list.append(istr)
    print_str = "\n".join(str_list)

    print(print_str)

    if not (files is None):
        for ifilen in files:
            ifilep = Path(ifilen)
            if ifilep.name[:5] != "gord_":
                raise NameError(
                    "File name does not follow the naming convention and will not be found and loaded as a plugin. It was not copied."
                )
            newfile = ISRSpectrum.__pluginpath__.joinpath(ifilep.name)
            shutil.copyfile(ifilep, newfile)
            print(f"{ifilen} copied to {str(newfile)}")


if __name__ == "__main__":
    args_commd = parse_command_line()
    arg_dict = {k: v for k, v in args_commd._get_kwargs() if v is not None}
    findthegord(**arg_dict)

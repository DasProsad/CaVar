"""
startup.py

Startup method for cavar
"""

import sys
import platform
from importlib.metadata import version, PackageNotFoundError
from argparse import Namespace
from textwrap import indent

def get_version(pkg: str) -> str:
    """
    Get package version from cavar
    """
    try:
        return version(pkg)
    except PackageNotFoundError:
        return "not installed"

def show_startup(args: Namespace, prog = "cavar", mods =  ("cyvcf2", "pysam")) -> None:
    """
    Print start-up banner with package versions
    """
    try:
        progv = version(prog)
    except PackageNotFoundError:
        progv = "unknown"
    print(f"{prog} {progv} starting with:")
    print(f"  Python    : {sys.version.split()[0]}")
    print(f"  Platform  : {platform.system()} {platform.release()}")
    print()
    print("External dependencies:")
    for mod in mods:
        print(f"  {mod:<10}: {get_version(mod)}")
    print()
    print("Arguments:")
    for k, v in vars(args).items():
        print(f"  {k:<10}: {v}")
    print()

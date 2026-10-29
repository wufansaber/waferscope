"""Locate portable solver assets without depending on the development checkout."""
from pathlib import Path
import sys


def distribution_root():
    if not getattr(sys, 'frozen', False):
        return Path(__file__).resolve().parents[1]
    folder = Path(sys.executable).resolve().parent
    return folder.parent if folder.name.lower() == 'fem' else folder


def default_engine(root):
    bundled = Path(root)/'fem/WaferScope_FEM.exe'
    return bundled if bundled.is_file() else Path(root)/'build/rcwa_env/Scripts/python.exe'


def default_s4(root=None):
    root = Path(root) if root is not None else distribution_root()
    bundled = root/'s4/S4.exe'
    return bundled if bundled.is_file() else root/'build/s4_runtime/bin_stripes/S4.exe'

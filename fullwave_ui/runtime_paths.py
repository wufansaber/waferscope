"""Paths shared by source runs and the two-executable Windows distribution."""
import os
from pathlib import Path
import sys


def application_root():
    return (Path(sys.executable).resolve().parent if getattr(sys, 'frozen', False)
            else Path(__file__).resolve().parents[1])


def default_engine(root):
    bundled = root/'fem/FarField_FEM.exe'
    return bundled if bundled.is_file() else root/'build/ngsolve_env/Scripts/python.exe'


def results_root(root):
    override = os.environ.get('FARFIELD_RESULTS_DIR')
    if override:
        return Path(override).expanduser().resolve()
    if getattr(sys, 'frozen', False):
        return Path(os.environ.get('LOCALAPPDATA', str(Path.home())))/'FarField_Simulation/results'
    return root/'build/fullwave_ui'

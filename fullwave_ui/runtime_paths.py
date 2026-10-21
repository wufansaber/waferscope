"""Paths shared by source runs and the two-executable Windows distribution."""
import os
from pathlib import Path
import sys


def application_root():
    return (Path(sys.executable).resolve().parent if getattr(sys, 'frozen', False)
            else Path(__file__).resolve().parents[1])


def default_engine(root):
    for name in ('WaferScope_FEM.exe', 'DFI_FEM.exe', 'FarField_FEM.exe'):
        bundled = root/'fem'/name
        if bundled.is_file():
            return bundled
    return root/'build/ngsolve_env/Scripts/python.exe'


def results_root(root):
    override = os.environ.get('FARFIELD_RESULTS_DIR')
    if override:
        return Path(override).expanduser().resolve()
    if getattr(sys, 'frozen', False):
        return Path(os.environ.get('LOCALAPPDATA', str(Path.home())))/'FarField_Simulation/results'
    return root/'build/fullwave_ui'

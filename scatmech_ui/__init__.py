"""Modular SCATMECH particle-on-wafer simulation interface."""

from .app import ScatmechStandaloneWindow, main, run_headless_test
from .backend import SCATMECH_IMPORT_ERROR, calculate_far_field
from .materials import PARTICLE_MATERIALS, WAFER_MATERIALS
from .parameters import SimulationParameters
from .plotting import PARTICLE_INTENSITY_CMAP, PlotCanvas
from .widgets import MaterialGroup, ScatmechSimulationWidget
from .worker import SimulationWorker

__all__ = [
    "MaterialGroup",
    "PARTICLE_INTENSITY_CMAP",
    "PARTICLE_MATERIALS",
    "PlotCanvas",
    "SCATMECH_IMPORT_ERROR",
    "ScatmechSimulationWidget",
    "ScatmechStandaloneWindow",
    "SimulationParameters",
    "SimulationWorker",
    "WAFER_MATERIALS",
    "calculate_far_field",
    "main",
    "run_headless_test",
]

"""Modular SCATMECH particle-on-wafer simulation interface."""

from .app import ScatmechStandaloneWindow, main, run_headless_test
from .backend import SCATMECH_IMPORT_ERROR, calculate_far_field
from .materials import COATING_MATERIALS, PARTICLE_MATERIALS, WAFER_MATERIALS
from .parameters import (
    CoatingLayerParameters,
    ParticleSizeSweepParameters,
    SimulationParameters,
)
from .plotting import PARTICLE_INTENSITY_CMAP, PlotCanvas
from .size_sweep import calculate_particle_size_sweep
from .widgets import MaterialGroup, ScatmechSimulationWidget
from .worker import ParticleSizeSweepWorker, SimulationWorker

__all__ = [
    "MaterialGroup",
    "COATING_MATERIALS",
    "CoatingLayerParameters",
    "PARTICLE_INTENSITY_CMAP",
    "PARTICLE_MATERIALS",
    "ParticleSizeSweepParameters",
    "ParticleSizeSweepWorker",
    "PlotCanvas",
    "SCATMECH_IMPORT_ERROR",
    "ScatmechSimulationWidget",
    "ScatmechStandaloneWindow",
    "SimulationParameters",
    "SimulationWorker",
    "WAFER_MATERIALS",
    "calculate_far_field",
    "calculate_particle_size_sweep",
    "main",
    "run_headless_test",
]

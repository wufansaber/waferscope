"""Compatibility entry point for the modular :mod:`scatmech_ui` package."""

from scatmech_ui import (
    MaterialGroup,
    PARTICLE_INTENSITY_CMAP,
    PARTICLE_MATERIALS,
    PlotCanvas,
    SCATMECH_IMPORT_ERROR,
    ScatmechSimulationWidget,
    ScatmechStandaloneWindow,
    SimulationParameters,
    SimulationWorker,
    WAFER_MATERIALS,
    calculate_far_field,
    main,
    run_headless_test,
)

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


if __name__ == "__main__":
    main()

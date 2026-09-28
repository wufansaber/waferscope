"""Shared plotting facade for far-field and objective-collection views."""

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure

from .far_field_plotting import plot_2d_far_field, plot_3d_far_field
from .objective_plotting import (
    plot_objective_collection_2d,
    plot_objective_collection_3d,
)
from .size_sweep_plotting import plot_particle_size_sweep
from .imaging_plotting import plot_imaging
from .detector_plotting import plot_detectors


PARTICLE_INTENSITY_CMAP = LinearSegmentedColormap.from_list(
    "particle_intensity",
    [
        (0.00, "#050b18"),
        (0.08, "#151746"),
        (0.18, "#d7191c"),
        (0.48, "#ff6d00"),
        (0.82, "#fff21a"),
        (1.00, "#fffde8"),
    ],
)


class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self):
        self.figure = Figure(figsize=(8, 6), constrained_layout=True)
        super().__init__(self.figure)


class ImagingCanvas(PlotCanvas):
    def __init__(self):
        self.imaging_result = None
        self.imaging_narrow = None
        super().__init__()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.imaging_result is not None and (self.width() < 950) != self.imaging_narrow:
            plot_imaging(self, self.imaging_result)


class DetectorCanvas(PlotCanvas):
    def __init__(self):
        self.detector_result = None
        self.detector_narrow = None
        super().__init__()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.detector_result is not None and (self.width() < 950) != self.detector_narrow:
            plot_detectors(self, self.detector_result)


def render_results(canvases, result):
    if result.get("plot_kind") == "coherent_imaging":
        plot_imaging(canvases["image_planes"], result)
        if "detectors" in canvases:
            plot_detectors(canvases["detectors"], result)
        return
    if result.get("plot_kind") == "particle_size_sweep":
        plot_particle_size_sweep(canvases["size_sweep"], result)
        return
    plot_3d_intensity(canvases["3d_intensity"], result)
    plot_2d_intensity(canvases["2d_intensity"], result)


def plot_3d_intensity(canvas, result):
    if result.get("plot_kind") == "objective_collection":
        plot_objective_collection_3d(canvas, result)
        return
    plot_3d_far_field(canvas, result)


def plot_2d_intensity(canvas, result):
    if result.get("plot_kind") == "objective_collection":
        plot_objective_collection_2d(canvas, result)
        return
    plot_2d_far_field(canvas, result, PARTICLE_INTENSITY_CMAP)

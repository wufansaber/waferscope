"""Shared plotting facade for far-field and objective-collection views."""

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure

from .far_field_plotting import plot_2d_far_field, plot_3d_far_field
from .objective_plotting import (
    plot_objective_collection_2d,
    plot_objective_collection_3d,
)


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


def render_results(canvases, result):
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

"""Matplotlib canvases and far-field visualization functions."""

import matplotlib.cm as cm
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Line3DCollection


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


def _sparse_indices(result):
    rows = np.arange(
        2,
        result["x"].shape[0] - 1,
        max(1, result["x"].shape[0] // 9),
    )
    columns = np.arange(
        0,
        result["x"].shape[1] - 1,
        max(1, result["x"].shape[1] // 18),
    )
    return np.ix_(rows, columns)


def plot_3d_intensity(canvas, result):
    canvas.figure.clear()
    axis = canvas.figure.add_subplot(111, projection="3d")
    floor = result["parameters"].display_floor_db
    norm = Normalize(floor, 0.0)
    colors = cm.turbo(norm(np.clip(result["intensity_db"], floor, 0.0)))
    axis.plot_surface(
        result["x"],
        result["y"],
        result["z"],
        facecolors=colors,
        linewidth=0,
        antialiased=True,
        shade=False,
    )
    _add_3d_polarization_axes(axis, result)
    _draw_3d_incidence_geometry(axis, result)
    _format_3d_axis(axis, "Far-field intensity with polarization directions")
    colorbar = canvas.figure.colorbar(
        cm.ScalarMappable(norm=norm, cmap=cm.turbo),
        ax=axis,
        shrink=0.70,
    )
    colorbar.set_label("Relative DSC (dB)")
    canvas.draw_idle()


def _add_3d_polarization_axes(axis, result):
    """Overlay Stokes polarization directions on the intensity surface."""
    index = _sparse_indices(result)
    length = 0.12
    segments = []
    for x, y, z, ex, ey, ez, dolp in zip(
        result["x"][index].ravel(),
        result["y"][index].ravel(),
        result["z"][index].ravel(),
        result["e_x"][index].ravel(),
        result["e_y"][index].ravel(),
        result["e_z"][index].ravel(),
        result["dolp"][index].ravel(),
    ):
        half = 0.5 * length * max(0.25, dolp)
        segments.append(
            [
                (x - half * ex, y - half * ey, z - half * ez),
                (x + half * ex, y + half * ey, z + half * ez),
            ]
        )
    axis.add_collection3d(
        Line3DCollection(segments, colors="#1236b5", linewidths=1.4)
    )


def plot_2d_intensity(canvas, result):
    canvas.figure.clear()
    axis = canvas.figure.add_subplot(111)
    relative = result["intensity"] / np.max(result["intensity"])
    display_intensity = np.clip(relative, 0.0, 1.0) ** 0.80
    image = axis.pcolormesh(
        result["x"],
        result["y"],
        display_intensity,
        cmap=PARTICLE_INTENSITY_CMAP,
        vmin=0.0,
        vmax=1.0,
        shading="gouraud",
    )
    _add_2d_polarization_axes(axis, result)
    _format_2d_axis(
        axis,
        "Particle Scattering\nIntensity with Polarization Directions",
    )
    colorbar = canvas.figure.colorbar(image, ax=axis)
    colorbar.set_label("Normalized particle scattering intensity")
    canvas.draw_idle()


def _add_2d_polarization_axes(axis, result):
    """Overlay projected Stokes polarization axes on a pupil plot."""
    index = _sparse_indices(result)
    px = result["x"][index].ravel()
    py = result["y"][index].ravel()
    ex = result["e_x"][index].ravel()
    ey = result["e_y"][index].ravel()
    dolp = result["dolp"][index].ravel()
    projected_norm = np.hypot(ex, ey)
    valid = projected_norm > 1e-8
    segments = []
    for x, y, vx, vy, degree, magnitude in zip(
        px[valid],
        py[valid],
        ex[valid],
        ey[valid],
        dolp[valid],
        projected_norm[valid],
    ):
        half = 0.06 * max(0.25, degree)
        vx, vy = vx / magnitude, vy / magnitude
        segments.append(
            [(x - half * vx, y - half * vy), (x + half * vx, y + half * vy)]
        )
    polarization_lines = LineCollection(
        segments, colors="#3157d5", linewidths=1.45
    )
    polarization_lines.set_clip_path(
        Circle((0.0, 0.0), 0.995, transform=axis.transData)
    )
    axis.add_collection(polarization_lines)


def _draw_3d_incidence_geometry(axis, result):
    """Draw wafer, surface normal, laser ray, and incidence-angle arc."""
    incidence_deg = result["parameters"].incidence_deg
    incidence = np.deg2rad(incidence_deg)
    diagram_y = 0.0

    plane_x, plane_y = np.meshgrid(
        np.linspace(-1.45, 1.10, 2),
        np.linspace(-1.10, 1.10, 2),
    )
    axis.plot_surface(
        plane_x,
        plane_y,
        np.zeros_like(plane_x),
        color="#8d9499",
        alpha=0.16,
        shade=False,
        zorder=0,
    )
    axis.text(0.90, -0.95, 0.015, "wafer plane", color="#50565a")

    axis.quiver(
        0.0,
        diagram_y,
        0.0,
        0.0,
        0.0,
        0.72,
        color="#666666",
        linewidth=2.0,
        arrow_length_ratio=0.10,
        zorder=20,
    )
    axis.text(
        0.03,
        diagram_y,
        0.76,
        "surface normal",
        color="#4f5559",
        zorder=20,
    )

    source_distance = 1.40
    source_x = -source_distance * np.sin(incidence)
    source_z = source_distance * np.cos(incidence)
    axis.quiver(
        source_x,
        diagram_y,
        source_z,
        -source_x,
        0.0,
        -source_z,
        color="#133e9a",
        linewidth=3.2,
        arrow_length_ratio=0.10,
        zorder=20,
    )
    axis.text(
        source_x - 0.02,
        diagram_y,
        source_z + 0.07,
        "incident laser",
        color="#133e9a",
        weight="bold",
        ha="center",
        zorder=20,
    )
    axis.scatter(
        [0.0],
        [0.0],
        [0.0],
        color="#133e9a",
        s=28,
        depthshade=False,
        zorder=21,
    )
    axis.text(
        0.04,
        0.02,
        0.03,
        "incidence point (0, 0, 0)",
        color="#133e9a",
        fontsize=8,
        zorder=21,
    )

    arc_angle = np.linspace(0.0, incidence, 80)
    arc_radius = 0.30
    arc_x = -arc_radius * np.sin(arc_angle)
    arc_z = arc_radius * np.cos(arc_angle)
    axis.plot(
        arc_x,
        np.full_like(arc_x, diagram_y),
        arc_z,
        color="#133e9a",
        linewidth=1.8,
        zorder=20,
    )
    label_angle = incidence * 0.52
    axis.text(
        -0.38 * np.sin(label_angle),
        diagram_y,
        0.38 * np.cos(label_angle),
        "{:.1f} deg".format(incidence_deg),
        color="#133e9a",
        ha="center",
        va="center",
        zorder=20,
    )


def _format_3d_axis(axis, title):
    axis.set(
        xlim=(-1.50, 1.10),
        ylim=(-1.25, 1.10),
        zlim=(0.0, 1.12),
        xlabel=r"$k_x/k$",
        ylabel=r"$k_y/k$",
        zlabel=r"$k_z/k$",
    )
    axis.set_title(title, fontsize=11, pad=4)
    axis.tick_params(labelsize=8, pad=1)
    axis.xaxis.label.set_size(9)
    axis.yaxis.label.set_size(9)
    axis.zaxis.label.set_size(9)
    axis.set_box_aspect((1, 1, 0.7))
    axis.view_init(elev=28, azim=-55)


def _format_2d_axis(axis, title):
    horizon = np.linspace(0.0, 2.0 * np.pi, 361)
    axis.plot(np.cos(horizon), np.sin(horizon), color="white", linewidth=1.0)
    axis.text(
        -0.97,
        0.94,
        "Particle",
        color="white",
        fontsize=13,
        weight="bold",
        ha="left",
        va="top",
        zorder=20,
    )
    axis.annotate(
        "",
        xy=(0.90, -1.20),
        xytext=(-0.90, -1.20),
        arrowprops=dict(arrowstyle="->", color="#123b91", lw=2.8),
        annotation_clip=False,
    )
    axis.text(
        0.0,
        -1.13,
        "Laser Incidence Direction",
        color="#123b91",
        fontsize=10,
        weight="bold",
        ha="center",
        va="bottom",
    )
    axis.set(
        xlim=(-1.05, 1.05),
        ylim=(-1.30, 1.05),
        aspect="equal",
        xlabel="",
        ylabel=r"$k_y/k$",
    )
    axis.set_title(title, fontsize=11, pad=4)
    axis.tick_params(labelsize=8, pad=1)
    axis.yaxis.label.set_size(9)
    axis.set_facecolor("#111111")

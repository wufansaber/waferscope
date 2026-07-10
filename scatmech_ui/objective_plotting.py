"""Matplotlib plots for ideal-objective collection results."""

import matplotlib.cm as cm
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Circle


def plot_objective_collection_3d(canvas, result):
    canvas.figure.clear()
    axis = canvas.figure.add_subplot(111, projection="3d")
    floor = result["parameters"].display_floor_db
    norm = Normalize(floor, 0.0)

    _draw_far_field_surface(axis, result, norm, floor)
    _draw_collected_surface(axis, result, norm, floor)
    _draw_collection_cones(axis, result)
    _format_collection_3d_axis(axis, result)

    colorbar = canvas.figure.colorbar(
        cm.ScalarMappable(norm=norm, cmap=cm.turbo),
        ax=axis,
        shrink=0.70,
    )
    colorbar.set_label("Far-field energy vs fixed reference (dB)")
    canvas.draw_idle()


def plot_objective_collection_2d(canvas, result):
    canvas.figure.clear()
    floor = result["parameters"].display_floor_db
    radius = result["parameters"].focal_length_mm * result["parameters"].numerical_aperture
    axes = canvas.figure.subplots(1, 3, sharex=True, sharey=True)
    images = [
        _plot_channel_pupil(axis, result, pupil, index, radius, floor)
        for index, (axis, pupil) in enumerate(
            zip(axes, result["channel_pupils"]),
            start=1,
        )
    ]
    canvas.figure.suptitle(
        "Three-channel collected pupil energy, union={:.2f}%".format(
            100.0 * result["collection_fraction"]
        ),
        fontsize=11,
    )
    colorbar = canvas.figure.colorbar(images[-1], ax=list(axes), shrink=0.82)
    colorbar.set_label("Collected energy vs fixed reference (dB)")
    canvas.draw_idle()


def _draw_far_field_surface(axis, result, norm, floor):
    colors = cm.turbo(norm(np.clip(result["far_intensity_db"], floor, 0.0)))
    colors[..., 3] = 0.24
    axis.plot_surface(
        result["far_x"],
        result["far_y"],
        result["far_z"],
        facecolors=colors,
        linewidth=0,
        antialiased=True,
        shade=False,
    )


def _draw_collected_surface(axis, result, norm, floor):
    collected_z = np.where(
        result["far_collection_mask"],
        result["far_z"],
        np.nan,
    )
    colors = cm.turbo(norm(np.clip(result["far_intensity_db"], floor, 0.0)))
    colors[..., 3] = 0.94
    axis.plot_surface(
        result["far_x"],
        result["far_y"],
        collected_z,
        facecolors=colors,
        linewidth=0,
        antialiased=True,
        shade=False,
    )


def _draw_collection_cones(axis, result):
    channel_colors = ("#101820", "#2457c5", "#13795b")
    for channel_index, (objective_axis, boundary, color) in enumerate(
        zip(result["objective_axes"], result["cone_boundaries"], channel_colors),
        start=1,
    ):
        _draw_cone_boundary(axis, boundary, color)
        _draw_objective_axis(axis, objective_axis, channel_index, color)

    horizon = np.linspace(0.0, 2.0 * np.pi, 361)
    axis.plot(
        np.cos(horizon),
        np.sin(horizon),
        np.zeros_like(horizon),
        color="#666666",
        linewidth=0.8,
        alpha=0.5,
    )


def _draw_cone_boundary(axis, boundary, color):
    boundary_x, boundary_y, boundary_z = boundary
    axis.plot(boundary_x, boundary_y, boundary_z, color=color, linewidth=2.1)
    for index in np.linspace(0, len(boundary_x) - 1, 10, dtype=int):
        axis.plot(
            [0.0, boundary_x[index]],
            [0.0, boundary_y[index]],
            [0.0, boundary_z[index]],
            color=color,
            alpha=0.16,
            linewidth=0.8,
        )


def _draw_objective_axis(axis, objective_axis, channel_index, color):
    axis.quiver(
        0.0,
        0.0,
        0.0,
        objective_axis[0],
        objective_axis[1],
        objective_axis[2],
        color=color,
        linewidth=2.4,
        arrow_length_ratio=0.08,
    )
    axis.text(
        1.06 * objective_axis[0],
        1.06 * objective_axis[1],
        1.06 * objective_axis[2],
        "CH{}".format(channel_index),
        color=color,
        weight="bold",
        ha="center",
    )


def _format_collection_3d_axis(axis, result):
    axis.set(
        xlim=(-1.05, 1.05),
        ylim=(-1.05, 1.05),
        zlim=(0.0, 1.08),
        xlabel=r"$k_x/k$",
        ylabel=r"$k_y/k$",
        zlabel=r"$k_z/k$",
    )
    axis.set_title(
        "Three-channel far-field collection cones\n"
        "NA={:.3g}, f={:.3g} mm, tilts={}, union collected={:.2f}%".format(
            result["parameters"].numerical_aperture,
            result["parameters"].focal_length_mm,
            "/".join(
                "{:.3g}".format(tilt)
                for tilt in result["parameters"].tilt_degrees
            ),
            100.0 * result["collection_fraction"],
        ),
        fontsize=11,
        pad=4,
    )
    axis.text2D(
        0.03,
        0.93,
        "Union cone / hemisphere energy: {:.2f}%".format(
            100.0 * result["collection_fraction"]
        ),
        transform=axis.transAxes,
        color="#101820",
        fontsize=10,
        weight="bold",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="#101820",
            alpha=0.82,
        ),
    )
    axis.tick_params(labelsize=8, pad=1)
    axis.set_box_aspect((1, 1, 0.75))
    axis.view_init(elev=28, azim=-52)


def _plot_channel_pupil(axis, result, pupil, index, radius, floor):
    display_intensity = np.ma.masked_invalid(
        np.clip(pupil["intensity_db"], floor, 0.0)
    )
    image = axis.pcolormesh(
        pupil["x"],
        pupil["y"],
        display_intensity,
        cmap=cm.turbo,
        vmin=floor,
        vmax=0.0,
        shading="gouraud",
    )
    _draw_pupil_outline(axis, radius)
    _add_objective_pupil_polarization(axis, pupil, radius)
    _format_channel_axis(axis, result, pupil, index, radius)
    return image


def _draw_pupil_outline(axis, radius):
    horizon = np.linspace(0.0, 2.0 * np.pi, 361)
    axis.plot(
        radius * np.cos(horizon),
        radius * np.sin(horizon),
        color="white",
        linewidth=1.1,
    )


def _format_channel_axis(axis, result, pupil, index, radius):
    channel_fraction = result["channel_fractions"][index - 1]
    mean_dolp = 100.0 * np.nanmean(pupil["dolp"])
    axis.set(
        xlim=(-1.05 * radius, 1.05 * radius),
        ylim=(-1.05 * radius, 1.05 * radius),
        aspect="equal",
        xlabel="Exit pupil x (mm)",
    )
    if index == 1:
        axis.set_ylabel("Exit pupil y (mm)")
    axis.set_title(
        "CH{} tilt={:.3g} deg\nenergy={:.2f}%, DoLP={:.1f}%".format(
            index,
            pupil["tilt_deg"],
            100.0 * channel_fraction,
            mean_dolp,
        ),
        fontsize=10,
        pad=4,
    )
    axis.tick_params(labelsize=7, pad=1)
    axis.set_facecolor("#111111")


def _add_objective_pupil_polarization(axis, pupil, radius):
    step = max(10, pupil["x"].shape[0] // 12)
    rows = np.arange(step // 2, pupil["x"].shape[0], step)
    columns = np.arange(step // 2, pupil["x"].shape[1], step)
    index = np.ix_(rows, columns)
    px = pupil["x"][index].ravel()
    py = pupil["y"][index].ravel()
    vx = pupil["pol_x"][index].ravel()
    vy = pupil["pol_y"][index].ravel()
    dolp = pupil["dolp"][index].ravel()
    valid = np.isfinite(dolp)
    segments = []
    for x, y, ex, ey, degree in zip(
        px[valid],
        py[valid],
        vx[valid],
        vy[valid],
        dolp[valid],
    ):
        half = 0.055 * radius * max(0.25, degree)
        segments.append(
            [(x - half * ex, y - half * ey), (x + half * ex, y + half * ey)]
        )
    polarization_lines = LineCollection(
        segments,
        colors="#f4f7ff",
        linewidths=1.05,
        alpha=0.88,
    )
    polarization_lines.set_clip_path(
        Circle((0.0, 0.0), radius, transform=axis.transData)
    )
    axis.add_collection(polarization_lines)

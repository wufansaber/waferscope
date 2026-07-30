"""Matplotlib chart for particle-size sweep results."""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from .coordinate_geometry import incidence_projection_text


CHANNEL_STYLES = (
    ("#1769aa", "o", "-"),
    ("#d18b00", "s", "--"),
    ("#c23b70", "^", "-."),
)


def plot_particle_size_sweep(canvas, result):
    canvas.figure.clear()
    intensity_axis, fraction_axis = canvas.figure.subplots(
        2,
        1,
        sharex=True,
        gridspec_kw={"height_ratios": (2.2, 1.0)},
    )
    diameters = result["diameters_nm"]
    collected = result["channel_collected_dsc_um2"]
    fractions = 100.0 * result["channel_energy_fractions"]
    objective = result["objective_parameters"]
    mark_every = max(1, len(diameters) // 18)
    legend_handles = []

    for index, (color, marker, linestyle) in enumerate(CHANNEL_STYLES):
        label = "CH{}  NA={:.3g}, tilt={:+.3g} deg".format(
            index + 1,
            objective.channel_numerical_apertures[index],
            objective.tilt_degrees[index],
        )
        if np.any(collected[:, index] > 0.0):
            intensity_axis.plot(
                diameters,
                collected[:, index],
                color=color,
                marker=marker,
                linestyle=linestyle,
                linewidth=2.0,
                markersize=5.0,
                markerfacecolor="white",
                markeredgewidth=1.2,
                markevery=mark_every,
                label=label,
            )
            fraction_axis.plot(
                diameters,
                fractions[:, index],
                color=color,
                marker=marker,
                linestyle=linestyle,
                linewidth=1.8,
                markersize=4.5,
                markerfacecolor="white",
                markeredgewidth=1.1,
                markevery=mark_every,
            )
        else:
            label += " (disabled)"
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                marker=marker,
                linestyle=linestyle,
                markerfacecolor="white",
                label=label,
            )
        )

    _format_sweep_axes(
        canvas,
        result,
        intensity_axis,
        fraction_axis,
        legend_handles,
    )
    canvas.draw_idle()


def _format_sweep_axes(canvas, result, intensity_axis, fraction_axis, handles):
    parameters = result["simulation_parameters"]
    diameters = result["diameters_nm"]
    canvas.figure.suptitle(
        "Particle-size sweep: three-channel collected scattering",
        fontsize=12,
    )
    intensity_axis.set_title(
        "{} nm, {} particle, {} wafer, {:.3g} deg {} polarization, "
        "illumination {}".format(
            parameters.wavelength_nm,
            parameters.particle_material,
            parameters.wafer_material,
            parameters.incidence_deg,
            parameters.polarization.upper(),
            incidence_projection_text(parameters.incidence_deg),
        ),
        fontsize=9,
        color="#4d5560",
        pad=6,
    )
    if np.any(result["channel_collected_dsc_um2"] > 0.0):
        intensity_axis.set_yscale("log")
    else:
        intensity_axis.text(
            0.5,
            0.5,
            "All collection channels are disabled",
            transform=intensity_axis.transAxes,
            ha="center",
            va="center",
        )
    intensity_axis.set_ylabel("Integrated collected DSC (um^2)")
    fraction_axis.set_ylabel("Energy fraction (%)")
    fraction_axis.set_xlabel("Particle diameter (nm)")
    intensity_axis.legend(handles=handles, frameon=False, loc="upper left", fontsize=8)

    for axis in (intensity_axis, fraction_axis):
        axis.grid(True, which="major", color="#d9dde3", linewidth=0.8)
        axis.grid(True, which="minor", axis="y", color="#eceff3", linewidth=0.5)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.set_xlim(float(diameters[0]), float(diameters[-1]))
    if len(diameters) <= 20:
        fraction_axis.set_xticks(diameters)
        fraction_axis.tick_params(axis="x", rotation=45, labelsize=8)
    else:
        fraction_axis.xaxis.set_major_locator(MaxNLocator(nbins=12))

"""Absolute, independent scales for three detector planes."""

import numpy as np
from matplotlib.colors import PowerNorm


def plot_imaging(canvas, result):
    figure = canvas.figure
    figure.clear()
    narrow = canvas.width() < 950
    canvas.imaging_result = result
    canvas.imaging_narrow = narrow
    axes = figure.subplots(3, 2, squeeze=False) if narrow else figure.subplots(2, 3, squeeze=False)
    imaging = result["imaging_parameters"]
    figure.suptitle("Independent images | M={:g}x | tube={:g} mm | defocus={:g} um".format(
        imaging.magnification, result["tube_focal_length_mm"], imaging.image_defocus_um), fontsize=11)
    for index, channel in enumerate(result["channels"]):
        axis = channel["image_axis_um"]
        values = channel["irradiance_w_um2"]
        dx = axis[1] - axis[0]
        extent = (axis[0]-dx/2, axis[-1]+dx/2, axis[0]-dx/2, axis[-1]+dx/2)
        ax = axes[index, 0] if narrow else axes[0, index]
        vmax = max(float(values.max()), 1e-30)
        im = ax.imshow(values, origin="lower", extent=extent, cmap="inferno",
                       norm=PowerNorm(0.35, vmin=0, vmax=vmax))
        # Show the central PSF while retaining the full computational array in exports.
        radius = min(axis[-1], 4*result["parameters"].wavelength_nm/1000
                     * imaging.magnification / max(channel["na"], 0.1))
        ax.set_xlim(-radius, radius)
        ax.set_ylim(-radius, radius)
        ax.set_xlabel("Image x (um)")
        ax.set_ylabel("Image y (um)")
        label = "OFF" if channel["na"] == 0 else "NA={:g}, tilt={:g} deg".format(
            channel["na"], channel["tilt_deg"])
        ax.set_title("CH{}  {}\nP={:.4g} W".format(index+1, label, channel["power_w"]), fontsize=8 if narrow else 9)
        figure.colorbar(im, ax=ax, shrink=0.7, label="W/um^2", fraction=0.045)
        ax = axes[index, 1] if narrow else axes[1, index]
        ax.plot(axis, values[len(axis)//2], color="#007f73")
        ax.set_xlim(-radius, radius)
        ax.set_xlabel("Image x (um), y=0")
        ax.set_ylabel("W/um^2")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.grid(alpha=0.2)
        ax.set_title("Flux error {:.1e}\nEdge {:.2%}{}".format(
            channel["flux_relative_error"], channel["edge_fraction"],
            "\nSpecular direction admitted" if channel["admits_specular"] else ""), fontsize=9)
    canvas.draw_idle()

"""Absolute expected and sampled detector electron counts."""

import numpy as np
from matplotlib.colors import PowerNorm


def plot_detectors(canvas, result):
    canvas.detector_result = result
    narrow = canvas.width() < 950
    canvas.detector_narrow = narrow
    figure = canvas.figure
    figure.clear()
    if not result["imaging_parameters"].detector.enabled:
        figure.text(0.5, 0.5, "Detectors disabled", ha="center", va="center")
        canvas.draw_idle()
        return
    settings = result["imaging_parameters"].detector
    axes = figure.subplots(3, 2) if narrow else figure.subplots(2, 3)
    for i, channel in enumerate(result["channels"]):
        detector = channel["detector"]
        x, y = detector["pixel_x_um"], detector["pixel_y_um"]
        half = settings.pixel_pitch_um/2
        extent = (x[0]-half, x[-1]+half, y[0]-half, y[-1]+half)
        for j, (key, title) in enumerate((("mean_signal_e", "Expected particle signal"),
                                          ("digital_dn", "Raw camera output"))):
            ax = axes[i, j] if narrow else axes[j, i]
            data = detector[key]
            # Retain signed readout noise; independent absolute scales per view.
            vmin = min(0.0, float(data.min()))
            vmax = max(1.0, float(data.max()))
            artist = ax.imshow(data, origin="lower", extent=extent, interpolation="nearest",
                               cmap="inferno" if j == 0 else "gray", norm=PowerNorm(0.5, vmin=vmin, vmax=vmax))
            ax.set_xlabel("Detector X (um)", fontsize=8)
            ax.set_ylabel("Detector Y (um)", fontsize=8)
            ax.tick_params(labelsize=7)
            if j == 0:
                detail = "{:.4g} e-; ROI {:.2%}".format(detector["signal_e"], detector["roi_power_fraction"])
            else:
                detail = "Total SNR {:.3g}; well/ADC sat {}/{} px".format(
                    detector["emva_total_roi_snr"], detector["saturated_pixels"], detector['adc_saturated_pixels'])
            warning = j == 1 and (detector["saturated_pixels"] or detector["expected_overfull_pixels"]
                                  or detector['adc_saturated_pixels'] or not detector['digital_snr_valid'].all())
            ax.set_title("CH{} | {}\n{}".format(i+1, title, detail), fontsize=9,
                         color="firebrick" if warning else "black")
            ticks = np.linspace(vmin, vmax, 5)
            colorbar = figure.colorbar(artist, ax=ax, pad=0.02, fraction=0.046,
                                      ticks=ticks, format="%.3g")
            colorbar.set_label("e- / pixel" if j == 0 else "DN / pixel", fontsize=8)
            colorbar.ax.tick_params(labelsize=7)
    figure.suptitle("Independent cameras | {:.4g} ms | QE {:.3g} | {:.4g} um | {} bit, {:.3g} e-/DN".format(
        settings.exposure_ms, settings.quantum_efficiency, settings.pixel_pitch_um,
        settings.bit_depth, settings.gain_e_dn), fontsize=10)
    canvas.draw_idle()

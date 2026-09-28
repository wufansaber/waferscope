"""Backend-independent area-integrating detector ROIs for sampled irradiance."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DetectorParameters:
    enabled: bool = False
    pixel_pitch_um: float = 6.5
    pixels_x: int = 32
    pixels_y: int = 32
    exposure_ms: float = 1.0
    quantum_efficiency: float = 0.6
    read_noise_e: float = 2.0
    dark_current_e_s: float = 0.1
    background_e_s: float = 0.0
    full_well_e: float = 30000.0
    seed: int = 0

    def validate(self):
        positive = (self.pixel_pitch_um, self.exposure_ms, self.full_well_e)
        nonnegative = (self.read_noise_e, self.dark_current_e_s, self.background_e_s)
        if not np.isfinite(positive + nonnegative + (self.quantum_efficiency,)).all():
            raise ValueError("Detector parameters must be finite.")
        if min(positive) <= 0 or min(nonnegative) < 0:
            raise ValueError("Pixel pitch, exposure and full well must be positive; noise rates nonnegative.")
        if not 0 <= self.quantum_efficiency <= 1:
            raise ValueError("Quantum efficiency must be between 0 and 1.")
        for count in (self.pixels_x, self.pixels_y):
            if not isinstance(count, int) or isinstance(count, bool) or not 1 <= count <= 512:
                raise ValueError("Detector ROI dimensions must be integers from 1 to 512.")
        if not isinstance(self.seed, int) or not 0 <= self.seed <= 2**32 - 1:
            raise ValueError("Noise seed must be an integer from 0 to 4294967295.")


def _overlap_weights(axis, count, pitch):
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or len(axis) < 2 or not np.isfinite(axis).all():
        raise ValueError("Image axis must be a finite, uniformly spaced array.")
    step = axis[1] - axis[0]
    if step <= 0 or not np.allclose(np.diff(axis), step, rtol=1e-9, atol=1e-12):
        raise ValueError("Image axis must be uniformly increasing.")
    edges = (np.arange(count + 1) - count / 2) * pitch
    tolerance = step * 1e-8
    if edges[0] < axis[0] - step/2 - tolerance or edges[-1] > axis[-1] + step/2 + tolerance:
        raise ValueError("Detector ROI exceeds the computed image field. Increase pupil samples "
                         "or reduce ROI size/pixel pitch; FFT padding does not enlarge the field.")
    overlap = np.maximum(0, np.minimum(edges[1:, None], axis[None, :] + step/2)
                         - np.maximum(edges[:-1, None], axis[None, :] - step/2))
    return (edges[:-1] + edges[1:])/2, overlap


def integrate_pixels(axis, irradiance, parameters):
    """Integrate a piecewise-constant image grid, without wrapping or extrapolation."""
    parameters.validate()
    irradiance = np.asarray(irradiance, dtype=float)
    if irradiance.shape != (len(axis), len(axis)) or not np.isfinite(irradiance).all():
        raise ValueError("Image irradiance must be a finite square array matching the image axis.")
    peak = np.max(np.abs(irradiance))
    if np.min(irradiance) < -1e-10 * peak:
        raise ValueError("Image has negative axial flux; this detector model requires nonnegative irradiance.")
    x, wx = _overlap_weights(axis, parameters.pixels_x, parameters.pixel_pitch_um)
    y, wy = _overlap_weights(axis, parameters.pixels_y, parameters.pixel_pitch_um)
    return x, y, wy @ np.maximum(irradiance, 0) @ wx.T


def calculate_detector(channel, wavelength_nm, parameters, rng=None):
    parameters.validate()
    if not np.isfinite(wavelength_nm) or wavelength_nm <= 0:
        raise ValueError("Detector wavelength must be positive and finite.")
    x, y, power = integrate_pixels(channel["image_axis_um"],
                                   channel["irradiance_w_um2"], parameters)
    exposure = parameters.exposure_ms * 1e-3
    photon_energy = 6.62607015e-34 * 299792458 / (wavelength_nm * 1e-9)
    photons = power * exposure / photon_energy
    signal = photons * parameters.quantum_efficiency
    background = exposure * (parameters.dark_current_e_s + parameters.background_e_s)
    total = signal + background
    if rng is None:
        rng = np.random.default_rng(parameters.seed)
    charge = rng.poisson(total)
    saturated = charge >= parameters.full_well_e
    measured = np.minimum(charge, parameters.full_well_e).astype(float)
    measured += rng.normal(0, parameters.read_noise_e, measured.shape)
    variance = total + parameters.read_noise_e**2
    snr = np.divide(signal, np.sqrt(variance), out=np.zeros_like(signal), where=variance > 0)
    variance_sum = float(variance.sum())
    roi_power = float(power.sum())
    return dict(pixel_x_um=x, pixel_y_um=y, pixel_power_w=power,
                mean_photons=photons, mean_signal_e=signal, mean_total_e=total,
                measured_e=measured, saturated=saturated, ideal_snr=snr,
                roi_power_w=roi_power,
                roi_power_fraction=roi_power/channel["power_w"] if channel["power_w"] > 0 else 0.0,
                signal_e=float(signal.sum()), background_e=float(background*signal.size),
                ideal_roi_snr=float(signal.sum()/np.sqrt(variance_sum)) if variance_sum > 0 else 0.0,
                saturated_pixels=int(saturated.sum()),
                expected_overfull_pixels=int(np.count_nonzero(total >= parameters.full_well_e)))

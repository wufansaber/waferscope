"""Backend-independent area-integrating detector ROIs for sampled irradiance."""

from dataclasses import dataclass

import numpy as np
from camera_readout import digitize, sensor_maps
from emva1288_noise import noise_metrics


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
    gain_e_dn: float = 0.5
    offset_dn: float = 100.0
    bit_depth: int = 16
    fill_factor: float = 1.0
    prnu_fraction: float = 0.0
    dsnu_e: float = 0.0
    dark_current_std_e_s: float = 0.0
    fixed_pattern_seed: int = 0
    reference_frames: int = 1
    haze_mismatch_fraction: float = 0.0

    def validate(self):
        positive = (self.pixel_pitch_um, self.exposure_ms, self.full_well_e, self.gain_e_dn, self.fill_factor)
        nonnegative = (self.read_noise_e, self.dark_current_e_s, self.background_e_s,
                       self.offset_dn, self.prnu_fraction, self.dsnu_e, self.dark_current_std_e_s,
                       self.haze_mismatch_fraction)
        if not np.isfinite(positive + nonnegative + (self.quantum_efficiency,)).all():
            raise ValueError("Detector parameters must be finite.")
        if min(positive) <= 0 or min(nonnegative) < 0:
            raise ValueError("Pixel pitch, exposure and full well must be positive; noise rates nonnegative.")
        if not 0 <= self.quantum_efficiency <= 1:
            raise ValueError("Quantum efficiency must be between 0 and 1.")
        if self.fill_factor > 1 or self.prnu_fraction > .5:
            raise ValueError('Fill factor must be in (0,1]; PRNU fraction in [0,0.5].')
        if isinstance(self.bit_depth, bool) or not isinstance(self.bit_depth, int) or not 8 <= self.bit_depth <= 16:
            raise ValueError('ADC bit depth must be an integer from 8 to 16')
        if self.offset_dn > 2**self.bit_depth-1:
            raise ValueError('Black level exceeds the ADC range')
        if (isinstance(self.reference_frames, bool) or not isinstance(self.reference_frames, int)
                or not 1 <= self.reference_frames <= 128):
            raise ValueError('Reference frames must be an integer from 1 to 128')
        if self.haze_mismatch_fraction > 10:
            raise ValueError('Residual haze mismatch fraction must be in [0,10]')
        for count in (self.pixels_x, self.pixels_y):
            if not isinstance(count, int) or isinstance(count, bool) or not 1 <= count <= 512:
                raise ValueError("Detector ROI dimensions must be integers from 1 to 512.")
        for seed in (self.seed, self.fixed_pattern_seed):
            if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2**32 - 1:
                raise ValueError("Noise seeds must be integers from 0 to 4294967295.")


def _overlap_weights(axis, count, pitch, fill_factor=1):
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
    centers = (edges[:-1] + edges[1:])/2
    half_active = pitch*np.sqrt(fill_factor)/2
    overlap = np.maximum(0, np.minimum(centers[:, None]+half_active, axis[None, :] + step/2)
                         - np.maximum(centers[:, None]-half_active, axis[None, :] - step/2))
    return centers, overlap


def integrate_pixels(axis, irradiance, parameters):
    """Integrate a piecewise-constant image grid, without wrapping or extrapolation."""
    parameters.validate()
    irradiance = np.asarray(irradiance, dtype=float)
    if irradiance.shape != (len(axis), len(axis)) or not np.isfinite(irradiance).all():
        raise ValueError("Image irradiance must be a finite square array matching the image axis.")
    peak = np.max(np.abs(irradiance))
    if np.min(irradiance) < -1e-10 * peak:
        raise ValueError("Image has negative axial flux; this detector model requires nonnegative irradiance.")
    x, wx = _overlap_weights(axis, parameters.pixels_x, parameters.pixel_pitch_um, parameters.fill_factor)
    y, wy = _overlap_weights(axis, parameters.pixels_y, parameters.pixel_pitch_um, parameters.fill_factor)
    return x, y, wy @ np.maximum(irradiance, 0) @ wx.T


def calculate_detector(channel, wavelength_nm, parameters, rng=None, sensor_index=0):
    parameters.validate()
    x, y, power = integrate_pixels(channel["image_axis_um"],
                                   channel["irradiance_w_um2"], parameters)
    result = calculate_pixel_detector(x, y, power, wavelength_nm, parameters, rng, sensor_index)
    result['roi_power_fraction'] = result['roi_power_w']/channel['power_w'] if channel['power_w'] > 0 else 0.
    return result


def calculate_pixel_detector(x, y, power, wavelength_nm, parameters, rng=None, sensor_index=0):
    """Shared readout for an already area-integrated, nonnegative pixel power map."""
    parameters.validate()
    if not np.isfinite(wavelength_nm) or wavelength_nm <= 0:
        raise ValueError("Detector wavelength must be positive and finite.")
    x, y, power = np.asarray(x, float), np.asarray(y, float), np.asarray(power, float)
    if (x.shape != (parameters.pixels_x,) or y.shape != (parameters.pixels_y,)
            or power.shape != (len(y), len(x)) or not np.isfinite(power).all()
            or not np.isfinite(x).all() or not np.isfinite(y).all() or np.any(power < 0)):
        raise ValueError('Finite pixel coordinates and nonnegative integrated pixel powers required')
    exposure = parameters.exposure_ms * 1e-3
    photon_energy = 6.62607015e-34 * 299792458 / (wavelength_nm * 1e-9)
    photons = power * exposure / photon_energy
    qe_map, dark_map, dark_offset = sensor_maps(parameters, photons.shape, sensor_index)
    signal = photons * qe_map
    background = exposure * (dark_map + parameters.background_e_s)
    total = signal + background
    if rng is None:
        rng = np.random.default_rng(parameters.seed)
    charge = rng.poisson(total)
    saturated = charge >= parameters.full_well_e
    measured = np.minimum(charge, parameters.full_well_e).astype(float)
    measured += rng.normal(0, parameters.read_noise_e, measured.shape)
    measured += dark_offset
    variance = total + parameters.read_noise_e**2
    snr = np.divide(signal, np.sqrt(variance), out=np.zeros_like(signal), where=variance > 0)
    variance_sum = float(variance.sum())
    roi_power = float(power.sum())
    result = dict(pixel_x_um=x, pixel_y_um=y, pixel_power_w=power,
                qe_map=qe_map, dark_current_map_e_s=dark_map, dsnu_offset_e=dark_offset, mean_background_e=background,
                mean_photons=photons, mean_signal_e=signal, mean_total_e=total,
                measured_e=measured, saturated=saturated, ideal_snr=snr,
                roi_power_w=roi_power,
                signal_e=float(signal.sum()), background_e=float(background.sum()),
                ideal_roi_snr=float(signal.sum()/np.sqrt(variance_sum)) if variance_sum > 0 else 0.0,
                saturated_pixels=int(saturated.sum()),
                expected_overfull_pixels=int(np.count_nonzero(total >= parameters.full_well_e)))
    result.update(digitize(measured, total, parameters, fixed_offset_e=dark_offset))
    result['nominal_signal_e'] = photons*parameters.quantum_efficiency
    result.update(noise_metrics(signal, background, parameters, result['nominal_signal_e']))
    result['digital_snr_valid'] &= ((np.sqrt(variance)/parameters.gain_e_dn >= .4)
        & ~saturated & ~result['adc_clipped_low'] & ~result['adc_saturated'])
    result['snr_valid_pixels'] = int(result['digital_snr_valid'].sum())
    result['sensor_index'] = sensor_index
    return result

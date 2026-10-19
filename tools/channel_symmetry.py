"""Numerical and statistical helpers for left/right optical symmetry audits."""
import numpy as np
from scipy.stats import t

from detector_model import integrate_pixels
from fullwave_prototype.imaging import image_from_amplitude


REFLECTION = np.array([-1., 1., 1.])


def reflect_centered_grid(values):
    """Reflect X on an even FFT grid whose origin is at index n/2."""
    return np.roll(values[:, ::-1], 1, axis=1)


def reflect_sources(data):
    result = dict(data)
    for name in ('with', 'without'):
        for prefix in ('source_points_', 'source_polarization_weights_'):
            result[prefix+name] = data[prefix+name]*REFLECTION
    result['local_spm_positions_um'] = data['local_spm_positions_um']*REFLECTION[:2]
    result['local_spm_weights'] = data['local_spm_weights']*REFLECTION
    result['whole_weighted_fields'] = reflect_centered_grid(data['whole_weighted_fields'])*REFLECTION
    return result


def relative_error(a, b):
    return float(np.linalg.norm(np.asarray(a)-b)/max(np.linalg.norm(b), 1e-300))


def paired_summary(values, confidence=.95, margin=.01):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 2 or values.shape[1] != 2 or len(values) < 2
            or not np.isfinite(values).all() or np.any(values < 0)):
        raise ValueError('At least two finite nonnegative paired powers are required')
    difference = values[:, 0]-values[:, 1]
    scale = float(values.mean())
    if scale <= 0 or not 0 < confidence < 1 or not 0 < margin < 1:
        raise ValueError('Positive power, confidence and equivalence margin are required')
    mean = float(difference.mean())
    half = float(t.ppf((1+confidence)/2, len(values)-1)*difference.std(ddof=1)/np.sqrt(len(values)))
    interval = [mean-half, mean+half]
    return dict(samples=len(values), means_w=values.mean(axis=0).tolist(),
        confidence=confidence, paired_difference_mean_w=mean,
        paired_difference_ci_w=interval,
        mean_difference_relative_to_pooled_mean=mean/scale,
        normalized_difference_ci=[v/scale for v in interval],
        zero_in_ci=bool(interval[0] <= 0 <= interval[1]),
        equivalence_margin=margin,
        ci_within_equivalence_margin=bool(max(abs(v/scale) for v in interval) <= margin),
        scope='Paired Student-t CI on independent surfaces; normalized by observed pooled mean. '
              'Zero in CI is not proof of equality. No mirrored pairs or optional stopping.')


def channel_image(amplitude, geometry, p, settings, camera):
    local, q, rotation, mask, du = geometry
    field = np.zeros(mask.shape+(2,), complex)
    field[mask] = amplitude
    image = image_from_amplitude(field, local, q, rotation, du, p.wavelength_nm/1000, settings)
    flux = image['flux']*settings.irradiance_w_um2
    _, _, pixels = integrate_pixels(image['axis_um'], flux, camera)
    metrics = dict(aperture_power_w=image['transmitted_dsc_um2']*settings.irradiance_w_um2,
                   roi_power_w=float(pixels.sum()))
    return metrics, dict(axis_um=image['axis_um'], flux=flux, pixels=pixels,
                         pupil_intensity=np.sum(abs(field)**2, axis=-1))

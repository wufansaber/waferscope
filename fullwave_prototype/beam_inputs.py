"""Finite coherent wafer envelopes and their exact pixel/rectangle Fourier transforms."""
from functools import lru_cache
import hashlib
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre

DEFINITIONS = ('wafer_gaussian', 'beam_cross_section', 'wafer_rectangle', 'measured_wafer')


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=8)
def measured_envelope(path, digest, phase_assumption):
    if fingerprint(path) != digest:
        raise ValueError('Measured illumination file changed; configure and solve again')
    with np.load(path, allow_pickle=False) as saved:
        x, y, intensity = (np.asarray(saved[key], float) for key in ('x_um', 'y_um', 'intensity'))
        phase = np.asarray(saved['phase_envelope_rad'], float) if 'phase_envelope_rad' in saved else None
    supplied_phase = phase is not None
    for axis in (x, y):
        if (axis.ndim != 1 or not 4 <= len(axis) <= 256 or not np.isfinite(axis).all() or
                np.any(np.diff(axis) <= 0) or not np.allclose(np.diff(axis), axis[1]-axis[0])):
            raise ValueError('Measured beam needs uniform increasing X/Y axes, 4..256 samples each')
    if (intensity.shape != (len(y), len(x)) or not np.isfinite(intensity).all() or
            np.any(intensity < 0) or intensity.max() <= 0):
        raise ValueError('Measured intensity must be finite, nonnegative and nonzero, indexed [Y,X]')
    if phase is None:
        if phase_assumption != 'flat':
            raise ValueError('Intensity does not determine phase; supply phase_envelope_rad or choose flat explicitly')
        phase = np.zeros_like(intensity)
    elif phase.shape != intensity.shape or not np.isfinite(phase).all():
        raise ValueError('Measured envelope phase must be finite and match intensity')
    amplitude = np.sqrt(intensity/intensity.max())*np.exp(1j*phase)
    return x, y, amplitude, dict(sha256=digest, phase='supplied envelope phase' if
        supplied_phase else 'assumed flat envelope phase (intensity only)',
        phase_assumption=phase_assumption, input_extent_x_um=len(x)*(x[1]-x[0]),
        input_extent_y_um=len(y)*(y[1]-y[0]), origin='Explicit wafer coordinates; carrier +kY is added separately',
        warnings=['Field is zero outside measured pixel support. Unknown tails are not reconstructed.',
                  'Fully coherent envelope and fixed transported s/p polarization are assumed.'])


def wafer_dimensions(p):
    if p.beam_definition == 'beam_cross_section':
        return p.beam_length_um/np.cos(np.deg2rad(p.incidence_deg)), p.beam_width_um
    if p.beam_definition == 'measured_wafer':
        _, _, _, info = measured_envelope(p.beam_measurement_path, p.beam_measurement_sha256,
                                         p.beam_phase_assumption)
        return info['input_extent_y_um'], info['input_extent_x_um']
    return p.beam_length_um, p.beam_width_um


def disk_quadrature(order):
    nodes, weights = roots_legendre(order)
    a = np.pi*nodes/2
    qx = np.sin(a)
    edge = np.cos(a)
    qy = edge[:, None]*np.sin(a)
    qz = -edge[:, None]*np.cos(a)
    w = (weights*edge*np.pi/2)[:, None]*(weights*np.cos(a)*np.pi/2)*edge[:, None]
    directions = np.stack((np.broadcast_to(qx[:, None], qy.shape), qy, qz), axis=-1)
    return directions, w


def envelope_fourier(p, q):
    kx, ky = p.k0*q[..., 0], p.k0*(q[..., 1]-np.sin(np.deg2rad(p.incidence_deg)))
    if p.beam_definition == 'wafer_rectangle':
        return p.beam_length_um*p.beam_width_um*np.sinc(kx*p.beam_width_um/(2*np.pi))*np.sinc(
            ky*p.beam_length_um/(2*np.pi))
    x, y, amplitude, _ = measured_envelope(p.beam_measurement_path, p.beam_measurement_sha256,
                                          p.beam_phase_assumption)
    dx, dy = x[1]-x[0], y[1]-y[0]
    result = np.empty(kx.shape, complex)
    for i in range(len(q)):
        ex = np.exp(-1j*kx[i, 0]*x)
        ey = np.exp(-1j*ky[i, :, None]*y)
        result[i] = (ey@(amplitude@ex))*dx*dy*np.sinc(kx[i]*dx/(2*np.pi))*np.sinc(ky[i]*dy/(2*np.pi))
    return result


def spectrum(p, order=None):
    if p.beam_definition == 'measured_wafer' and fingerprint(p.beam_measurement_path) != p.beam_measurement_sha256:
        raise ValueError('Measured illumination file changed; configure and solve again')
    if p.beam_definition in ('wafer_gaussian', 'beam_cross_section'):
        from dataclasses import replace
        from .finite_illumination import gaussian_spectrum
        length, width = wafer_dimensions(p)
        return gaussian_spectrum(replace(p, beam_length_um=length, beam_width_um=width), order)
    q, w = disk_quadrature(order or p.beam_order)
    amplitudes = envelope_fourier(p, q)*w*(p.k0/(2*np.pi))**2
    return q.reshape(-1, 3), amplitudes.ravel()


def raw_area(p, order=128, flux=True):
    q, w = disk_quadrature(order)
    values = abs(envelope_fourier(p, q))**2*w*(p.k0/(2*np.pi))**2
    return float(np.sum(values*(-q[..., 2] if flux else 1)))


def truncation_report(p):
    if p.beam_definition in ('wafer_gaussian', 'beam_cross_section'):
        # Scalar-envelope Parseval power before the propagating-disk cutoff.
        length, width = wafer_dimensions(p)
        total = np.pi*length*width/8
        from dataclasses import replace
        from .finite_illumination import _raw_power_area
        effective = replace(p, beam_length_um=length, beam_width_um=width)
        # Beam-scaled integration resolves narrow spectra of large illuminated fields.
        retained = _raw_power_area(effective, 192, flux=False)/total
        check = _raw_power_area(effective, 128, flux=False)/total
        if not np.isfinite(retained) or abs(retained-check) > 1e-5:
            raise ValueError('Gaussian scalar spectral coverage integral did not converge')
        metadata = {}
    else:
        if p.beam_definition == 'measured_wafer':
            x, y, amplitude, metadata = measured_envelope(p.beam_measurement_path,
                p.beam_measurement_sha256, p.beam_phase_assumption)
            total = np.sum(abs(amplitude)**2)*(x[1]-x[0])*(y[1]-y[0])
        else:
            total, metadata = p.beam_length_um*p.beam_width_um, {}
        retained = raw_area(p, 192, flux=False)/total
    return dict(definition=p.beam_definition, declared_length_y_um=p.beam_length_um,
        declared_width_x_um=p.beam_width_um, wafer_length_y_um=wafer_dimensions(p)[0],
        wafer_width_x_um=wafer_dimensions(p)[1], retained_scalar_spectral_fraction=retained,
        excluded_scalar_spectral_fraction=max(0., 1-retained), measurement=metadata,
        scope='Scalar envelope Parseval spectral weight, NOT transmitted/scattered energy. '
              'Evanescent components are omitted; vector reconstruction may differ from the target intensity.')

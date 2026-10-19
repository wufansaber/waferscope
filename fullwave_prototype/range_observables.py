"""Solid-angle-weighted complex field differences, with unweighted probes separate."""
import numpy as np

from .aperture_integrals import cap_quadrature


def sample_caps(radiation, apertures, tilts, order=64):
    return [(radiation.amplitude(q), w) for q, w in
            (cap_quadrature(float(na), float(tilt), order) for na, tilt in zip(apertures, tilts))]


def cap_powers(samples):
    return np.array([float(w@np.sum(abs(a)**2, axis=-1)) for a, w in samples])


def cap_field_changes(current, reference):
    return [float(np.sqrt((w@np.sum(abs(a-b)**2, axis=-1))/
                         max(float(w@np.sum(abs(b)**2, axis=-1)), 1e-300)))
            for (a, w), (b, _) in zip(current, reference)]

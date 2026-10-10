"""Independent high-order surface integration preserves the original vector field."""
from dataclasses import replace

import numpy as np
import pytest

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.planar import background_vectors
from fullwave_prototype.illumination import matched_ray
from fullwave_prototype.finite_illumination import gaussian_spectrum
from fullwave_prototype.surface_scattering import surface_fields
from fullwave_prototype.surface_beam import plane_ray_vectors, surface_fields_order, checked_surface_order


@pytest.mark.parametrize('polarization', ['s', 'p'])
@pytest.mark.parametrize('incidence', [0, 75])
def test_vectorized_matching_and_surface_phase_agree(polarization, incidence):
    p = Parameters(beam_profile='gaussian', wavelength_nm=266, beam_order=21,
                   beam_length_um=2, beam_width_um=2,
                   polarization=polarization, incidence_deg=incidence)
    q, _ = gaussian_spectrum(p)
    q = np.vstack((q[::17], [0, 0, -1]))
    incoming, reflected, magnetic = plane_ray_vectors(p, q)
    down0, _, _, e0, _, _ = background_vectors(p)
    rays = [matched_ray(p, down0, e0, down, 1) for down in q]
    np.testing.assert_allclose(incoming, [r[3] for r in rays], atol=2e-14)
    np.testing.assert_allclose(reflected, [r[4] for r in rays], atol=2e-14)
    np.testing.assert_allclose(magnetic, np.cross(q, incoming), atol=2e-14)
    x, y = np.linspace(-2, 2, 7), np.linspace(-3, 3, 9)
    for actual, expected in zip(surface_fields_order(x, y, p, 21), surface_fields(x, y, p)):
        np.testing.assert_allclose(actual, expected, atol=2e-13)


def test_high_order_hermite_is_finite_and_does_not_change_physical_parameters():
    p = Parameters(beam_profile='gaussian', wavelength_nm=266, beam_order=41,
                   beam_length_um=2, beam_width_um=2)
    q, weights = gaussian_spectrum(p, 401)
    assert np.isfinite(q).all() and np.isfinite(weights).all()
    np.testing.assert_allclose(np.linalg.norm(q, axis=1), 1, atol=2e-14)
    for field in surface_fields_order(np.array([0., 1.]), np.array([0., 1.]), p, 401):
        assert np.isfinite(field).all()
    assert p.beam_order == 41
    with pytest.raises(ValueError, match='Gaussian beam'):
        surface_fields_order([0], [0], replace(p, beam_profile='plane'), 21)


def test_order_acceptance_requires_two_stable_steps(monkeypatch):
    p = Parameters(beam_profile='gaussian', beam_order=41)
    used = []

    def fields(x, y, parameters, order):
        used.append(order)
        return (np.ones((1, 1, 3)),)*3

    monkeypatch.setattr('fullwave_prototype.surface_beam.surface_fields_order', fields)
    order, checks, met = checked_surface_order(p, 28)
    assert met and order == 81 and used == [41, 61, 81] and len(checks) == 2

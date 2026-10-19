"""Footprint definitions, measured phase, spectral cutoff and normalization."""
from dataclasses import replace
import numpy as np
import pytest
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.beam_inputs import (wafer_dimensions, spectrum, disk_quadrature, raw_area,
                                          measured_envelope, fingerprint, truncation_report)
from fullwave_prototype.finite_illumination import power_area
from fullwave_prototype.surface_beam import surface_fields_order
from fullwave_prototype.surface_scattering import surface_fields
from fullwave_prototype.background_sampling import sampled_background
from fullwave_prototype.planar import background_fields


def beam(**values):
    return replace(Parameters(beam_profile='gaussian', wavelength_nm=266, beam_length_um=2,
                              beam_width_um=2, beam_order=61), **values)


def test_cross_section_projection_changes_physical_beam_and_power():
    p = beam(beam_definition='beam_cross_section', incidence_deg=60)
    assert wafer_dimensions(p) == pytest.approx((4, 2))
    same = beam(beam_length_um=4, incidence_deg=60)
    np.testing.assert_allclose(spectrum(p)[0], spectrum(same)[0])
    assert power_area(p) == pytest.approx(power_area(same), rel=1e-12)
    assert power_area(p) > power_area(beam(incidence_deg=60))


def test_disk_jacobian_and_rectangle_integrated_power():
    q, w = disk_quadrature(64)
    assert w.sum() == pytest.approx(np.pi, rel=1e-12)
    np.testing.assert_allclose(np.linalg.norm(q, axis=-1), 1, atol=2e-14)
    p = beam(beam_definition='wafer_rectangle', incidence_deg=0)
    assert 0 < raw_area(p)/4 < 1
    assert raw_area(p, 192) == pytest.approx(raw_area(p, 128), rel=1e-5)
    axis = np.array([-.1, 0, .2])
    for actual, expected in zip(surface_fields_order(axis, axis, p, 61), surface_fields(axis, axis, p)):
        np.testing.assert_allclose(actual, expected, atol=2e-13)
    info = truncation_report(p)
    assert 0 < info['excluded_scalar_spectral_fraction'] < 1


def test_measured_phase_is_explicit_and_fingerprint_is_enforced(tmp_path):
    axis = (np.arange(24)-12)*.15
    x, y = np.meshgrid(axis, axis)
    intensity = np.exp(-2*(x*x+y*y))
    path = tmp_path/'beam.npz'
    np.savez(path, x_um=axis, y_um=axis, intensity=intensity)
    sha = fingerprint(path)
    with pytest.raises(ValueError, match='does not determine phase'):
        measured_envelope(str(path), sha, 'reject')
    _, _, amplitude, report = measured_envelope(str(path), sha, 'flat')
    assert 'assumed' in report['phase']
    np.testing.assert_allclose(abs(amplitude)**2, intensity)
    p = beam(beam_definition='measured_wafer', beam_measurement_path=str(path),
             beam_measurement_sha256=sha, incidence_deg=0)
    p.validate()
    assert np.isfinite(power_area(p)) and power_area(p) > 0
    phi = .2*x
    np.savez(path, x_um=axis, y_um=axis, intensity=intensity, phase_envelope_rad=phi)
    with pytest.raises(ValueError, match='unchanged'):
        p.validate()
    with pytest.raises(ValueError, match='changed'):
        spectrum(p)
    _, _, field, report = measured_envelope(str(path), fingerprint(path), 'reject')
    np.testing.assert_allclose(field, np.sqrt(intensity)*np.exp(1j*phi))
    assert report['phase'] == 'supplied envelope phase'


@pytest.mark.parametrize('polarization', ['s', 'p'])
def test_batch_magnetic_background_matches_native_ray_reference(polarization):
    p = beam(polarization=polarization, beam_order=21)
    xyz = np.array([[.1, .2, 0], [0, 0, -.05], [.2, -.1, .1]])
    actual = sampled_background(xyz, p, magnetic=True, block=2)
    expected = background_fields(xyz, p)
    for a, e in zip(actual, expected):
        np.testing.assert_allclose(a, e, atol=2e-13)


def test_local_general_envelope_rejects_unresolved_quadrature(monkeypatch):
    from fullwave_prototype.beam_validation import check_local_beam
    monkeypatch.setattr('fullwave_prototype.beam_validation.electric_background',
                        lambda points, p: np.ones_like(points)*p.beam_order)
    with pytest.raises(ValueError, match='increase beam order'):
        check_local_beam(beam(beam_definition='wafer_rectangle', beam_order=31))
    assert check_local_beam(beam()) is None


def test_measured_whole_planner_uses_independent_high_order(monkeypatch):
    from fullwave_prototype.surface_planning import coverage_at_width
    monkeypatch.setattr('fullwave_prototype.surface_beam.surface_fields_order',
        lambda x, y, p, order: tuple(np.ones((len(y), len(x), 3), complex) for _ in range(3)))
    monkeypatch.setattr('fullwave_prototype.surface_planning.power_area', lambda p: 1)
    candidate, plan = coverage_at_width(beam(beam_definition='wafer_rectangle'), 2, .03, .05)
    assert candidate.beam_order == 61 and plan['surface_beam_order'] == 61


def test_blocked_surface_variant_preserves_complex_fields_and_power(monkeypatch):
    from fullwave_prototype.hybrid_surface_audit import surface_variant
    def fields(x, y, p):
        xx, yy = np.meshgrid(x, y)
        incident = np.zeros(xx.shape+(3,), complex)
        incident[..., 0] = np.exp(1j*(xx+2*yy))
        magnetic = np.cross(np.array([0, 0, -1]), incident)
        return 1.2*incident, incident, magnetic
    monkeypatch.setattr('fullwave_prototype.hybrid_surface_audit.surface_fields', fields)
    n, spacing = 300, .01
    axis = (np.arange(n)-n//2)*spacing
    x, y = np.meshgrid(axis, axis)
    height = .002*np.cos(x)*np.sin(y)
    p = beam()
    reference = surface_variant(axis, height, n, spacing, .05, p, 2)
    actual = surface_variant(axis, height, n, spacing, .05, p, 2, fields_provider=fields)
    assert actual[1] == pytest.approx(reference[1], rel=1e-13)
    np.testing.assert_array_equal(actual[3], reference[3])
    directions = np.array([[0, 0, 1], [.1, .2, np.sqrt(.95)]])
    np.testing.assert_allclose(actual[0].amplitude(directions), reference[0].amplitude(directions), atol=1e-15)


def test_measured_power_coverage_does_not_alias_pixel_beat_fringes(monkeypatch):
    from fullwave_prototype.surface_planning import coverage_at_width
    def fringes(x, y, p, order):
        xx, yy = np.meshgrid(x, y)
        e = np.zeros(xx.shape+(3,), complex)
        e[..., 0] = np.sqrt(1+.8*np.cos(2*np.pi*yy/(2/96)))
        h = np.cross(np.array([0, 0, -1]), e)
        return e, e, h
    monkeypatch.setattr('fullwave_prototype.surface_beam.surface_fields_order', fringes)
    monkeypatch.setattr('fullwave_prototype.whole_surface.window',
                        lambda x, y, width, taper: np.ones((len(y), len(x))))
    monkeypatch.setattr('fullwave_prototype.surface_planning.power_area', lambda p: 4)
    _, plan = coverage_at_width(beam(beam_definition='wafer_rectangle'), 2, .02, .05)
    assert plan['power_coverage_samples_per_axis'] == 128
    assert plan['planning_power_coverage'] == pytest.approx(1., abs=1e-13)

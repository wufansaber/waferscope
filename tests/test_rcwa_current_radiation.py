"""Solid-angle quadrature and provenance/fail-closed diagnostic checks."""
import hashlib
import json
from dataclasses import asdict, replace

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from rcwa_prototype.current_radiation import cone_quadrature, integrate_cones, lattice_diagnostic, load_source, matching_scenes, run


def test_scene_matching_allows_serialization_roundoff_but_not_physical_changes():
    scene = Scene(lx_um=2/3)
    assert matching_scenes(scene, replace(scene, lx_um=.666666666666667))
    assert not matching_scenes(scene, replace(scene, lx_um=.667))
    assert not matching_scenes(scene, replace(scene, substrate_k=0.))
    assert not matching_scenes(scene, replace(scene, slices=scene.slices+1))


@pytest.mark.parametrize('tilt', [45., -45., 0.])
def test_cone_has_exact_solid_angle_and_correct_first_moment(tilt):
    na = .6
    rays, weights = cone_quadrature(na, tilt, 12, 36)
    axis = np.array([np.sin(np.deg2rad(tilt)), 0., np.cos(np.deg2rad(tilt))])
    np.testing.assert_allclose(np.linalg.norm(rays, axis=1), 1., atol=1e-14)
    assert np.min(rays@axis) > np.sqrt(1-na*na)
    assert weights.sum() == pytest.approx(2*np.pi*(1-np.sqrt(1-na*na)), rel=1e-14)
    np.testing.assert_allclose(weights@rays, np.pi*na*na*axis, atol=1e-14)


def test_absolute_constant_dsc_integrates_without_normalizing(monkeypatch):
    import rcwa_prototype.current_radiation as radiation
    monkeypatch.setattr(radiation, 'project_ray', lambda *a, **kw: np.array([1e-8, 2e-8j, 0.]))
    result = integrate_cones({}, Scene(), .6, (45., -45., 0.), 8, 24, 20)
    np.testing.assert_allclose(result['collected_cross_section_um2'], 5e-16*2*np.pi*.2, rtol=1e-14)


def test_disabled_cone_collects_zero():
    _, weights = cone_quadrature(0., 0., 2, 4)
    np.testing.assert_array_equal(weights, 0.)


def test_lattice_weights_match_modal_cross_section_normalization(monkeypatch):
    import rcwa_prototype.current_radiation as radiation
    from rcwa_prototype.comparison import modal_power
    scene = Scene()
    rays = np.array([[.3, 0., np.sqrt(.91)], [0., 0., 1.]])
    far = np.array([1e-4, 2e-4j, 0.])
    monkeypatch.setattr(radiation, 'project_ray', lambda *a, **kw: far)
    grid = dict(directions=rays, propagating=np.ones(2, bool), orders=np.array([[0, 0], [1, 0]]))
    result = lattice_diagnostic({}, scene, grid, scene, tilts=(0.,))
    area = scene.lx_um*scene.ly_um
    k0 = 2*np.pi/scene.wavelength_um
    modal = 2j*np.pi/(area*k0*rays[:, 2, None])*far
    expected = area*np.cos(np.deg2rad(scene.incidence_deg))*modal_power(modal, rays, grid['propagating'], scene.incidence_deg).sum()
    assert result['collected_cross_section_um2'][0] == pytest.approx(expected, rel=1e-14)


@pytest.mark.parametrize('na,tilt', [(1., 0.), (.6, 60.), (-.1, 0.), (np.nan, 0.)])
def test_invalid_or_horizon_crossing_cone_is_rejected(na, tilt):
    with pytest.raises(ValueError):
        cone_quadrature(na, tilt, 8, 24)


def test_source_raw_hash_is_verified(tmp_path):
    raw = tmp_path/'modes.npz'
    np.savez_compressed(raw, electric=np.ones((2, 3), complex))
    info = dict(scene=asdict(Scene()), raw_npz_sha256=hashlib.sha256(raw.read_bytes()).hexdigest())
    (tmp_path/'report.json').write_text(json.dumps(info), encoding='utf-8')
    data, _, _, provenance = load_source(tmp_path)
    assert data['electric'].dtype == complex
    assert provenance['sha256'] == info['raw_npz_sha256']
    np.savez_compressed(raw, electric=np.zeros((2, 3)))
    with pytest.raises(ValueError, match='hash'):
        load_source(tmp_path)


def test_duplicate_levels_rejected_before_any_computation(tmp_path):
    with pytest.raises(ValueError, match='Distinct integer'):
        run('missing', 'missing', tmp_path/'new', [(8, 24, 20), (8, 24, 20)])
    assert not (tmp_path/'new').exists()


def test_failed_modal_gate_emits_no_continuous_results(tmp_path, monkeypatch):
    import rcwa_prototype.current_radiation as radiation
    scene = Scene(grating_method='analytic_li')
    data = dict(source_electric_fourier=np.zeros((1, 1, 3)),
                electric_reflected_z0=np.ones((1, 3)), directions=np.array([[0., 0., 1.]]),
                orders=np.array([[0, 0]]))
    entries = iter([(data, scene, dict(particle=True, current_z_order=8), {}),
                    (data, scene, dict(particle=False), {})])
    monkeypatch.setattr(radiation, 'load_source', lambda _: next(entries))
    gate = dict(comparable=True, relative_power_weighted_complex_field_change=.01,
                relative_total_power_difference=.02)
    monkeypatch.setattr(radiation, 'modal_reconstruction', lambda *a: (np.zeros((1, 3)), gate))
    with pytest.raises(RuntimeError, match='gate failed'):
        run('source', 'reference', tmp_path/'result', [(8, 24, 20)])
    assert (tmp_path/'result'/'gate.json').exists()
    assert not list((tmp_path/'result').glob('cone*.npz'))

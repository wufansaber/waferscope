"""Analytic grating convolution and grcwa adapter regression tests."""
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from rcwa_prototype.grating_fourier import line_convolution, material_matrices


@pytest.mark.parametrize('periods', [1, 2, 3])
def test_exact_mean_spectrum_and_forbidden_couplings(periods):
    scene = Scene(periods_y=periods, ridge_width_um=.2)
    orders = np.array([[0, 0], [0, periods], [0, -periods], [1, 0], [0, 1]])
    matrix = line_convolution(orders, scene)
    eta = line_convolution(orders, scene, reciprocal=True)
    fill = .4
    np.testing.assert_allclose(np.diag(matrix), 1+fill*(scene.eps_substrate-1))
    np.testing.assert_allclose(np.diag(eta), 1+fill*(1/scene.eps_substrate-1))
    assert matrix[0, 3] == 0
    if periods > 1:
        assert matrix[0, 4] == 0
    y = (np.arange(50000)+.5)*scene.ly_um/50000
    local = (y-scene.ly_um/2+scene.pitch_um/2) % scene.pitch_um-scene.pitch_um/2
    epsilon = np.where(abs(local) < scene.ridge_width_um/2, scene.eps_substrate, 1.)
    expected = np.mean(epsilon*np.exp(2j*np.pi*periods*y/scene.ly_um))
    assert matrix[0, 1] == pytest.approx(expected, abs=2e-3)


def test_li_blocks_reduce_to_uniform_material_without_ridges():
    scene = Scene(ridge_width_um=.5)
    orders = np.array([[0, 0], [0, 1], [0, -2]])
    epinv, ep2 = material_matrices(orders, scene)
    np.testing.assert_allclose(epinv, np.eye(3)/scene.eps_substrate, atol=1e-15)
    np.testing.assert_allclose(ep2, np.eye(6)*scene.eps_substrate, atol=1e-12)


def test_grating_rule_is_explicit_and_old_scene_defaults_to_fft():
    assert Scene().grating_method == 'fft'
    with pytest.raises(ValueError, match='material rule'):
        replace(Scene(), grating_method='unknown').validate()


def test_analytic_background_independent_of_unused_pixel_grid():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    scene = Scene(modes=31, lx_um=.005, periods_y=1, grating_method='analytic_li')
    a, info = solve(scene)
    b, _ = solve(replace(scene, grid_x=3, grid_y=5))
    assert info['grid_layers'] == 0
    assert info['analytic_grating_layer'] is not None
    np.testing.assert_allclose(a['electric_reflected_z0'], b['electric_reflected_z0'], atol=1e-12)


def test_li_background_against_saved_nist_order_powers():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    path = Path(__file__).resolve().parents[1]/'output/rcwa_trial_20261021/nist245.npz'
    if not path.is_file():
        pytest.skip('Saved independent NIST grating oracle not available')
    scene = Scene(modes=121, lx_um=.005, periods_y=1, grating_method='analytic_li')
    data, _ = solve(scene)
    with np.load(path, allow_pickle=False) as saved:
        for order, expected in zip(saved['orders'], saved['efficiency']):
            if expected > 1e-8:
                match = np.all(data['orders'] == [0, order], axis=1)
                assert data['reflectance_by_order'][match].item() == pytest.approx(expected, rel=.015)


def test_exact_block_eigensolve_preserves_dense_complex_fields(monkeypatch):
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    from rcwa_prototype import grating_fourier
    scene = Scene(modes=51, grating_method='analytic_li')
    a, _ = solve(scene, near=True)
    original = grating_fourier.setup_layer
    monkeypatch.setattr(grating_fourier, 'setup_layer',
        lambda model, parameters, index: original(model, parameters, index, block_diagonal=False))
    b, _ = solve(scene, near=True)
    np.testing.assert_allclose(a['electric_reflected_z0'], b['electric_reflected_z0'], atol=2e-10)
    np.testing.assert_allclose(a['near_electric'], b['near_electric'], atol=2e-10)


def test_rectangular_orders_and_memory_budget_count_actual_basis():
    from rcwa_prototype.basis import rectangular_orders
    scene = Scene(harmonics_x=2, harmonics_y=6)
    orders = rectangular_orders(scene)
    assert scene.requested_modes == 65
    assert len(set(map(tuple, orders))) == 65
    np.testing.assert_array_equal(orders[0], [0, 0])
    np.testing.assert_array_equal(np.max(abs(orders), axis=0), [2, 6])
    for change in (dict(harmonics_x=2), dict(harmonics_x=True, harmonics_y=3),
                   dict(harmonics_x=-1, harmonics_y=3)):
        with pytest.raises(ValueError):
            replace(Scene(), **change).validate()


def test_same_incident_block_in_1d_and_3d_rectangular_basis():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    scene = Scene(grating_method='analytic_li', harmonics_x=0, harmonics_y=14)
    a, _ = solve(scene)
    b, _ = solve(replace(scene, harmonics_x=2))
    for i, order in enumerate(a['orders']):
        j = np.flatnonzero(np.all(b['orders'] == order, axis=1)).item()
        np.testing.assert_allclose(a['electric_reflected_z0'][i], b['electric_reflected_z0'][j], atol=1e-11)


def test_rectangular_basis_guard_rejects_before_allocating(monkeypatch, tmp_path):
    pytest.importorskip('psutil')
    import sys
    from tools.rcwa_trial import main
    monkeypatch.setattr(sys, 'argv', ['rcwa_trial', '--output', str(tmp_path/'blocked'),
                                     '--harmonics', '1000', '1000'])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
    assert not (tmp_path/'blocked').exists()

from dataclasses import replace

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene, layers, matrix_memory_estimate


@pytest.mark.parametrize('change', [dict(modes=0), dict(slices=1.5), dict(periods_y=True),
    dict(grid_x=0), dict(lx_um=-1), dict(particle_k=-1), dict(wavelength_um=float('nan')),
    dict(ridge_width_um=.6), dict(kind='unknown')])
def test_invalid_scene(change):
    with pytest.raises(ValueError):
        replace(Scene(), **change).validate()


def test_sphere_layer_volume_and_reference_are_explicit():
    p = Scene(slices=16, grid_x=256, grid_y=512)
    stack, info = layers(p, True)
    assert abs(info['voxel_particle_volume_um3']/info['exact_particle_volume_um3']-1) < .02
    assert sum(d for d, _ in stack) == pytest.approx(p.cap_top_um+p.depth_um)
    _, reference = layers(p, False)
    assert reference['particle_present'] is False
    assert reference['voxel_particle_volume_um3'] == 0


def test_particle_cannot_overlap_periodic_copy():
    with pytest.raises(ValueError, match='copies'):
        layers(Scene(lx_um=.1), True)


def test_memory_estimate_scales_with_total_modes_not_axis_order():
    a, b = [matrix_memory_estimate(m, 9) for m in (100, 200)]
    assert b['one_complex_matrix_mib']/a['one_complex_matrix_mib'] == pytest.approx(4)


def test_flat_fresnel_phase_and_power_plus_near_field():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    from fullwave_prototype.planar import fresnel, background_fields
    from fullwave_prototype.parameters import Parameters
    p = Scene(kind='flat', modes=21, slices=2, grid_x=32, grid_y=64)
    data, info = solve(p, False, near=True)
    cosine = np.cos(np.deg2rad(p.incidence_deg))
    rp = fresnel(cosine, p.eps_substrate)[1]
    index = np.flatnonzero(np.all(data['orders'] == 0, axis=1))[0]
    expected = rp*np.array((0, -cosine, np.sin(np.deg2rad(p.incidence_deg))))
    np.testing.assert_allclose(data['electric_reflected_z0'][index], expected, atol=1e-10)
    assert info['reflectance'] == pytest.approx(abs(rp)**2, abs=1e-10)
    actual = data['near_electric'][40, 40]
    expected = background_fields(np.array([[0, 0, .02]]), Parameters(
        wavelength_nm=p.wavelength_um*1000))[0][0]
    np.testing.assert_allclose(actual, expected, atol=1e-10)


def test_identical_background_and_sample_without_contrast():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    scene = Scene(modes=21, slices=2, grid_x=32, grid_y=64, particle_n=1.)
    a, _ = solve(scene, False)
    b, _ = solve(scene, True)
    np.testing.assert_allclose(a['electric_reflected_z0'], b['electric_reflected_z0'], atol=1e-12)


def test_pair_report_keeps_interference_and_zero_na_separate():
    from rcwa_prototype.comparison import compare_pair, field_change
    scene = Scene(incidence_deg=0.)
    r = dict(orders=np.array([[0, 0]]), propagating=np.array([True]),
             directions=np.array([[0., 0., 1.]]),
             electric_reflected_z0=np.array([[1., 0., 0.]], complex),
             reflectance_by_order=np.array([1.]))
    s = dict(r, electric_reflected_z0=np.array([[2., 0., 0.]], complex),
             reflectance_by_order=np.array([4.]))
    channel = compare_pair(r, s, scene)['channels'][2]
    assert channel['actual_power_difference'] == pytest.approx(3.)
    assert channel['difference_field_efficiency'] == pytest.approx(1.)
    assert channel['interference_efficiency'] == pytest.approx(2.)
    assert compare_pair(r, s, scene, na=0)['channels'][2]['collected_orders'] == 0
    assert field_change(r, s, 0)['relative_power_weighted_complex_field_change'] == pytest.approx(.5)
    assert field_change(r, r, 0)['relative_power_weighted_complex_field_change'] == 0


def test_pair_report_rejects_mismatched_orders():
    from rcwa_prototype.comparison import compare_pair
    r = dict(orders=np.array([[0, 0]]), propagating=np.array([True]))
    with pytest.raises(ValueError, match='identical'):
        compare_pair(r, dict(r, orders=np.array([[0, 1]])), Scene())


def test_oversized_trial_rejected_before_worker_start(monkeypatch, tmp_path, capsys):
    pytest.importorskip('psutil')
    import sys
    from tools.rcwa_trial import main
    output = tmp_path/'must_not_start'
    monkeypatch.setattr(sys, 'argv', ['rcwa_trial', '--output', str(output), '--modes', '1000000'])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
    assert 'Planning allowance' in capsys.readouterr().err
    assert not output.exists()

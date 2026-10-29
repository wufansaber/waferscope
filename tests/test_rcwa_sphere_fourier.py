"""Check analytic disk geometry before using it as an implementation crosscheck."""
from dataclasses import replace

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from rcwa_prototype.sphere_fourier import disk_convolution, sphere_stack


@pytest.mark.parametrize('loss', [0., .02])
def test_mean_and_center_phase_against_independent_disk_quadrature(loss):
    scene = Scene(particle_k=loss)
    radius = .043
    orders = np.array([[0, 0], [1, 0], [0, -2], [1, 3]])
    matrix = disk_convolution(orders, scene, radius)
    fill = np.pi*radius**2/(scene.lx_um*scene.ly_um)
    np.testing.assert_allclose(np.diag(matrix), 1+(scene.eps_particle-1)*fill)
    # Integrate over polar coordinates without the Bessel transform formula.
    nodes, weights = np.polynomial.legendre.leggauss(120)
    r, w = (nodes+1)*radius/2, weights*radius/2
    phi = np.arange(600)*2*np.pi/600
    for i in range(len(orders)):
        for j in range(len(orders)):
            gx, gy = orders[i]-orders[j]
            x = scene.lx_um/2+r[:, None]*np.cos(phi)
            y = scene.ly_um/2+r[:, None]*np.sin(phi)
            integral = np.sum(w*r*np.mean(np.exp(-2j*np.pi*(gx*x/scene.lx_um+
                gy*y/scene.ly_um)), axis=1))*2*np.pi/(scene.lx_um*scene.ly_um)
            expected = (i == j)+(scene.eps_particle-1)*integral
            assert matrix[i, j] == pytest.approx(expected, abs=3e-14)


@pytest.mark.parametrize('slices', [1, 5, 8, 12])
def test_stack_volume_and_midpoint_slices_match_s4(slices):
    from s4_prototype.lua_scene import script
    scene = Scene(slices=slices)
    stack, info = sphere_stack(scene)
    _, expected = script(scene, particle=True)
    assert info['sliced_particle_volume_um3'] == pytest.approx(expected, rel=1e-14)
    assert info['sliced_particle_volume_um3']/info['exact_particle_volume_um3'] == pytest.approx(
        1+1/(2*slices**2), rel=1e-14)
    assert info['voxel_particle_volume_um3'] is None
    assert len(stack) == slices+2


def test_analytic_sphere_is_independent_of_unused_pixel_grid():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    scene = Scene(harmonics_x=1, harmonics_y=4, slices=3, grating_method='analytic_li')
    a, info = solve(scene, particle=True, particle_method='analytic')
    b, _ = solve(replace(scene, grid_x=3, grid_y=5), particle=True, particle_method='analytic')
    assert info['grid_layers'] == 0
    assert info['analytic_particle_layers'] == 3
    np.testing.assert_allclose(a['electric_reflected_z0'], b['electric_reflected_z0'], atol=1e-12)


def test_analytic_sphere_does_not_silently_switch_factorization():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    with pytest.raises(ValueError, match='Laurent'):
        solve(Scene(factorization='pol'), particle=True, particle_method='analytic')


def test_worker_refuses_to_overwrite_existing_raw_results(monkeypatch, tmp_path):
    import sys
    from tools.rcwa_trial import main
    saved = tmp_path/'modes.npz'
    saved.write_bytes(b'previous raw result')
    monkeypatch.setattr(sys, 'argv', ['rcwa_trial', '--worker', '--output', str(tmp_path)])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
    assert saved.read_bytes() == b'previous raw result'

"""Absolute prescribed-source and same-basis reciprocity checks, without fitting."""
from dataclasses import replace

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from rcwa_prototype.contrast_current import disk_transform, modal_reconstruction, project_ray
from rcwa_prototype.reciprocal_background import spectrum


def cylinder_source():
    nodes, weights = np.polynomial.legendre.leggauss(8)
    field = np.array([1.+.2j, .3, .1j])
    return dict(source_z_um=.06+.005*nodes, source_weight_um=.005*weights,
        source_radius_um=np.full(8, .03), source_layer=np.ones(8, int),
        source_electric_fourier=np.tile(field, (8, 1, 1)), source_kxy_per_um=np.zeros((1, 2)))


def test_free_space_absolute_complex_radiation_of_prescribed_cylinder():
    pytest.importorskip('grcwa')
    scene = Scene(kind='flat', substrate_n=1., substrate_k=0.)
    n = np.array([.3, -.4, np.sqrt(.75)])
    source = cylinder_source()
    result = project_ray(source, scene, n, reference_halfwidth=2)
    k0 = 2*np.pi/scene.wavelength_um
    field = source['source_electric_fourier'][0, 0]
    projected = field-n*np.dot(n, field)
    xy = disk_transform(np.array([-k0*n[0]]), np.array([-k0*n[1]]), .03).item()
    axial = .01*np.sinc(k0*n[2]*.01/(2*np.pi))*np.exp(-1j*k0*n[2]*.06)
    expected = k0*k0/(4*np.pi)*(scene.eps_particle-1)*projected*xy*axial
    np.testing.assert_allclose(result, expected, rtol=2e-11, atol=1e-13)


def test_reciprocal_planar_reference_matches_fresnel_vectors():
    pytest.importorskip('grcwa')
    scene = Scene(kind='flat')
    n = np.array([.3, -.4, np.sqrt(.75)])
    result = spectrum(scene, n, np.array([[0, 0]]))
    pol = result['polarizations']
    root = np.sqrt(scene.eps_substrate-(1-n[2]**2))
    rs = (n[2]-root)/(n[2]+root)
    rp = (scene.eps_substrate*n[2]-root)/(scene.eps_substrate*n[2]+root)
    expected_p = np.r_[-rp*pol[:2, 1], rp*pol[2, 1]]
    np.testing.assert_allclose(result['electric'][1, :, 0], rs*pol[:, 0], atol=2e-12)
    np.testing.assert_allclose(result['electric'][1, :, 1], expected_p, atol=2e-12)


@pytest.mark.parametrize('kind,particle_k', [('grating', 0.), ('flat', 0.), ('grating', .01)])
def test_contrast_current_reproduces_direct_discrete_complex_modes(kind, particle_k):
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    scene = Scene(harmonics_x=1, harmonics_y=4, slices=3, grating_method='analytic_li',
                  kind=kind, particle_k=particle_k)
    sample, _ = solve(scene, particle=True, particle_method='analytic', current_z_order=8)
    reference, _ = solve(scene, particle_method='analytic')
    _, result = modal_reconstruction(sample, scene, reference)
    assert result['relative_power_weighted_complex_field_change'] < 1e-10
    assert abs(result['relative_total_power_difference']) < 1e-10


def test_lossless_eigenvalue_branch_ignores_only_roundoff_imaginary_parts():
    from rcwa_prototype.sphere_fourier import lossless_branch
    q = np.array([-3.+1e-15j, 3.+1e-15j, 1e-15+4j, -2.+.01j])
    np.testing.assert_allclose(lossless_branch(q, 5.), [3., 3., 4j, -2.+.01j], atol=1e-13)


def test_particle_stack_tangential_fields_are_continuous(monkeypatch):
    pytest.importorskip('grcwa')
    import rcwa_prototype.contrast_current as current
    from rcwa_prototype.solver import solve
    captured = {}
    original = current.sample_current

    def capture(model, *args):
        captured['model'] = model
        return original(model, *args)

    monkeypatch.setattr(current, 'sample_current', capture)
    scene = Scene(harmonics_x=1, harmonics_y=4, slices=3, grating_method='analytic_li')
    solve(scene, particle=True, particle_method='analytic', current_z_order=8)
    model = captured['model']
    for layer in range(model.Layer_N-1):
        left = np.asarray(model.Solve_FieldFourier(layer, float(model.thickness_list[layer])))[0]
        right = np.asarray(model.Solve_FieldFourier(layer+1, 0.))[0]
        np.testing.assert_allclose(left[:, :2], right[:, :2], rtol=1e-10, atol=1e-12)


def test_current_export_rejects_a_mismatched_grating_background():
    pytest.importorskip('grcwa')
    from rcwa_prototype.solver import solve
    with pytest.raises(ValueError, match='matching analytic Li background'):
        solve(Scene(grating_method='fft'), particle=True, particle_method='analytic', current_z_order=8)


def test_zero_contrast_does_not_radiate():
    pytest.importorskip('grcwa')
    scene = Scene(kind='flat', substrate_n=1., substrate_k=0., particle_n=1.)
    value = project_ray(cylinder_source(), scene, np.array([0., 0., 1.]), reference_halfwidth=1)
    np.testing.assert_array_equal(value, np.zeros(3))

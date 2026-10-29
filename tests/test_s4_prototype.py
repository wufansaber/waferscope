"""Protocol, geometry and optional native-runtime physical checks."""
import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from s4_prototype.lua_scene import script
from s4_prototype.output import parse
from s4_prototype.solver import DEFAULT_EXE, solve


def test_sphere_is_not_allowed_in_1d():
    with pytest.raises(ValueError, match='sphere'):
        script(Scene(), particle=True, one_dimensional=True)


def test_sliced_analytic_volume_converges():
    exact = np.pi*.1**3/6
    errors = []
    for slices in (8, 16, 32):
        source, volume = script(Scene(slices=slices), particle=True)
        assert source.count('SetLayerPatternCircle') == slices
        errors.append(abs(volume/exact-1))
    assert errors[2] < errors[1] < errors[0]
    assert errors[-1] < .001


def test_protocol_rejects_incomplete_data():
    with pytest.raises(ValueError, match='Missing'):
        parse('WS4\tCOUNT\t1\n')


@pytest.mark.skipif(not DEFAULT_EXE.exists(), reason='Private S4 not built')
@pytest.mark.parametrize('polarization', ['s', 'p'])
@pytest.mark.parametrize('angle', [0., 45., 75.])
@pytest.mark.parametrize('index', [1.+0j, 1.5+0j, 5.57+.387j])
@pytest.mark.parametrize('one_dimensional', [False, True])
def test_fresnel_complex_amplitude_and_power(tmp_path, polarization, angle, index, one_dimensional):
    scene = Scene(kind='flat', modes=1, incidence_deg=angle,
                  substrate_n=index.real, substrate_k=index.imag)
    data, info = solve(scene, tmp_path/'trial', polarization=polarization, one_dimensional=one_dimensional,
                      near_points=[(.017, -.023, scene.cap_top_um+.02)])
    theta = np.deg2rad(angle)
    c, s = np.cos(theta), np.sin(theta)
    qt = np.sqrt(index**2-s*s)
    rs = (c-qt)/(c+qt)
    rp = (index**2*c-qt)/(index**2*c+qt)
    expected = np.array([-rs, 0., 0.]) if polarization == 's' else np.array([0., -c*rp, s*rp])
    np.testing.assert_allclose(data['electric_reflected_z0'][0], expected, atol=2e-12)
    np.testing.assert_allclose(info['reflectance'], abs(rs if polarization == 's' else rp)**2, atol=2e-12)
    expected_transmission = (1-info['reflectance'])*np.exp(
        -2*(2*np.pi/scene.wavelength_um)*qt.imag*scene.depth_um)
    np.testing.assert_allclose(info['flux_into_substrate'], expected_transmission, atol=2e-12)
    if not index.imag:
        np.testing.assert_allclose(info['reflectance']+info['flux_into_substrate'], 1., atol=2e-12)
    point = data['near_points_um'][0]
    k0 = 2*np.pi/scene.wavelength_um
    incident = np.array([-1., 0., 0.]) if polarization == 's' else np.array([0., c, s])
    phase = np.exp(1j*k0*s*point[1])
    field = phase*(incident*np.exp(-1j*k0*c*point[2])+expected*np.exp(1j*k0*c*point[2]))
    np.testing.assert_allclose(data['near_electric'][0], field, atol=2e-12)
    hi = np.cross([0., s, -c], incident)
    hr = np.cross([0., s, c], expected)
    h = phase*(hi*np.exp(-1j*k0*c*point[2])+hr*np.exp(1j*k0*c*point[2]))
    np.testing.assert_allclose(data['near_magnetic_impedance_scaled'][0], h, atol=2e-12)


@pytest.mark.skipif(not DEFAULT_EXE.exists(), reason='Private S4 not built')
def test_zero_contrast_particle_does_not_change_complex_reference(tmp_path):
    scene = Scene(modes=31, particle_n=1., particle_k=0.)
    reference, _ = solve(scene, tmp_path/'reference')
    sample, _ = solve(scene, tmp_path/'sample', particle=True)
    np.testing.assert_array_equal(reference['orders'], sample['orders'])
    np.testing.assert_allclose(reference['electric_reflected_z0'], sample['electric_reflected_z0'], atol=2e-12)


@pytest.mark.skipif(not DEFAULT_EXE.exists(), reason='Private S4 not built')
def test_weak_contrast_particle_has_continuous_zero_contrast_limit(tmp_path):
    scene = Scene(modes=31, particle_n=1.)
    reference, _ = solve(scene, tmp_path/'reference')
    changes = []
    for index in (1.001, 1.0001):
        from dataclasses import replace
        sample, _ = solve(replace(scene, particle_n=index), tmp_path/str(index), particle=True)
        changes.append(np.linalg.norm(sample['electric_reflected_z0']-reference['electric_reflected_z0']))
    assert changes[0] > 0
    assert .08 < changes[1]/changes[0] < .12

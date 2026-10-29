"""Physical invariance checks for the separate, explicitly enabled S4 adapter."""
from dataclasses import replace

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from s4_prototype.lua_scene import script
from s4_prototype.solver import DEFAULT_EXE, ROOT, solve

ADAPTER = ROOT/'build/s4_runtime/bin_stripes/S4.exe'
native = pytest.mark.skipif(not ADAPTER.exists(), reason='S4 stripe adapter not built')


def test_rectangular_basis_is_explicit_and_requires_2d():
    scene = Scene(harmonics_x=1, harmonics_y=12)
    text, _ = script(scene, stripe_rule=True)
    assert 'S:SetRectangularBasis(1,12)' in text
    assert 'S:UseAnalyticStripeRule()' in text
    with pytest.raises(ValueError, match='2D'):
        script(scene, one_dimensional=True)


@native
@pytest.mark.parametrize('polarization', ['s', 'p'])
@pytest.mark.parametrize('periods', [1, 2, 3])
def test_2d_stripe_matches_1d_complex_field_at_same_cutoff(tmp_path, polarization, periods):
    primitive_halfwidth = 6
    scene = Scene(modes=2*primitive_halfwidth+1, periods_y=periods)
    points = [(.017, .03, .02), (-.05, .21, -.05)]
    reference, a = solve(scene, tmp_path/'one_d', one_dimensional=True,
                         polarization=polarization, near_points=points)
    rectangle = replace(scene, harmonics_x=1, harmonics_y=primitive_halfwidth*periods)
    actual, b = solve(rectangle, tmp_path/'two_d', executable=ADAPTER,
                      stripe_rule=True, polarization=polarization, near_points=points)
    lookup = {tuple(g): i for i, g in enumerate(actual['orders'])}
    take = [lookup[tuple(g)] for g in reference['orders']]
    np.testing.assert_allclose(actual['electric_reflected_z0'][take],
                               reference['electric_reflected_z0'], atol=2e-9, rtol=2e-9)
    np.testing.assert_allclose(a['reflectance'], b['reflectance'], atol=2e-10)
    keep = (actual['orders'][:, 0] == 0) & (actual['orders'][:, 1] % periods == 0)
    assert actual['reflectance_by_order'][~keep].sum() < 1e-20
    assert b['actual_modes'] == rectangle.requested_modes
    assert b['retained_order_extent'] == [1, primitive_halfwidth*periods]
    np.testing.assert_allclose(actual['near_electric'], reference['near_electric'], atol=2e-9, rtol=2e-9)
    np.testing.assert_allclose(actual['near_magnetic_impedance_scaled'],
                               reference['near_magnetic_impedance_scaled'], atol=2e-9, rtol=2e-9)


@native
def test_off_switch_preserves_unmodified_circular_result(tmp_path):
    scene = Scene(modes=31)
    reference, _ = solve(scene, tmp_path/'baseline')
    actual, _ = solve(scene, tmp_path/'adapter_off', executable=ADAPTER)
    np.testing.assert_array_equal(actual['orders'], reference['orders'])
    np.testing.assert_allclose(actual['electric_reflected_z0'],
                               reference['electric_reflected_z0'], atol=2e-12)


@native
def test_sphere_layer_is_not_replaced_by_stripe_rule(tmp_path):
    scene = Scene(kind='flat', modes=31, slices=4)
    reference, _ = solve(scene, tmp_path/'baseline', particle=True)
    actual, _ = solve(scene, tmp_path/'adapter_on', particle=True,
                      executable=ADAPTER, stripe_rule=True)
    np.testing.assert_allclose(actual['electric_reflected_z0'],
                               reference['electric_reflected_z0'], atol=2e-12)

from dataclasses import asdict

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from tools.rcwa_crosscheck_report import check_methods, field_difference


def infos():
    a = dict(backend='S4', scene=asdict(Scene()), polarization='p', one_dimensional=False,
             factorization='default', stripe_rule=True, epsilon_method='analytic')
    b = dict(backend='grcwa', scene=asdict(Scene(grating_method='analytic_li')), particle_method='analytic')
    return a, b


def test_only_unused_pixel_grids_may_differ_in_analytic_crosscheck():
    a, b = infos()
    b['scene']['grid_x'] = 3
    check_methods(a, b)
    b['scene']['slices'] = 12
    with pytest.raises(ValueError, match='Matched'):
        check_methods(a, b)


def test_fft_particles_cannot_be_mislabeled_as_matched_analytic_geometry():
    a, b = infos()
    b['particle_method'] = 'fft'
    with pytest.raises(ValueError, match='analytic'):
        check_methods(a, b)


def test_s4_backend_cannot_be_replaced_by_grcwa_metadata():
    a, b = infos()
    a['backend'] = 'grcwa'
    with pytest.raises(ValueError, match='relabeled'):
        check_methods(a, b)


def test_matching_indices_do_not_hide_different_ray_directions():
    a = dict(orders=np.array([[0, 0]]), propagating=np.array([True]),
             directions=np.array([[0., 0., 1.]]), electric_reflected_z0=np.array([[1., 0., 0.]], complex))
    b = dict(a, directions=np.array([[.1, 0., np.sqrt(.99)]]))
    with pytest.raises(AssertionError):
        field_difference(a, b, 75.)


def test_equal_power_with_a_phase_difference_is_not_fitted_to_agree():
    a = dict(orders=np.array([[0, 0]]), propagating=np.array([True]),
             directions=np.array([[0., 0., 1.]]), electric_reflected_z0=np.array([[1., 0., 0.]], complex))
    b = dict(a, electric_reflected_z0=a['electric_reflected_z0']*np.exp(.1j))
    difference = field_difference(a, b, 75.)['power_weighted_relative_complex_field_difference']
    assert difference == pytest.approx(abs(1-np.exp(.1j)))

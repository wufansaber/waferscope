"""Convergence claims must use compatible inputs and one-variable steps."""
from copy import deepcopy
from dataclasses import asdict, replace
import json

import numpy as np
import pytest

from rcwa_prototype.geometry import Scene
from s4_prototype.solver import sha256
from tools.s4_convergence_report import compare_cases, compare_joint_cases, load_checked, prepare_pair, same_settings


def pair(scene=None, scale=1.):
    scene = scene or Scene(harmonics_x=3, harmonics_y=40)
    directions = np.array([[.707, 0., np.sqrt(1-.707**2)],
                           [-.707, 0., np.sqrt(1-.707**2)], [0., 0., 1.]])
    electric = scale*np.array([[0., .01, 0.], [0., .01, 0.], [.01, 0., 0.]], complex)
    sample = dict(orders=np.array([[1, 0], [-1, 0], [0, 0]]),
        directions=directions, propagating=np.ones(3, bool), electric_reflected_z0=electric,
        reflectance_by_order=np.full(3, .001))
    reference = dict(sample, electric_reflected_z0=np.zeros((3, 3), complex),
                     reflectance_by_order=np.zeros(3))
    info = dict(backend='S4', scene=asdict(scene), particle=True, polarization='p',
        one_dimensional=False, factorization='default', resolution=8, executable_sha256='test',
        actual_modes=scene.requested_modes, reflectance=.003,
        geometry=dict(sliced_particle_volume_um3=1., exact_particle_volume_um3=1.))
    reference_info = dict(deepcopy(info), particle=False, reflectance=0.)
    return reference_info, reference, info, sample


def test_one_variable_change_uses_complex_field_and_channel_power():
    coarse = prepare_pair(*pair(scale=.99))
    fine = prepare_pair(*pair(Scene(harmonics_x=4, harmonics_y=40)))
    result = compare_cases(coarse, fine, 'harmonics_x')
    assert result['relative_power_weighted_complex_field_change'] == pytest.approx(.01)
    assert result['channel_difference_field_power_relative_changes'] == pytest.approx([.99**2-1]*3)
    assert result['meets_local_tolerance'] is True


def test_phase_change_cannot_be_hidden_by_equal_power():
    coarse = prepare_pair(*pair(scale=np.exp(.1j)))
    fine = prepare_pair(*pair(Scene(harmonics_x=4, harmonics_y=40)))
    result = compare_cases(coarse, fine, 'harmonics_x')
    assert result['relative_power_weighted_complex_field_change'] > .09
    assert result['channel_difference_field_power_relative_changes'] == pytest.approx([0.]*3)
    assert result['meets_local_tolerance'] is False


@pytest.mark.parametrize('changes', [dict(harmonics_x=4, harmonics_y=48),
                                    dict(harmonics_x=4, substrate_n=4.), dict(lx_um=1.)])
def test_rejects_multi_parameter_or_physical_changes(changes):
    scene = Scene(harmonics_x=3, harmonics_y=40)
    with pytest.raises(ValueError, match='exactly one'):
        compare_cases(prepare_pair(*pair(scene)),
                       prepare_pair(*pair(replace(scene, **changes))), 'harmonics_x')


def test_bare_reference_reuse_is_only_allowed_for_slices():
    a, ref, b, sample = pair()
    a['scene']['slices'] = 12
    result = prepare_pair(a, ref, b, sample)
    assert result['summary']['reference_reused_across_slice_count'] is True
    a['scene']['gap_um'] = .02
    with pytest.raises(ValueError, match='geometry'):
        prepare_pair(a, ref, b, sample)


def test_different_solver_binary_cannot_be_called_a_truncation_change():
    _, _, info, _ = pair()
    other = deepcopy(info)
    other['executable_sha256'] = 'other'
    with pytest.raises(ValueError, match='executable'):
        same_settings(info, other)


@pytest.mark.parametrize('tamper', ['input', 'npz'])
def test_raw_data_and_input_hashes_are_verified(tmp_path, tamper):
    source = tmp_path/'input.lua'
    archive = tmp_path/'modes.npz'
    source.write_text('original', encoding='ascii')
    np.savez(archive, electric_reflected_z0=np.ones((1, 3), complex))
    info = dict(raw_npz_sha256=sha256(archive), input_sha256=sha256(source))
    (tmp_path/'report.json').write_text(json.dumps(info), encoding='utf-8')
    assert load_checked(tmp_path)[0] == info
    if tamper == 'input':
        source.write_text('changed', encoding='ascii')
    else:
        np.savez(archive, electric_reflected_z0=np.zeros((1, 3), complex))
    with pytest.raises(ValueError, match='hash mismatch'):
        load_checked(tmp_path)


def test_missing_propagating_orders_cannot_pass_local_tolerance():
    coarse = prepare_pair(*pair())
    fine = prepare_pair(*pair(Scene(harmonics_x=4, harmonics_y=40)))
    coarse['delta']['propagating'] = np.array([False, True, True])
    result = compare_cases(coarse, fine, 'harmonics_x')
    assert result['comparable'] is False
    assert result['meets_local_tolerance'] is False


def test_joint_step_is_explicit_and_does_not_relax_single_variable_checks():
    coarse = prepare_pair(*pair())
    fine = prepare_pair(*pair(Scene(harmonics_x=4, harmonics_y=48, slices=10)))
    result = compare_joint_cases(coarse, fine)
    assert result['parameter'] == 'joint'
    assert set(result['changed_numerical_parameters']) == {'harmonics_x', 'harmonics_y', 'slices'}
    with pytest.raises(ValueError, match='exactly one'):
        compare_cases(coarse, fine, 'harmonics_x')


@pytest.mark.parametrize('changes', [dict(harmonics_x=4, slices=6),
                                    dict(harmonics_x=4, slices=10, diameter_um=.2)])
def test_joint_refinement_rejects_decreasing_or_physical_changes(changes):
    scene = Scene(harmonics_x=3, harmonics_y=40)
    with pytest.raises(ValueError):
        compare_joint_cases(prepare_pair(*pair(scene)), prepare_pair(*pair(replace(scene, **changes))))

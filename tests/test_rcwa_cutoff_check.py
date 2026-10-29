from dataclasses import asdict
import hashlib
import json

import numpy as np
import pytest

from rcwa_prototype.cutoff_check import compare_cutoff
from rcwa_prototype.geometry import Scene


def case(tmp_path, name, mx=4, my=40, scale=1., phase=1., **changes):
    scene = Scene(harmonics_x=mx, harmonics_y=my, lx_um=2/3, **changes)
    folder = tmp_path/name
    source, reference = folder/'source', folder/'reference'
    source.mkdir(parents=True)
    reference.mkdir()
    np.savez_compressed(source/'modes.npz',
        source_electric_fourier=np.full((2, 1, 3), scale*phase, complex),
        source_kxy_per_um=np.array([[0., 0.]]), source_z_um=np.array([.02, .08]),
        source_weight_um=np.array([.05, .05]), source_radius_um=np.array([.03, .03]))
    np.savez_compressed(reference/'modes.npz', electric_reflected_z0=np.array([[1., 0., 0.]], complex),
        reflectance_by_order=np.array([1.]), orders=np.array([[0, 0]]),
        propagating=np.array([True]), directions=np.array([[0., 0., 1.]]))
    provenance = []
    for path, particle in ((source, True), (reference, False)):
        digest = hashlib.sha256((path/'modes.npz').read_bytes()).hexdigest()
        report = dict(scene=asdict(scene), raw_npz_sha256=digest, particle=particle,
            particle_method='analytic', current_z_order=2 if particle else 0,
            backend='grcwa', backend_version='test', analytic_particle_branch='test',
            coordinates='global z=0 reference; no fitted phase')
        (path/'report.json').write_text(json.dumps(report), encoding='utf-8')
        provenance.append(dict(folder=str(path), sha256=digest))
    weights = np.tile(np.array([.1, .2]), (3, 1))
    fields = np.zeros((3, 2, 3), complex)
    fields[:, :, 0] = scale*phase
    cross = np.sum(weights*np.sum(abs(fields)**2, axis=-1), axis=-1)
    np.savez_compressed(folder/'cone.npz', directions=np.tile([0., 0., 1.], (3, 2, 1)),
        solid_angle_weights=weights, far_electric_um=fields, collected_cross_section_um2=cross)
    radiation = dict(scene=asdict(scene), na=.6, tilts_deg=[45., -45., 0.],
        source=provenance[0], reference=provenance[1],
        status='FROZEN_PERIODIC_CURRENT_DIAGNOSTIC_NOT_ISOLATED_DEFECT',
        same_basis_gate=dict(comparable=True, relative_power_weighted_complex_field_change=1e-12,
                             relative_total_power_difference=1e-12),
        levels=[dict(raw='cone.npz', reference_halfwidth=40, cross_section_um2=cross.tolist())])
    (folder/'report.json').write_text(json.dumps(radiation), encoding='utf-8')
    return folder


def test_absolute_amplitude_and_local_pass_do_not_certify_isolation(tmp_path):
    result = compare_cutoff(case(tmp_path, 'a'), case(tmp_path, 'b', mx=5, scale=1.01))
    assert result['local_step_passed']
    assert not result['isolated_defect_certified']
    assert result['refined_axis'] == 'x'
    assert result['channels'][0]['relative_collected_power_change'] == pytest.approx(.0201)
    assert result['channels'][0]['relative_weighted_complex_field_change'] == pytest.approx(.01)
    assert result['internal_source']['internal_e_relative_l2'] == pytest.approx(.01)
    assert result['background_control']['relative_power_weighted_complex_field_change'] == 0.


def test_phase_is_not_fitted(tmp_path):
    result = compare_cutoff(case(tmp_path, 'a'), case(tmp_path, 'b', mx=5, phase=np.exp(.1j)))
    assert not result['local_step_passed']
    assert result['channels'][0]['relative_collected_power_change'] == pytest.approx(0.)
    assert result['channels'][0]['relative_weighted_complex_field_change'] == pytest.approx(abs(np.exp(.1j)-1))


def test_y_only_refinement_is_supported(tmp_path):
    result = compare_cutoff(case(tmp_path, 'a'), case(tmp_path, 'b', my=48))
    assert result['refined_axis'] == 'y'
    assert result['local_step_passed']


@pytest.mark.parametrize('change', [dict(lx_um=.7), dict(gap_um=0.), dict(slices=10), dict(particle_n=1.6)])
def test_noncutoff_changes_are_rejected(tmp_path, change):
    a = case(tmp_path, 'a')
    if 'lx_um' in change:
        b = case(tmp_path, 'b', mx=5)
        for path in (b/'report.json', b/'source/report.json', b/'reference/report.json'):
            report = json.loads(path.read_text())
            report['scene'].update(change)
            path.write_text(json.dumps(report), encoding='utf-8')
    else:
        b = case(tmp_path, 'b', mx=5, **change)
    with pytest.raises(ValueError, match='Fixed cell'):
        compare_cutoff(a, b)


@pytest.mark.parametrize('mx,my', [(4, 40), (3, 40), (5, 48), (5, 39)])
def test_joint_unchanged_or_decreasing_orders_rejected(tmp_path, mx, my):
    with pytest.raises(ValueError, match='exactly one'):
        compare_cutoff(case(tmp_path, 'a'), case(tmp_path, 'b', mx=mx, my=my))


def test_reference_integrity_is_checked(tmp_path):
    a, b = case(tmp_path, 'a'), case(tmp_path, 'b', mx=5)
    path = b/'report.json'
    report = json.loads(path.read_text())
    report['reference']['sha256'] = 'incorrect'
    path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match='Reference hash'):
        compare_cutoff(a, b)


def test_source_solver_version_must_match(tmp_path):
    a, b = case(tmp_path, 'a'), case(tmp_path, 'b', mx=5)
    path = b/'source/report.json'
    report = json.loads(path.read_text())
    report['backend_version'] = 'changed'
    path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match='backend_version'):
        compare_cutoff(a, b)


def test_reference_scene_must_match(tmp_path):
    a, b = case(tmp_path, 'a'), case(tmp_path, 'b', mx=5)
    path = b/'reference/report.json'
    report = json.loads(path.read_text())
    report['scene']['incidence_deg'] = 70.
    path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match='both solved scenes'):
        compare_cutoff(a, b)


@pytest.mark.parametrize('key', ['analytic_particle_branch', 'coordinates'])
def test_branch_and_phase_conventions_must_match(tmp_path, key):
    a, b = case(tmp_path, 'a'), case(tmp_path, 'b', mx=5)
    path = b/'source/report.json'
    report = json.loads(path.read_text())
    report[key] = 'changed'
    path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match=key):
        compare_cutoff(a, b)


def test_report_preserves_values_and_refuses_overwrite(tmp_path):
    from tools.rcwa_cutoff_report import save_report
    a, b = case(tmp_path, 'a'), case(tmp_path, 'b', mx=5, scale=1.1)
    output = tmp_path/'report'
    result = save_report(a, b, output)
    assert not result['local_step_passed']
    assert result['channels'][0]['relative_collected_power_change'] == pytest.approx(.21)
    assert (output/'cutoff_sensitivity.png').stat().st_size > 10000
    assert (output/'COMPARISON.md').is_file()
    assert json.loads((output/'comparison.json').read_text())['channels'] == result['channels']
    with pytest.raises(ValueError, match='empty output'):
        save_report(a, b, output)


def test_nonfinite_reference_fields_are_rejected_even_with_updated_hash(tmp_path):
    a, b = case(tmp_path, 'a'), case(tmp_path, 'b', mx=5)
    raw_path = b/'reference/modes.npz'
    with np.load(raw_path) as raw:
        values = dict(raw)
    values['electric_reflected_z0'] *= np.nan
    np.savez_compressed(raw_path, **values)
    digest = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    for path in (b/'report.json', b/'reference/report.json'):
        report = json.loads(path.read_text())
        if 'levels' in report:
            report['reference']['sha256'] = digest
        else:
            report['raw_npz_sha256'] = digest
        path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match='Finite bare-reference'):
        compare_cutoff(a, b)

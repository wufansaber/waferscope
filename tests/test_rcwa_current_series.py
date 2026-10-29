from dataclasses import asdict
import hashlib
import json

import numpy as np
import pytest

from rcwa_prototype.current_series import build_series
from rcwa_prototype.geometry import Scene


def case(tmp_path, mx, scale=1., phase=1., **changes):
    scene = Scene(harmonics_x=mx, harmonics_y=40, lx_um=mx/6, **changes)
    source = tmp_path/('source_'+str(mx))
    folder = tmp_path/('radiation_'+str(mx))
    source.mkdir()
    folder.mkdir()
    np.savez_compressed(source/'modes.npz',
        source_electric_fourier=np.full((2, 1, 3), scale*phase, complex),
        source_kxy_per_um=np.array([[0., 0.]]), source_z_um=np.array([.02, .08]),
        source_weight_um=np.array([.05, .05]), source_radius_um=np.array([.03, .03]))
    digest = hashlib.sha256((source/'modes.npz').read_bytes()).hexdigest()
    (source/'report.json').write_text(json.dumps(dict(scene=asdict(scene),
        raw_npz_sha256=digest, particle=True, current_z_order=2)), encoding='utf-8')
    weights = np.tile(np.array([.1, .2]), (3, 1))
    fields = np.zeros((3, 2, 3), complex)
    fields[:, :, 0] = scale*phase
    cross = np.sum(weights*np.sum(abs(fields)**2, axis=-1), axis=-1)
    np.savez_compressed(folder/'cone.npz', directions=np.tile([0., 0., 1.], (3, 2, 1)),
        solid_angle_weights=weights, far_electric_um=fields, collected_cross_section_um2=cross)
    report = dict(scene=asdict(scene), na=.6, tilts_deg=[45., -45., 0.],
        source=dict(folder=str(source), sha256=digest),
        status='FROZEN_PERIODIC_CURRENT_DIAGNOSTIC_NOT_ISOLATED_DEFECT',
        same_basis_gate=dict(comparable=True, relative_power_weighted_complex_field_change=1e-12,
                             relative_total_power_difference=1e-12),
        levels=[dict(raw='cone.npz', reference_halfwidth=40, cross_section_um2=cross.tolist())])
    (folder/'report.json').write_text(json.dumps(report), encoding='utf-8')
    return folder


def test_three_cells_keep_absolute_amplitude_changes(tmp_path):
    folders = [case(tmp_path, mx, scale) for mx, scale in ((3, 1.), (4, 1.1), (5, 1.2))]
    result = build_series(folders)
    assert result['adjacent_steps'][0]['channels'][0]['relative_collected_power_change'] == pytest.approx(.21)
    assert not result['last_two_x_steps_passed']
    assert not result['isolated_defect_certified']


def test_two_local_passes_do_not_certify_an_isolated_defect(tmp_path):
    result = build_series([case(tmp_path, mx) for mx in (3, 4, 5)])
    assert result['last_two_x_steps_passed']
    assert not result['isolated_defect_certified']


def test_global_phase_is_not_fitted(tmp_path):
    result = build_series([case(tmp_path, 3), case(tmp_path, 4, phase=1j), case(tmp_path, 5)])
    row = result['adjacent_steps'][0]['channels'][0]
    assert row['relative_collected_power_change'] == pytest.approx(0.)
    assert row['relative_weighted_complex_field_change'] == pytest.approx(np.sqrt(2))


def test_two_cells_and_unordered_cells_are_rejected(tmp_path):
    folders = [case(tmp_path, mx) for mx in (3, 4, 5)]
    with pytest.raises(ValueError, match='At least three'):
        build_series(folders[:2])
    with pytest.raises(ValueError, match='strictly increasing'):
        build_series(folders[::-1])


def test_changed_particle_gap_is_not_a_cell_size_step(tmp_path):
    folders = [case(tmp_path, 3), case(tmp_path, 4, gap_um=0.), case(tmp_path, 5)]
    with pytest.raises(ValueError, match='Only X cell size'):
        build_series(folders)


def test_changed_raw_source_is_rejected(tmp_path):
    folders = [case(tmp_path, mx) for mx in (3, 4, 5)]
    np.savez_compressed(tmp_path/'source_4'/'modes.npz', modified=np.array(1))
    with pytest.raises(ValueError, match='hash'):
        build_series(folders)


def test_stale_integral_is_rejected(tmp_path):
    folders = [case(tmp_path, mx) for mx in (3, 4, 5)]
    path = folders[1]/'report.json'
    report = json.loads(path.read_text())
    report['levels'][0]['cross_section_um2'][0] *= 2
    path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match='Saved collected power'):
        build_series(folders)


def test_failed_modal_gate_is_not_accepted(tmp_path):
    folders = [case(tmp_path, mx) for mx in (3, 4, 5)]
    path = folders[1]/'report.json'
    report = json.loads(path.read_text())
    report['same_basis_gate']['relative_power_weighted_complex_field_change'] = 1e-3
    path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError, match='modal gate'):
        build_series(folders)


def test_changed_phase_is_detected_when_radiation_hash_is_recorded(tmp_path):
    folders = [case(tmp_path, mx) for mx in (3, 4, 5)]
    folder = folders[1]
    path = folder/'report.json'
    report = json.loads(path.read_text())
    report['levels'][0]['raw_sha256'] = hashlib.sha256((folder/'cone.npz').read_bytes()).hexdigest()
    path.write_text(json.dumps(report), encoding='utf-8')
    with np.load(folder/'cone.npz') as raw:
        values = dict(raw)
    values['far_electric_um'] *= 1j
    np.savez_compressed(folder/'cone.npz', **values)
    with pytest.raises(ValueError, match='Radiation raw NPZ hash'):
        build_series(folders)

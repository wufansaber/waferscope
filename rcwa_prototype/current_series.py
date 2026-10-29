"""Auditable X-cell sequences; never certify an isolated defect from X alone."""
import hashlib
import json
from pathlib import Path

import numpy as np

from .current_compare import compare_arrays, compare_folders, load
from .current_radiation import load_source, matching_scenes
from .current_overlap import compare_source_fields
from .geometry import Scene


def checked_level(folder):
    report, data, provenance = load(folder)
    gate = report['same_basis_gate']
    if (report.get('status') != 'FROZEN_PERIODIC_CURRENT_DIAGNOSTIC_NOT_ISOLATED_DEFECT'
            or not gate.get('comparable')):
        raise ValueError('An explicit periodic-current diagnostic with a passing modal gate is required')
    for key, tolerance in (('relative_power_weighted_complex_field_change', 1e-7),
                           ('relative_total_power_difference', 2e-7)):
        value = gate.get(key)
        if value is None or not np.isfinite(value) or abs(value) > tolerance:
            raise ValueError('Same-basis modal gate did not pass')
    _, source_scene, _, source = load_source(report['source']['folder'])
    if source['sha256'] != report['source']['sha256']:
        raise ValueError('The source current has changed since its radiation calculation')
    scene = Scene(**report['scene'])
    if not matching_scenes(scene, source_scene):
        raise ValueError('Radiation scene does not match the solved current')
    rows = compare_arrays(data, data)
    cross = np.array([row['first_cross_section_um2'] for row in rows])
    for saved in (data['collected_cross_section_um2'], report['levels'][-1]['cross_section_um2']):
        if not np.allclose(cross, saved, rtol=2e-13, atol=0., equal_nan=False):
            raise ValueError('Saved collected power does not match the raw complex fields')
    provenance['report_sha256'] = hashlib.sha256((Path(folder)/'report.json').read_bytes()).hexdigest()
    return report, scene, cross, provenance


def build_series(folders, power_tolerance=.05, field_tolerance=.02):
    if len(folders) < 3:
        raise ValueError('At least three separately solved increasing X cells are required')
    if not all(np.isfinite(v) and 0 < v < 1 for v in (power_tolerance, field_tolerance)):
        raise ValueError('Finite relative tolerances in (0,1) required')
    reports, scenes, values, provenance = zip(*(checked_level(folder) for folder in folders))
    if any(b.lx_um <= a.lx_um for a, b in zip(scenes, scenes[1:])):
        raise ValueError('X cell lengths must be strictly increasing')
    rows = []
    for i in range(1, len(folders)):
        compared = compare_folders(folders[i-1], folders[i])
        first_source = load_source(reports[i-1]['source']['folder'])[0]
        second_source = load_source(reports[i]['source']['folder'])[0]
        internal = compare_source_fields(first_source, second_source, scenes[i].eps_particle-1)
        channel_rows = compared['channels']
        passed = all(row['relative_collected_power_change'] is not None and
            abs(row['relative_collected_power_change']) < power_tolerance and
            row['relative_weighted_complex_field_change'] is not None and
            row['relative_weighted_complex_field_change'] < field_tolerance for row in channel_rows)
        rows.append(dict(first_index=i-1, second_index=i, channels=channel_rows, internal_source=internal,
                         local_step_passed=bool(passed)))
    return dict(status='X_CELL_SEQUENCE_NOT_ISOLATED_CERTIFICATION', isolated_defect_certified=False,
        physics_scope='Periodic currents radiated into the same continuous cones; no fitted amplitude or phase.',
        cells=[dict(lx_um=scene.lx_um, ly_um=scene.ly_um,
                    harmonics_x=scene.harmonics_x, harmonics_y=scene.harmonics_y,
                    physical_cutoff_cycles_um=[scene.harmonics_x/scene.lx_um,
                                               scene.harmonics_y/scene.ly_um],
                    collected_cross_section_um2=cross.tolist(), provenance=origin)
               for scene, cross, origin in zip(scenes, values, provenance)],
        na=reports[0]['na'], tilts_deg=reports[0]['tilts_deg'],
        reference_halfwidth=reports[0]['levels'][-1]['reference_halfwidth'],
        tolerances=dict(power_relative=power_tolerance, complex_field_relative=field_tolerance),
        adjacent_steps=rows, last_two_x_steps_passed=bool(all(row['local_step_passed'] for row in rows[-2:])),
        limitations=['An X sequence alone does not test Y-cell separation.',
            'Fixed maximum Fourier frequency does not establish Fourier-cutoff convergence.',
            'The last-step differences are sensitivities, not rigorous error bounds.',
            'No quantitative sample/reference image or camera SNR acceptance from this diagnostic.'])


def save_report(folders, output):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use an empty output directory to preserve prior comparisons')
    result = build_series(folders)
    output.mkdir(parents=True, exist_ok=True)
    (output/'series.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    return result

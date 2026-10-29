"""Single-axis Fourier-cutoff control on a fixed periodic particle cell."""
from dataclasses import replace

import numpy as np

from .current_compare import compare_arrays, load
from .comparison import field_change
from .current_overlap import compare_source_fields
from .current_radiation import load_source, matching_scenes
from .current_series import checked_level


def _checked_source(report, scene):
    source, actual_scene, info, _ = load_source(report['source']['folder'])
    reference, reference_scene, reference_info, reference_origin = load_source(report['reference']['folder'])
    if (reference_origin['sha256'] != report['reference']['sha256']
            or not matching_scenes(scene, actual_scene)
            or not matching_scenes(scene, reference_scene)):
        raise ValueError('Reference hash and both solved scenes must match the radiation report')
    if (not info.get('particle') or reference_info.get('particle')
            or info.get('particle_method') != 'analytic' or not info.get('current_z_order')):
        raise ValueError('Analytic particle current and matching particle-free reference required')
    if (not np.all(np.isfinite(reference['electric_reflected_z0']))
            or not np.all(np.isfinite(reference['reflectance_by_order']))):
        raise ValueError('Finite bare-reference complex fields and powers required')
    return source, info, reference


def compare_cutoff(baseline, refined, power_tolerance=.05, field_tolerance=.02):
    if not all(np.isfinite(v) and 0 < v < 1 for v in (power_tolerance, field_tolerance)):
        raise ValueError('Finite relative tolerances in (0,1) required')
    ra, sa, cross_a, pa = checked_level(baseline)
    rb, sb, cross_b, pb = checked_level(refined)
    if (sa.harmonics_x is None or sb.harmonics_x is None
            or not matching_scenes(sa, replace(sb, harmonics_x=sa.harmonics_x,
                                               harmonics_y=sa.harmonics_y))):
        raise ValueError('Fixed cell, geometry, materials and explicit rectangular bases required')
    orders_a = np.array([sa.harmonics_x, sa.harmonics_y])
    orders_b = np.array([sb.harmonics_x, sb.harmonics_y])
    step = orders_b-orders_a
    if np.any(step < 0) or np.count_nonzero(step) != 1:
        raise ValueError('Increase exactly one Fourier halfwidth, with the other unchanged')
    if ra['na'] != rb['na'] or ra['tilts_deg'] != rb['tilts_deg']:
        raise ValueError('Matching collection cones required')
    if ra['levels'][-1]['reference_halfwidth'] != rb['levels'][-1]['reference_halfwidth']:
        raise ValueError('Matching reciprocal-background truncation required')
    source_a, info_a, reference_a = _checked_source(ra, sa)
    source_b, info_b, reference_b = _checked_source(rb, sb)
    for key in ('backend', 'backend_version', 'particle_method', 'current_z_order',
                'analytic_particle_branch', 'coordinates'):
        if key not in info_a or key not in info_b or info_a[key] != info_b[key]:
            raise ValueError('Matching source solver metadata required: '+key)
    _, field_a, _ = load(baseline)
    _, field_b, _ = load(refined)
    channels = compare_arrays(field_a, field_b)
    passed = all(row['relative_collected_power_change'] is not None and
        abs(row['relative_collected_power_change']) < power_tolerance and
        row['relative_weighted_complex_field_change'] is not None and
        row['relative_weighted_complex_field_change'] < field_tolerance for row in channels)

    def describe(scene, cross, origin, report):
        return dict(harmonics=[scene.harmonics_x, scene.harmonics_y], modes=scene.requested_modes,
            physical_cutoff_cycles_um=[scene.harmonics_x/scene.lx_um, scene.harmonics_y/scene.ly_um],
            cross_section_um2=cross.tolist(), provenance=origin, reference=report['reference'],
            same_basis_gate=report['same_basis_gate'])

    background = field_change(reference_a, reference_b, sa.incidence_deg)
    background.update(baseline_reflectance=float(reference_a['reflectance_by_order'].sum()),
        refined_reflectance=float(reference_b['reflectance_by_order'].sum()),
        norm_denominator='Refined bare-reference field, unlike the baseline denominator for continuous channels.',
        scope='Bare-reference control only. X-invariant grating excitation is unchanged by X cutoff; '
              'Y cutoff can change the background. Neither is an isolated-particle validation.')
    return dict(status='SINGLE_AXIS_CUTOFF_SENSITIVITY_NOT_FULL_CONVERGENCE',
        refined_axis='x' if step[0] else 'y', cell_um=[sa.lx_um, sa.ly_um],
        baseline=describe(sa, cross_a, pa, ra), refined=describe(sb, cross_b, pb, rb),
        channels=channels, na=ra['na'], tilts_deg=ra['tilts_deg'],
        background_control=background,
        internal_source=compare_source_fields(source_a, source_b, sa.eps_particle-1),
        tolerances=dict(power_relative=power_tolerance, complex_field_relative=field_tolerance),
        local_step_passed=bool(passed), isolated_defect_certified=False,
        limitations=['One single-axis step is not joint or asymptotic Fourier convergence.',
            'Periodic cell, geometry, slicing and source z quadrature are held fixed.',
            'No fitted amplitude or global phase; angular quadrature remains finite.',
            'Do not subtract scalar cutoff and cell sensitivities to infer image coupling.',
            'No isolated-defect image, calibrated material or camera SNR acceptance.'])

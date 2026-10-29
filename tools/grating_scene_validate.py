"""Bounded native S4 cutoff/slicing controls on the same periodic scene."""
import argparse
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import psutil

from fullwave_prototype.guarded_run import run_guarded
from grating_scene.pipeline import DEFAULT_EXECUTABLE, digest, load_scene, write_json
from grating_scene.diagnostics import scene_diagnostics
from rcwa_prototype.geometry import matrix_memory_estimate
from rcwa_prototype.current_radiation import matching_scenes


def compare_scenes(baseline, refined):
    ca, ra, da = load_scene(baseline)
    cb, rb, db = load_scene(refined)
    if not matching_scenes(ca.scene, replace(cb.scene, harmonics_x=ca.scene.harmonics_x,
            harmonics_y=ca.scene.harmonics_y, slices=ca.scene.slices)):
        raise ValueError('Only modal cutoff and sphere slicing may change in this fixed-cell control')
    options_a, options_b = ca.to_dict(), cb.to_dict()
    for key in ('scene', 'image_samples', 'near_enabled'):
        options_a.pop(key)
        options_b.pop(key)
    if options_a != options_b:
        raise ValueError('Identical illumination, collection, observation and relay controls required')
    a, b = da['far'], db['far']
    ta = {tuple(g): i for i, g in enumerate(a['orders']) if a['propagating'][i]}
    tb = {tuple(g): i for i, g in enumerate(b['orders']) if b['propagating'][i]}
    if ta.keys() != tb.keys():
        raise ValueError('All propagating orders must be represented in both truncations')
    keys = sorted(ta)
    ia, ib = [ta[k] for k in keys], [tb[k] for k in keys]
    np.testing.assert_allclose(a['directions'][ia], b['directions'][ib], atol=1e-12, rtol=1e-12)
    np.testing.assert_array_equal(a['channel_masks'][:, ia], b['channel_masks'][:, ib])
    ea, eb = a['reflected_electric_delta'][ia], b['reflected_electric_delta'][ib]
    weights = a['directions'][ia, 2]/np.cos(np.deg2rad(ca.scene.incidence_deg))
    channels = []
    for i, mask in enumerate(a['channel_masks'][:, ia]):
        na = float(np.sum(weights[mask, None]*abs(ea[mask])**2))
        nb = float(np.sum(weights[mask, None]*abs(eb[mask])**2))
        delta = float(np.sum(weights[mask, None]*abs(eb[mask]-ea[mask])**2))
        channels.append(dict(channel=i+1, enabled=ca.numerical_apertures[i] > 0,
            baseline_difference_field_efficiency=na, refined_difference_field_efficiency=nb,
            power_change=nb/na-1 if na else None,
            weighted_complex_field_change=float(np.sqrt(delta/na)) if na else None))
    near = {}
    for plane in ('xy', 'xz', 'yz') if da['near'] and db['near'] else ():
        np.testing.assert_allclose(da['near'][plane+'_points_um'], db['near'][plane+'_points_um'], atol=1e-13)
        x, y = da['near'][plane+'_delta_electric'], db['near'][plane+'_delta_electric']
        norm = np.sum(abs(x)**2)
        near[plane] = float(np.sqrt(np.sum(abs(y-x)**2)/norm)) if norm else None
    mirror_checks = dict(baseline=scene_diagnostics(a, ra['native'], ca)['mirror_check_passed'],
                         refined=scene_diagnostics(b, rb['native'], cb)['mirror_check_passed'])
    passed = all(mirror_checks.values()) and any(row['enabled'] for row in channels) and all(not row['enabled'] or (row['power_change'] is not None and
        abs(row['power_change']) < .05 and row['weighted_complex_field_change'] < .02) for row in channels)
    return dict(baseline=str(Path(baseline).resolve()), refined=str(Path(refined).resolve()),
        baseline_report_sha256=digest(Path(baseline)/'scene_report.json'),
        refined_report_sha256=digest(Path(refined)/'scene_report.json'),
        channels=channels, mirror_checks_passed=mirror_checks,
        near_delta_pointwise_relative_l2=near, local_far_step_passed=bool(passed),
        isolated_defect_certified=False,
        scope='Same periodic cell, discrete-order field/power changes; near metric is a sampled-plane norm, not stored energy.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=900.)
    parser.add_argument('--with-near-fields', action='store_true')
    parser.add_argument('--cancel-file', type=Path)
    parser.add_argument('--executable', type=Path, default=DEFAULT_EXECUTABLE)
    parser.add_argument('--max-memory-gib', type=float, default=8.)
    parser.add_argument('--memory-profile', type=Path)
    args = parser.parse_args()
    config, native_report, _ = load_scene(args.baseline)
    if any(not np.isfinite(v) or v <= 0 for v in (args.timeout, args.max_memory_gib)):
        parser.error('Positive finite timeout and memory cap required')
    if any(digest(args.executable) != native_report['native'][key]['executable_sha256']
           for key in ('sample', 'reference')):
        parser.error('Use the same S4 executable as the saved scene for a cutoff/slicing control')
    root = args.output
    if root.exists() and any(root.iterdir()):
        parser.error('Use an empty validation directory')
    root.mkdir(parents=True, exist_ok=True)
    s = config.scene
    cases = dict(x_cutoff=replace(s, harmonics_x=s.harmonics_x+1),
        y_cutoff=replace(s, harmonics_y=s.harmonics_y+8),
        slices=replace(s, slices=2*s.slices),
        joint=replace(s, harmonics_x=s.harmonics_x+1, harmonics_y=s.harmonics_y+8, slices=2*s.slices))
    report = dict(baseline=str(args.baseline.resolve()), status='running', cases={},
        targets=dict(power_relative=.05, complex_field_relative=.02), isolated_defect_certified=False)
    for name, scene in cases.items():
        samples = max(config.image_samples, 2**(4*max(scene.harmonics_x, scene.harmonics_y)).bit_length())
        trial = replace(config, scene=scene, image_samples=samples,
                        near_enabled=config.near_enabled if args.with_near_fields else False)
        trial.validate()
        path = root/(name+'.input.json')
        write_json(path, trial.to_dict())
        estimate = matrix_memory_estimate(scene.requested_modes, scene.slices+1)
        if args.memory_profile:
            from grating_scene.resources import calibrated_estimate
            estimate = calibrated_estimate(args.memory_profile, args.executable,
                scene.requested_modes, scene.slices+1, trial.image_samples)
        available = psutil.virtual_memory().available/2**30
        row = dict(planning=estimate, available_gib=available, system_reserve_gib=3.,
                   maximum_gib=args.max_memory_gib)
        if estimate['planning_gib'] > min(args.max_memory_gib, available-3.):
            row['status'] = 'blocked_memory_preflight'
        else:
            resource = run_guarded([sys.executable, '-m', 'grating_scene.run', '--worker',
                '--config', str(path), '--output', str(root/name), '--timeout', str(args.timeout),
                '--executable', str(args.executable)],
                root/(name+'.resources.json'), maximum=args.max_memory_gib, reserve=3., timeout=args.timeout,
                cancel_file=args.cancel_file)
            row['status'] = resource['status']
            if resource['status'] == 'completed':
                row['comparison'] = compare_scenes(args.baseline, root/name)
        report['cases'][name] = row
        write_json(root/'validation_report.json', report)
        print('VALIDATION_CASE '+name+' '+row['status'], flush=True)
        if row['status'] == 'stopped_cancelled':
            break
    report['status'] = ('cancelled' if any(r['status'] == 'stopped_cancelled' for r in report['cases'].values())
                        else 'completed_controls_not_isolated_certification')
    report['all_cases_solved'] = len(report['cases']) == len(cases) and all(
        r['status'] == 'completed' for r in report['cases'].values())
    report['all_local_far_steps_passed'] = report['all_cases_solved'] and all(
        r['comparison']['local_far_step_passed'] for r in report['cases'].values())
    write_json(root/'validation_report.json', report)
    print('VALIDATION_FINISHED', flush=True)


if __name__ == '__main__':
    main()

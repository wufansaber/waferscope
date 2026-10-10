"""Resource-bounded one-factor haze scans with paired numerical error estimates."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from .far_field import angular_grid, radiation_amplitude
from .guarded_run import run_guarded
from .haze_statistics import confidence_interval, statistical_report, estimated_sample_budget
from .haze_sample_cache import completed_sample
from .parameters import Parameters
from .rough_surface import Roughness, taper_effective_area


AXES = ('fem_mesh_scale', 'pml_thickness_um', 'pml_strength', 'domain_radius_um',
        'height_grid_points', 'patch_width_um', 'taper_width_um', 'angular_step_deg',
        'quadrature_order', 'fem_order', 'rough_mesh_scale', 'pml_mesh_scale', 'pml_order')
PAIRED_AXES = ('fem_mesh_scale', 'pml_thickness_um', 'pml_strength', 'domain_radius_um',
               'angular_step_deg', 'quadrature_order', 'fem_order', 'rough_mesh_scale',
               'pml_mesh_scale', 'pml_order')


class SampleLimitError(RuntimeError):
    def __init__(self, path, status):
        super().__init__('Sample '+status+': '+str(path))
        self.resource_path, self.status = str(path), status


def default_plan(parameters, roughness, samples=8, apertures=(.5, .5, .5), tilts=(45, -45, 0)):
    spacing = roughness.width_um/(roughness.grid_points-1)
    taper = roughness.width_um*roughness.taper_fraction
    plan = dict(parameters=asdict(parameters), roughness=asdict(roughness), samples=samples,
        numerical_apertures=list(apertures), tilt_degrees=list(tilts),
        comparison_area_um2=.04, angular_step_deg=5, quadrature_order=5,
        tolerance=.05, target_relative_ci=.05,
        scans=dict(fem_mesh_scale=[1., .85, .7],
            pml_thickness_um=[parameters.pml_thickness_um*.75,
                              parameters.pml_thickness_um, parameters.pml_thickness_um*1.25],
            pml_strength=[parameters.pml_strength*.75, parameters.pml_strength,
                          parameters.pml_strength*1.25],
            domain_radius_um=[parameters.domain_radius_um,
                              parameters.domain_radius_um*1.15, parameters.domain_radius_um*1.3],
            height_grid_points=sorted(set([roughness.grid_points,
                min(41, roughness.grid_points+4), min(41, roughness.grid_points+8)])),
            patch_width_um=[roughness.width_um, roughness.width_um*1.1, roughness.width_um*1.2],
            taper_width_um=sorted(set([max(spacing, taper*.75), taper, min(roughness.width_um/2, taper*1.25)])),
            angular_step_deg=[10, 5, 2.5], quadrature_order=[5, 7]))
    if roughness.model == 'measured_height':
        plan['samples'] = 1
        plan['scans']['patch_width_um'] = [roughness.width_um*.8,
                                          roughness.width_um*.9, roughness.width_um]
    if parameters.order < 3:
        plan['scans']['fem_order'] = [parameters.order, parameters.order+1]
    plan['scans']['rough_mesh_scale'] = [1., .9, .8]
    plan['scans']['pml_mesh_scale'] = [parameters.pml_mesh_scale,
                                     parameters.pml_mesh_scale*.85]
    if parameters.order > 1 and parameters.pml_order not in (0, parameters.order):
        plan['scans']['pml_order'] = [parameters.pml_order, parameters.order]
    return plan


def case_settings(plan, axis=None, value=None):
    p, r = Parameters(**plan['parameters']), Roughness(**plan['roughness'])
    angular, quadrature = plan.get('angular_step_deg', 5), plan.get('quadrature_order', 5)
    if axis == 'fem_mesh_scale':
        p = replace(p, air_mesh_um=p.air_mesh_um*value,
                    substrate_mesh_um=p.substrate_mesh_um*value)
        r = replace(r, mesh_um=r.mesh_um*value)
    elif axis == 'rough_mesh_scale':
        r = replace(r, mesh_um=r.mesh_um*value)
    elif axis == 'pml_mesh_scale':
        p = replace(p, pml_mesh_scale=value)
    elif axis == 'pml_order':
        p = replace(p, pml_order=int(value))
    elif axis in ('pml_thickness_um', 'pml_strength', 'domain_radius_um'):
        p = replace(p, **{axis: value})
    elif axis == 'height_grid_points':
        r = replace(r, grid_points=int(value))
        if r.model in ('gaussian', 'abc'):
            r = replace(r, sampling_grid_levels=tuple(sorted(set(
                [r.grid_points]+[int(n) for n in plan['scans'][axis]]))))
    elif axis == 'fem_order':
        p = replace(p, order=int(value))
    elif axis == 'patch_width_um':
        count = max(4, int(np.ceil(value/(r.width_um/(r.grid_points-1))))+1)
        r = replace(r, width_um=value, grid_points=count,
                    taper_fraction=r.width_um*r.taper_fraction/value)
    elif axis == 'taper_width_um':
        r = replace(r, taper_fraction=value/r.width_um)
    elif axis == 'angular_step_deg':
        angular = value
    elif axis == 'quadrature_order':
        quadrature = value
    p.validate()
    r.validate(p, particle_enabled=False)
    if r.model not in ('gaussian', 'abc', 'measured_psd', 'measured_height') or r.rms_nm <= 0:
        raise ValueError('Convergence requires nonzero statistical or measured roughness')
    if p.beam_waist_um != 0:
        raise ValueError('Use plane or Gaussian footprint mode, not legacy beam waist')
    if not np.isfinite(angular) or angular not in (2.5, 5, 10):
        raise ValueError('Angular steps must be 2.5, 5 or 10 degrees')
    if quadrature not in (5, 7):
        raise ValueError('Source quadrature must be 5 or 7')
    return p, r, angular, quadrature


def validate_plan(plan):
    count = plan.get('samples', 8)
    if not isinstance(count, int) or not 1 <= count <= 128:
        raise ValueError('Use 1..128 independent surfaces per convergence case')
    if plan['roughness']['model'] == 'measured_height' and count != 1:
        raise ValueError('A measured height map is deterministic: use samples=1, not repeated seeds')
    for key in ('comparison_area_um2', 'tolerance', 'target_relative_ci'):
        value = plan.get(key, .04 if key == 'comparison_area_um2' else .05)
        if not np.isfinite(value) or value <= 0 or (key != 'comparison_area_um2' and value >= 1):
            raise ValueError(key+' is outside its valid range')
    case_settings(plan)
    for axis, levels in plan.get('scans', {}).items():
        if axis not in AXES or not 2 <= len(levels) <= 6 or len(set(levels)) != len(levels):
            raise ValueError('Scan axes require 2..6 unique levels: '+axis)
        for level in levels:
            if not np.isfinite(level) or level <= 0 or (axis in ('height_grid_points', 'fem_order', 'pml_order') and int(level) != level):
                raise ValueError('Invalid scan level for '+axis)
            case_settings(plan, axis, level)
    if not plan.get('scans'):
        raise ValueError('At least one convergence axis is required')
    if plan['roughness']['seed']+count-1 > 1000000:
        raise ValueError('Last seed exceeds 1000000')


def compare_cases(first, second, axis, tolerance, deterministic=False):
    a, b = np.array(first['statistics']['sample_observables_um2']), np.array(
        second['statistics']['sample_observables_um2'])
    same_surface = first['surface_sha256'] == second['surface_sha256']
    parents = first.get('realization_sha256', [])
    same_realization = bool(parents and all(parents) and parents == second.get('realization_sha256'))
    paired = (axis in PAIRED_AXES and same_surface) or (axis == 'height_grid_points' and same_realization)
    differences = {}
    for i, name in enumerate(('hemisphere', 'ch1', 'ch2', 'ch3')):
        denominator = abs(float(b[:, i].mean()))
        delta = float(b[:, i].mean()-a[:, i].mean())
        interval = confidence_interval(b[:, i]-a[:, i]) if paired else None
        half = ((interval['ci95_um2'][1]-delta) if interval and interval['ci95_um2'] else
                0. if deterministic and paired else None)
        bound = (abs(delta)+half)/denominator if half is not None and denominator > 0 else None
        inactive = name.startswith('ch') and second['statistics']['numerical_apertures'][i-1] == 0
        if inactive:
            bound = 0.
        differences[name] = dict(delta_um2=delta,
            relative_change=abs(delta)/denominator if denominator > 0 else None,
            paired_delta_ci95_um2=interval['ci95_um2'] if interval else None,
            relative_error_upper_bound=bound,
            target_met=bool(bound is not None and bound <= tolerance and
                            (len(a) >= 3 or deterministic or inactive)))
    return dict(paired_same_surface=bool(paired and same_surface),
        paired_same_realization=bool(paired and same_realization), observables=differences,
        status=('within_tested_tolerance' if all(v['target_met'] for v in differences.values()) else
                'unresolved'),
        scope='Last two supplied levels only; no universal accuracy certification. '
              'Height refinement is paired only for an exact shared Gaussian-process parent; '
              'patch and taper remain unpaired.')


def _sample(parameters, roughness, root, quadrature, memory_gib, timeout_seconds):
    measurement = None
    if roughness.model.startswith('measured_'):
        from .experimental_surface import input_fingerprint
        measurement = input_fingerprint(roughness.measurement_path)
    request = dict(parameters=asdict(parameters), roughness=asdict(roughness),
                   measurement_sha256=measurement, angular_step_deg=5, quadrature_order=quadrature)
    # Source-code identity prevents stale resumption after model implementation changes.
    code_hash = hashlib.sha256()
    source_names = ('parameters.py', 'rough_surface.py', 'experimental_surface.py',
        'rough_statistics.py', 'solver.py', 'absorbing_boundary.py', 'field_space.py',
        'planar.py', 'illumination.py', 'volume_sources.py', 'far_field.py', 'haze_ensemble.py',
        'coupled_surface.py', 'rough_geometry.py', 'finite_illumination.py')
    for source in sorted(Path(__file__).parent/name for name in source_names):
        if not source.is_file():
            continue
        code_hash.update(source.name.encode('ascii'))
        code_hash.update(source.read_bytes())
    request['solver_code_sha256'] = code_hash.hexdigest()
    signature = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    folder = root/'samples'/signature[:20]
    if completed_sample(folder, request):
        return folder
    folder.mkdir(parents=True, exist_ok=True)
    # Invalidate old success before replacing inputs, even if the guard cannot start.
    (folder/'resources.json').write_text(json.dumps(dict(status='pending')), encoding='utf-8')
    (folder/'input.json').write_text(json.dumps(request, indent=2), encoding='utf-8')
    result = run_guarded([sys.executable, '-m', 'fullwave_prototype.haze_ensemble', 'sample',
        '--input', str(folder/'input.json'), '--seed', str(roughness.seed), '--output', str(folder)],
        folder/'resources.json', maximum=memory_gib, timeout=timeout_seconds)
    if result['status'] != 'completed':
        if result['status'] in ('stopped_memory_limit', 'stopped_system_memory_reserve', 'stopped_timeout', 'failed'):
            raise SampleLimitError(folder/'resources.json', result['status'])
        raise RuntimeError('Convergence sample stopped: '+str(folder/'resources.json'))
    return folder


def run_plan(plan, output, memory_gib=10, timeout_seconds=1800, sample_cache=None):
    validate_plan(plan)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    cache = Path(sample_cache).resolve() if sample_cache is not None else output
    (output/'convergence_inputs.json').write_text(json.dumps(plan, indent=2), encoding='utf-8')
    report = dict(status='running', axes={}, baseline=None, inputs=plan, sample_cache=str(cache),
                  physical_convergence='unverified',
                  method='One-factor scans; shared seeds and exact-height pairing for numerical axes; '
                         'joint Gaussian-process parents for Gaussian/ABC height refinement; '
                         'intensity averaging; sequential resource-bounded samples; cached source regridding.')
    def save():
        (output/'haze_convergence_report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    def calculate(axis=None, value=None):
        p, rough, step, quadrature = case_settings(plan, axis, value)
        maps, sample_info, resources, sample_folders = [], [], [], []
        theta, phi, directions = angular_grid(step)
        for offset in range(plan['samples']):
            print('Convergence {}={} | sample {}/{}'.format(axis or 'baseline', value,
                  offset+1, plan['samples']), flush=True)
            r = replace(rough, seed=rough.seed+offset)
            folder = _sample(p, r, cache, quadrature, memory_gib, timeout_seconds)
            sample_folders.append(str(folder))
            with np.load(folder/'haze_sources.npz', allow_pickle=False) as data:
                amplitude = radiation_amplitude(data['points'], data['weights'], directions.reshape(-1, 3), p)
            maps.append(np.sum(abs(amplitude.reshape(directions.shape[:2]+(2,)))**2, axis=-1))
            sample_info.append(json.loads((folder/'sample_report.json').read_text(encoding='utf-8')))
            resources.append(json.loads((folder/'resources.json').read_text(encoding='utf-8')))
        scale = plan.get('comparison_area_um2', .04)/taper_effective_area(rough)
        if rough.model == 'measured_height' or p.beam_profile == 'gaussian':
            scale = 1.
        statistics = statistical_report(maps, theta, phi, scale,
            plan.get('numerical_apertures', (.5, .5, .5)), plan.get('tilt_degrees', (45, -45, 0)),
            plan.get('target_relative_ci', .05))
        statistics['sample_budget_estimate'] = estimated_sample_budget(
            statistics, plan.get('target_relative_ci', .05))
        from .illumination_output import illumination_metadata
        return dict(level=value, parameters=asdict(p), roughness=asdict(rough),
            illumination=illumination_metadata(p, rough), area_extrapolation_applied=scale != 1.,
            angular_step_deg=step, quadrature_order=quadrature, statistics=statistics,
            sample_folders=sample_folders,
            surface_sha256=[item['surface_sha256'] for item in sample_info],
            realization_sha256=[item.get('realization_sha256') for item in sample_info],
            maximum_dofs=max(item['solver']['dofs'] for item in sample_info),
            mean_sample_elapsed_seconds=float(np.mean([item['elapsed_seconds'] for item in sample_info])),
            peak_accounted_gib=max(item['peak_accounted_gib'] for item in resources))
    save()
    try:
        report['baseline'] = calculate()
        save()
        for axis, levels in plan['scans'].items():
            report['axes'][axis] = dict(cases=[], limits=[], comparison=None)
            for level in levels:
                try:
                    report['axes'][axis]['cases'].append(calculate(axis, level))
                except SampleLimitError as exc:
                    report['axes'][axis]['limits'].append(dict(level=level, status=exc.status,
                        resources=exc.resource_path,
                        note='Whole case excluded; partial successful surfaces would bias the ensemble'))
                save()
            cases = report['axes'][axis]['cases']
            if len(cases) >= 2 and not report['axes'][axis]['limits']:
                report['axes'][axis]['comparison'] = compare_cases(cases[-2], cases[-1], axis,
                    plan.get('tolerance', .05), deterministic=plan['roughness']['model'] == 'measured_height')
            else:
                report['axes'][axis]['comparison'] = dict(status=(
                    'failed_levels' if any(v['status'] == 'failed' for v in report['axes'][axis]['limits'])
                    else 'resource_limited'),
                    note='Not all requested levels completed; no convergence claim')
            save()
        report['status'] = ('complete_with_limits' if any(v['limits'] for v in report['axes'].values())
                            else 'complete')
    except BaseException as exc:
        report.update(status='interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed', error=str(exc))
        save()
        raise
    save()
    from .haze_convergence_plotting import plot_convergence
    plot_convergence(output, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='JSON convergence plan')
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-memory-gib', type=float, default=10)
    parser.add_argument('--sample-timeout-seconds', type=float, default=1800)
    parser.add_argument('--sample-cache', help='Study root whose completed samples may be reused')
    args = parser.parse_args()
    run_plan(json.loads(Path(args.input).read_text(encoding='utf-8')), args.output,
             args.max_memory_gib, args.sample_timeout_seconds, args.sample_cache)


if __name__ == '__main__':
    main()

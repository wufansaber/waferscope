"""Serial, resource-guarded NGSolve thread trials on a reproducible rough mesh."""
import argparse
from dataclasses import replace, asdict
import json
from pathlib import Path
import sys
import time

import numpy as np


def worker(source, output, threads):
    import ngsolve as ng
    from .parameters import Parameters
    from .rough_surface import Roughness, build_rough_mesh
    from .solver import solve, evaluate_fields
    from .volume_sources import polarization_sources
    from .far_field import radiation_amplitude

    original = json.loads(Path(source).read_text(encoding='utf-8'))
    p = replace(Parameters(**original['parameters']), threads=threads)
    p.validate()
    rough = Roughness(**original['roughness'])
    # Mesh construction is deliberately single-threaded in every trial.
    ng.SetNumThreads(1)
    start = time.perf_counter()
    mesh = build_rough_mesh(p, rough)
    mesh_seconds = time.perf_counter()-start
    _, field, metrics = solve(p, mesh=mesh)
    angles = np.deg2rad(np.arange(-60, 61, 15))
    probes = np.column_stack((.13*np.sin(angles), np.full(len(angles), .035),
                             .13*np.cos(angles)))
    e, h, _ = evaluate_fields(mesh, field, probes, p)
    start = time.perf_counter()
    points, sources, _ = polarization_sources(mesh, field, p, True, 7)
    directions = np.column_stack((np.sin(angles), np.zeros(len(angles)), np.cos(angles)))
    amplitude = radiation_amplitude(points, sources, directions, p)
    metrics.update(mesh_seconds=mesh_seconds, radiation_seconds=time.perf_counter()-start)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'fields.npz', e=e, h=h, amplitude=amplitude)
    (output/'result.json').write_text(json.dumps(dict(parameters=asdict(p), roughness=asdict(rough),
        metrics=metrics), indent=2), encoding='utf-8')


def summarize(output):
    output = Path(output)
    trials = []
    baseline = None
    identity = None
    for directory in sorted(output.glob('trial_*')):
        resource = json.loads((directory/'resources.json').read_text(encoding='utf-8'))
        entry = dict(trial=directory.name, resources=resource)
        if resource['status'] == 'completed':
            result = json.loads((directory/'result.json').read_text(encoding='utf-8'))
            entry.update(result)
            with np.load(directory/'fields.npz') as data:
                fields = {key: data[key] for key in data.files}
            if baseline is None:
                baseline = fields
                identity = tuple(result['metrics'][k] for k in ('elements', 'dofs', 'factored_dofs'))
            entry['same_mesh_sizes'] = identity == tuple(result['metrics'][k]
                for k in ('elements', 'dofs', 'factored_dofs'))
            entry['field_relative_l2'] = {key: float(np.linalg.norm(value-baseline[key]) /
                max(np.linalg.norm(baseline[key]), 1e-300)) for key, value in fields.items()}
        trials.append(entry)
    groups = []
    completed = [t for t in trials if 'parameters' in t]
    valid = bool(trials) and len(completed) == len(trials) and all(
        t['same_mesh_sizes'] and all(np.isfinite(v) and v < 1e-8 for v in t['field_relative_l2'].values())
        and np.isfinite(t['metrics']['solve_seconds']) and t['metrics']['solve_seconds'] > 0
        for t in completed)
    for threads in sorted({t['parameters']['threads'] for t in completed}):
        group = [t for t in completed if t['parameters']['threads'] == threads]
        times = [t['metrics']['solve_seconds'] for t in group]
        groups.append(dict(threads=threads, repeats=len(times), median_seconds=float(np.median(times)),
            min_seconds=min(times), max_seconds=max(times),
            peak_accounted_gib=max(t['resources']['peak_accounted_gib'] for t in group),
            worst_field_relative_l2=max(max(t['field_relative_l2'].values()) for t in group)))
    recommendation = None
    if valid and groups and all(g['repeats'] >= 2 for g in groups):
        fastest = min(g['median_seconds'] for g in groups)
        recommendation = min(g['threads'] for g in groups if g['median_seconds'] <= 1.05*fastest)
    summary = dict(trials=trials, groups=groups, correctness_passed=valid,
        recommended_threads=recommendation,
        selection='Fewest threads within 5% of fastest median; >=2 trials each and relative field error <1e-8.',
        scope='One with-particle solve per process; same single-thread mesh '
        'construction. BLAS/OpenMP limited to one; NGSolve threads vary. No concurrent trials.')
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, help='Existing scattering report.json')
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, nargs='+', default=(1, 2, 4, 8))
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    if args.repeats < 1 or any(t < 1 for t in args.threads):
        parser.error('Positive threads and repeats required')
    if args.worker:
        if len(args.threads) != 1:
            parser.error('Worker requires exactly one thread count')
        worker(args.source, args.output, args.threads[0])
        return
    from .guarded_run import run_guarded
    output = Path(args.output)
    if output.exists() and any(output.iterdir()):
        parser.error('Use an empty output directory to avoid mixing trials')
    for repeat in range(args.repeats):
        # Reverse the second pass to reduce systematic warm-up/order bias.
        counts = args.threads if repeat % 2 == 0 else args.threads[::-1]
        for threads in counts:
            trial = output/'trial_{}_{}t'.format(repeat+1, threads)
            result = run_guarded([sys.executable, '-m', __package__+'.thread_benchmark',
                '--worker', '--source', args.source, '--output', str(trial),
                '--threads', str(threads)], trial/'resources.json', maximum=12, reserve=3,
                timeout=900)
            summarize(output)
            if result['status'] != 'completed':
                raise RuntimeError('Benchmark stopped: '+result['status'])


if __name__ == '__main__':
    main()

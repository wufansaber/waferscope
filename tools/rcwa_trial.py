"""Guarded experimental periodic RCWA trial, isolated from the UI runtime."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from rcwa_prototype.geometry import Scene, matrix_memory_estimate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--particle', action='store_true')
    parser.add_argument('--near', action='store_true')
    parser.add_argument('--kind', choices=('flat', 'grating'), default='grating')
    parser.add_argument('--modes', type=int, default=121)
    parser.add_argument('--slices', type=int, default=8)
    parser.add_argument('--lx', type=float, default=.5)
    parser.add_argument('--periods-y', type=int, default=2)
    parser.add_argument('--grid', type=int, nargs=2, default=(128, 256))
    parser.add_argument('--factorization', choices=('pol', 'laurent'), default='laurent')
    parser.add_argument('--max-memory-gib', type=float, default=8.)
    parser.add_argument('--timeout', type=float, default=600.)
    args = parser.parse_args()
    scene = Scene(kind=args.kind, modes=args.modes, slices=args.slices, lx_um=args.lx,
        periods_y=args.periods_y, grid_x=args.grid[0], grid_y=args.grid[1], factorization=args.factorization)
    scene.validate()
    if args.worker:
        from rcwa_prototype.solver import solve
        args.output.mkdir(parents=True, exist_ok=True)
        result, info = solve(scene, args.particle, args.near)
        info['memory_planning'] = matrix_memory_estimate(scene.modes, info['patterned_layers'])
        np.savez_compressed(args.output/'modes.npz', **result)
        (args.output/'report.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
    else:
        import psutil
        from fullwave_prototype.guarded_run import run_guarded
        if not np.isfinite(args.max_memory_gib) or args.max_memory_gib <= 0:
            parser.error('Memory limit must be finite and positive')
        if not np.isfinite(args.timeout) or args.timeout <= 0:
            parser.error('Timeout must be finite and positive')
        estimate = matrix_memory_estimate(scene.modes, scene.slices+1 if args.particle else 1)
        available = psutil.virtual_memory().available/2**30
        budget = min(args.max_memory_gib, available-3.)
        if estimate['planning_gib'] > budget:
            parser.error('Planning allowance {:.2f} GiB exceeds {:.2f} GiB budget '
                         '(3 GiB system reserve). Reduce modes/layers.'.format(
                             estimate['planning_gib'], max(0., budget)))
        result = run_guarded([sys.executable, '-m', 'tools.rcwa_trial', '--worker']+sys.argv[1:],
            args.output.with_suffix('.resources.json'), maximum=args.max_memory_gib,
            reserve=3., timeout=args.timeout)
        if result['status'] != 'completed':
            raise SystemExit(1)


if __name__ == '__main__':
    main()

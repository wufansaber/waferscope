"""Run an isolated S4 trial with process-tree memory and timeout monitoring."""
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
    parser.add_argument('--one-dimensional', action='store_true')
    parser.add_argument('--kind', choices=('flat', 'grating'), default='grating')
    parser.add_argument('--modes', type=int, default=121)
    parser.add_argument('--harmonics', type=int, nargs=2, metavar=('MX', 'MY'))
    parser.add_argument('--stripe-rule', action='store_true')
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--slices', type=int, default=8)
    parser.add_argument('--lx', type=float, default=.5)
    parser.add_argument('--periods-y', type=int, default=2)
    parser.add_argument('--factorization', choices=('normal', 'default', 'vector', 'jones'), default='default',
                        help='All 2D formulations remain diagnostic until grating/particle convergence passes')
    parser.add_argument('--resolution', type=int, default=8)
    parser.add_argument('--epsilon-method', choices=('analytic', 'fft'), default='analytic')
    parser.add_argument('--polarization', choices=('s', 'p'), default='p')
    parser.add_argument('--max-memory-gib', type=float, default=8.)
    parser.add_argument('--timeout', type=float, default=600.)
    args = parser.parse_args()
    scene = Scene(kind=args.kind, modes=args.modes, slices=args.slices,
                  lx_um=args.lx, periods_y=args.periods_y,
                  harmonics_x=args.harmonics[0] if args.harmonics else None,
                  harmonics_y=args.harmonics[1] if args.harmonics else None)
    scene.validate()
    if args.particle and args.one_dimensional:
        parser.error('A sphere is 3D and requires a 2D supercell')
    if args.worker:
        from s4_prototype.solver import DEFAULT_EXE, solve
        _, info = solve(scene, args.output, particle=args.particle,
            polarization=args.polarization, one_dimensional=args.one_dimensional,
            factorization=args.factorization, resolution=args.resolution, timeout=args.timeout,
            epsilon_method=args.epsilon_method, stripe_rule=args.stripe_rule,
            executable=args.executable or DEFAULT_EXE)
        print(json.dumps({key: info[key] for key in (
            'actual_modes', 'elapsed_seconds', 'reflectance', 'flux_into_substrate',
            'field_vs_library_power_absolute_difference')}, indent=2))
        return
    import psutil
    from fullwave_prototype.guarded_run import run_guarded
    if not np.isfinite(args.max_memory_gib) or args.max_memory_gib <= 0:
        parser.error('Memory limit must be finite and positive')
    if not np.isfinite(args.timeout) or args.timeout <= 0:
        parser.error('Timeout must be finite and positive')
    estimate = matrix_memory_estimate(scene.requested_modes, scene.slices+1 if args.particle else 1)
    budget = min(args.max_memory_gib, psutil.virtual_memory().available/2**30-3.)
    if estimate['planning_gib'] > budget:
        parser.error('Planning allowance exceeds memory budget; reduce modes/layers')
    result = run_guarded([sys.executable, '-m', 'tools.s4_trial', '--worker']+sys.argv[1:],
        args.output.with_suffix('.resources.json'), maximum=args.max_memory_gib,
        reserve=3., timeout=args.timeout)
    if result['status'] != 'completed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()

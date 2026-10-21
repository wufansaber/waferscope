"""Run bounded, reproducible FEM probes for an infinite periodic grating."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np


def worker(args):
    from fullwave_prototype.periodic_grating import GratingInputs, solve, extract
    p = GratingInputs(ridge_width_um=args.width)
    mesh, field, hbg, dhbg, inv_eps, info = solve(p, args.order, args.factor, args.pml,
        args.top, args.bottom, args.bonus)
    arrays, physics = extract(p, mesh, field, hbg, dhbg, inv_eps,
        heights=(args.top*.4, args.top*.65), bottom=-args.bottom+1e-7)
    args.output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output/'modes.npz', **arrays)
    report = dict(inputs=asdict(p), solver=info, physics=physics,
        scope='Independent 2D TM benchmark, not validation of the finite 3D H(curl) solver.',
        normalization='Each efficiency is reflected normal power / incident normal power per period and unit line length.')
    (args.output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--order', type=int, choices=range(1, 7), default=3)
    parser.add_argument('--factor', type=float, default=1.)
    parser.add_argument('--pml', type=float, default=.3)
    parser.add_argument('--top', type=float, default=.4)
    parser.add_argument('--bottom', type=float, default=.3)
    parser.add_argument('--bonus', type=int, default=8)
    parser.add_argument('--width', type=float, default=.25)
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([sys.executable, '-m', 'tools.periodic_grating_audit', '--worker']+
            sys.argv[1:], args.output.with_suffix('.resources.json'), maximum=8, reserve=3, timeout=600)
        if result['status'] != 'completed' or result['returncode']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()

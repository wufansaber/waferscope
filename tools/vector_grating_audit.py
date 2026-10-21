"""Check the existing 3D vector solver on the NIST periodic grating case."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np


def worker(args):
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.periodic_grating import GratingInputs
    from fullwave_prototype.periodic_vector import build_cell, extract_modes
    from fullwave_prototype.solver import solve
    g = GratingInputs()
    p = Parameters(wavelength_nm=266, gap_nm=10, domain_radius_um=1.5,
        order=args.order, pml_order=args.pml_order, air_mesh_um=.065*args.factor,
        substrate_mesh_um=.03*args.factor, pml_mesh_scale=args.pml_scale,
        pml_thickness_um=args.pml, direct_solver='umfpack', integration_bonus=10,
        static_condensation=True, max_dofs=500000)
    mesh, space = build_cell(g, p, line_period=args.line_period, corner_mesh=args.corner_mesh)
    if args.probe:
        print(json.dumps(dict(dofs=space.ndof, elements=mesh.ne)), flush=True)
        return
    _, field, info = solve(p, mesh, particle_enabled=False, space=space)
    arrays, checks = extract_modes(g, p, mesh, field, line_period=args.line_period)
    from fullwave_prototype.periodic_vector import volume_modes, energy_check
    arrays.update(volume_modes(g, p, mesh, field, args.line_period))
    checks.update(energy_check(g, p, mesh, field, arrays, args.line_period))
    args.output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output/'modes.npz', **arrays)
    (args.output/'report.json').write_text(json.dumps(dict(inputs=asdict(g),
        parameters=asdict(p), line_period_um=args.line_period, corner_mesh_um=args.corner_mesh,
        solver=info, checks=checks,
        scope='Same 3D vector solver on an infinite periodic grating. Not a finite-pattern or particle validation.'), indent=2), encoding='utf-8')
    print(json.dumps(checks), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--probe', action='store_true')
    parser.add_argument('--order', type=int, default=2)
    parser.add_argument('--pml-order', type=int, default=1)
    parser.add_argument('--factor', type=float, default=1.)
    parser.add_argument('--pml', type=float, default=.3)
    parser.add_argument('--line-period', type=float, default=.08)
    parser.add_argument('--pml-scale', type=float, default=1.8)
    parser.add_argument('--corner-mesh', type=float)
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([sys.executable, '-m', 'tools.vector_grating_audit', '--worker']+
            sys.argv[1:], args.output.with_suffix('.resources.json'), maximum=8, reserve=3, timeout=600)
        if result['status'] != 'completed' or result['returncode']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()

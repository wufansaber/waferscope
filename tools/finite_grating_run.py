"""Bounded, checkpointed 3D finite-trench FEM; not an infinite periodic model."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--output', type=Path, required=True)
    result.add_argument('--worker', action='store_true')
    result.add_argument('--replay', action='store_true', help='Replay only your own trusted local checkpoint')
    result.add_argument('--preview-only', action='store_true')
    result.add_argument('--probe', action='store_true')
    result.add_argument('--count', type=int, default=2)
    for name, default in [('pitch-nm', 500), ('top-width-nm', 250), ('bottom-width-nm', 250),
                          ('depth-nm', 100), ('length-um', .8), ('margin-um', .15),
                          ('air-height-um', .35), ('substrate-depth-um', .18),
                          ('corner-mesh-um', 0), ('factor', 1.), ('pml-um', .2)]:
        result.add_argument('--'+name, type=float, default=default)
    result.add_argument('--order', type=int, choices=(1, 2, 3), default=2)
    result.add_argument('--pml-order', type=int, default=0,
                        help='0 keeps the physical order; any reduction must be validated explicitly')
    result.add_argument('--source-order', type=int, choices=(7, 9, 11), default=9)
    result.add_argument('--surface-orders', type=int, nargs='+', default=[48, 96])
    result.add_argument('--disk-samples', type=int, default=48)
    result.add_argument('--near-samples', type=int, default=100)
    result.add_argument('--max-memory-gib', type=float, default=8.)
    result.add_argument('--timeout-seconds', type=float, default=900.)
    return result


def inputs(args):
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.finite_grating import FiniteGrating
    g = FiniteGrating(trench_count=args.count, pitch_um=args.pitch_nm/1000,
        top_width_um=args.top_width_nm/1000, bottom_width_um=args.bottom_width_nm/1000,
        depth_um=args.depth_nm/1000, length_um=args.length_um, margin_um=args.margin_um,
        air_height_um=args.air_height_um, substrate_depth_um=args.substrate_depth_um,
        corner_mesh_um=args.corner_mesh_um)
    g.validate()
    if not np.isfinite(args.factor) or args.factor <= 0:
        raise ValueError('Positive finite mesh factor required')
    p = Parameters(wavelength_nm=266., domain_radius_um=max(1.5, g.half_width_um),
        order=args.order, pml_order=args.pml_order,
        air_mesh_um=.09*args.factor, substrate_mesh_um=.035*args.factor,
        pml_thickness_um=args.pml_um, pml_mesh_scale=1.5, integration_bonus=10,
        static_condensation=args.order > 1, direct_solver='sparsecholesky', max_dofs=650000)
    p.validate()
    return p, g


def worker(args):
    from fullwave_prototype.finite_grating import build_mesh, FiniteGrating
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.grating_checkpoint import save, load_local
    from fullwave_prototype.grating_plotting import plot_geometry, plot_fields
    from fullwave_prototype.grating_analysis import cap_checks, near_planes, far_disk
    import ngsolve as ng
    if (args.near_samples < 16 or args.disk_samples < 16 or not args.surface_orders or
            min(args.surface_orders) < 4 or args.surface_orders != sorted(set(args.surface_orders))):
        raise ValueError('At least 16 image samples and increasing surface orders >= 4 required')
    args.output.mkdir(parents=True, exist_ok=True)
    if args.replay:
        mesh, field, checkpoint = load_local(args.output)
        p, g = Parameters(**checkpoint['parameters']), FiniteGrating(**checkpoint['geometry'])
        info = checkpoint['solver']
    else:
        p, g = inputs(args)
        plot_geometry(args.output, g)
        if args.preview_only:
            (args.output/'geometry.json').write_text(json.dumps(asdict(g), indent=2), encoding='utf-8')
            return
        ng.SetNumThreads(p.threads)
        mesh = build_mesh(p, g)
        if args.probe:
            from fullwave_prototype.field_space import make_space
            space = make_space(mesh, p)
            print(json.dumps(dict(elements=mesh.ne, dofs=space.ndof,
                factored_dofs=sum(space.FreeDofs(coupling=p.static_condensation)))), flush=True)
            return
        from fullwave_prototype.solver import solve
        _, field, info = solve(p, mesh, particle_enabled=False)
    ng.SetNumThreads(p.threads)
    info = dict(info, defect_geometry=dict(shape='finite_etched_grating',
        parameters=asdict(g), expected_volume_um3=g.trench_volume_um3, particle_present=False))
    if not args.replay:
        save(args.output, mesh, field, p, g, info)
    plot_geometry(args.output, g)
    data, checks, points, weights = cap_checks(mesh, field, p, g, args.source_order, args.surface_orders)
    report = dict(parameters=asdict(p), geometry=asdict(g), solver=info, checks=checks,
        scope='Finite 3D bare etched grating; plane illumination; no particle or random roughness.',
        material_caveat='Legacy Si n=5.57 k=0.387 retained for numerical comparisons, not calibrated at 266 nm.',
        power_caveat='DSC is differential cross section, not a finite-beam collection efficiency.',
        certification='Independent extraction diagnostics only; mesh and PML convergence still required.')
    (args.output/'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    np.savez_compressed(args.output/'cap_checks.npz', **data)
    np.savez_compressed(args.output/'volume_sources.npz', points=points, weights=weights)
    data.update(near_planes(mesh, field, p, g, args.near_samples))
    data.update(far_disk(points, weights, p, args.disk_samples))
    np.savez_compressed(args.output/'fields.npz', **data)
    plot_fields(args.output, data, report)
    print(json.dumps(checks, indent=2), flush=True)


def main():
    args = parser().parse_args()
    if args.worker:
        worker(args)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([sys.executable, '-m', 'tools.finite_grating_run', '--worker']+sys.argv[1:],
            args.output.with_suffix('.resources.json'), maximum=args.max_memory_gib,
            reserve=3., timeout=args.timeout_seconds)
        if result['status'] != 'completed':
            raise SystemExit(1)


if __name__ == '__main__':
    main()

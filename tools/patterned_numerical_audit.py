"""Resource-bounded CH3 convergence probes; no fitted intensity or phase."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys

import numpy as np


def relative(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Finite matching fields required for a convergence comparison')
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-300))


def worker(args):
    import ngsolve as ng
    from tools.patterned_wafer_case import build_parser, case_inputs
    from fullwave_prototype.patterned_wafer import build_mesh, mirror_sources
    from fullwave_prototype.solver import solve
    from fullwave_prototype.volume_sources import polarization_sources
    from fullwave_prototype.radiation import radiation_amplitude
    from fullwave_prototype.direct_imaging import polar_cap
    base = build_parser().parse_args(['--pitch-nm', '500', '--line-width-nm', '250',
        '--line-length-um', '.8', '--half-width-um', '.55', '--air-height-um', '.35',
        '--magnification', '40', '--na', '.6', '--exposure-ms', '.001'])
    p, pattern, settings = case_inputs(base)
    p = replace(p, order=args.order, pml_order=args.pml_order,
        static_condensation=args.order > 1, air_mesh_um=.09*args.factor,
        substrate_mesh_um=.035*args.factor, particle_mesh_um=.015*args.factor,
        pml_thickness_um=args.pml, pml_mesh_scale=args.pml_scale,
        max_dofs=1500000)
    pattern = replace(pattern, half_width_um=.55+args.expansion,
        air_height_um=.35+args.expansion, substrate_depth_um=.18+args.expansion)
    ng.SetNumThreads(2)
    args.output.mkdir(parents=True, exist_ok=True)
    mesh = build_mesh(p, pattern, True, patterned=not args.flat)
    if args.probe:
        from fullwave_prototype.field_space import make_space
        space = make_space(mesh, p)
        print(json.dumps(dict(elements=mesh.ne, dofs=space.ndof,
            factored_dofs=sum(space.FreeDofs(coupling=p.static_condensation)))), flush=True)
        return
    _, area, directions = polar_cap(.6, 0, 4, 12)
    saved = dict(directions=directions, angular_weights=area)
    rows = []
    for bonus in args.bonus:
        q = replace(p, integration_bonus=bonus)
        row = dict(integration_bonus=bonus, states={})
        for enabled, name in ((False, 'without'), (True, 'with')):
            if args.flat and not enabled:
                continue
            print('Solving order={}, bonus={}, {}'.format(q.order, bonus, name), flush=True)
            _, field, info = solve(q, mesh, particle_enabled=enabled)
            for order in (args.source_order, args.source_order+2):
                print('Checking volume radiation quadrature {}'.format(order), flush=True)
                points, weights, _ = polarization_sources(mesh, field, q, enabled, order=order)
                points, weights = mirror_sources(points, weights)
                amplitude = radiation_amplitude(points, weights, directions, q)
                saved['b{}_{}_q{}'.format(bonus, name, order)] = amplitude
                if args.save_sources and bonus == args.bonus[-1] and order == args.source_order+2:
                    np.savez_compressed(args.output/('sources_'+name+'.npz'), points=points, weights=weights)
                del points, weights
            low, high = [saved['b{}_{}_q{}'.format(bonus, name, order)]
                         for order in (args.source_order, args.source_order+2)]
            info['quadrature_relative_l2'] = relative(low, high)
            info['cap_dsc_um2'] = float(np.sum(area/directions[:, 2]*np.sum(abs(high)**2, axis=1)))
            if args.huygens:
                from fullwave_prototype.reciprocity_surface import box_quadrature, amplitude, evaluate_half_field
                bounds = [[-(pattern.line_length_um/2+pattern.half_width_um)/2,
                            (pattern.line_length_um/2+pattern.half_width_um)/2],
                          [-(pattern.pitch_um/2+pattern.trench_width_um/2+pattern.half_width_um)/2,
                            (pattern.pitch_um/2+pattern.trench_width_um/2+pattern.half_width_um)/2],
                          [-(pattern.trench_depth_um+pattern.substrate_depth_um)/2,
                            ((q.diameter_nm+q.gap_nm)/1000+pattern.air_height_um)/2]]
                checks = []
                previous = None
                for quad in (32, 48, 72):
                    surface, normals, areas = box_quadrature(bounds, quad)
                    electric, magnetic = evaluate_half_field(mesh, field, surface, q)
                    surface_amplitude = amplitude(surface, normals, areas, electric, magnetic, directions, q)
                    checks.append(dict(order=quad, volume_vs_surface_relative_l2=relative(surface_amplitude, high),
                        surface_change_relative_l2=relative(previous, surface_amplitude) if previous is not None else None))
                    saved['b{}_{}_surface_q{}'.format(bonus, name, quad)] = surface_amplitude
                    previous = surface_amplitude
                info['reciprocity_surface_checks'] = checks
                info['reciprocity_box_um'] = bounds
            row['states'][name] = info
            del field
        rows.append(row)
        report = dict(parameters=asdict(q), pattern=asdict(pattern), settings=asdict(settings),
            flat=args.flat, rows=rows, scope='CH3 selected-direction convergence probe, not a certification')
        (args.output/'audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        np.savez_compressed(args.output/'amplitudes.npz', **saved)
        print(json.dumps(row), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--flat', action='store_true')
    parser.add_argument('--order', type=int, default=1)
    parser.add_argument('--pml-order', type=int, default=0)
    parser.add_argument('--bonus', type=int, nargs='+', default=[0, 6, 10])
    parser.add_argument('--factor', type=float, default=1.)
    parser.add_argument('--pml', type=float, default=.2)
    parser.add_argument('--pml-scale', type=float, default=1.8)
    parser.add_argument('--expansion', type=float, default=0.)
    parser.add_argument('--source-order', type=int, choices=(7, 9, 11), default=9)
    parser.add_argument('--save-sources', action='store_true')
    parser.add_argument('--probe', action='store_true')
    parser.add_argument('--huygens', action='store_true')
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        result = run_guarded([sys.executable, '-m', 'tools.patterned_numerical_audit',
            '--worker']+sys.argv[1:], args.output.with_suffix('.resources.json'),
            maximum=9., reserve=3., timeout=900.)
        if result['status'] != 'completed' or result['returncode']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()

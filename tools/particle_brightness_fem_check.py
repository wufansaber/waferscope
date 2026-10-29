"""Recompute selected actual FEM spheres; never relabel the archived line as a sphere."""
import argparse
import json
import math
from pathlib import Path

from fullwave_prototype.guarded_run import run_guarded


def mesh_arguments(diameter, scale):
    if not math.isfinite(diameter) or diameter <= 0:
        raise ValueError('Positive finite sphere diameter required')
    if not math.isfinite(scale) or not .4 <= scale <= 1:
        raise ValueError('Mesh scale must be finite and in [0.4,1]')
    return ['--mesh-scale', str(scale), '--particle-mesh-um', str(min(.01, diameter/6000))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--worker-python', required=True)
    parser.add_argument('--diameters', type=float, nargs='+', default=(10., 50., 100.))
    parser.add_argument('--order', type=int, choices=(1, 2, 3), default=2)
    parser.add_argument('--pml-mesh-scale', type=float, default=1.)
    parser.add_argument('--domain-radius-um', type=float)
    parser.add_argument('--rough-width-um', type=float)
    parser.add_argument('--pml-order', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--max-memory-gib', type=float, default=8.)
    parser.add_argument('--incidence-deg', type=float)
    parser.add_argument('--pml-power', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--pml-strength', type=float, default=2.)
    parser.add_argument('--mesh-scale', type=float, default=1.,
                        help='Refine air/substrate mesh without changing particle-region mesh')
    args = parser.parse_args()
    for diameter in args.diameters:
        mesh_arguments(diameter, args.mesh_scale)
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Use an empty output root to preserve prior FEM fields')
    request = json.loads((args.source_run/'request.json').read_text(encoding='utf-8'))
    saved_args = request['arguments']
    imaging_args = saved_args[saved_args.index('imaging_run')+1:]
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for diameter in args.diameters:
        folder = args.output/('d{:g}'.format(diameter))
        override = ['--output', str(folder), '--defect-shape', 'sphere', '--diameter-nm', str(diameter),
            '--beam-profile', 'plane', '--illumination-fit', 'manual', '--beam-waist-um', '0',
            '--wavelength-nm', '266', '--gap-nm', '0', '--rms-nm', '0', '--order', str(args.order),
            '--accuracy', 'balanced', '--direct-solver', 'sparsecholesky',
            '--pml-strength', str(args.pml_strength), '--pml-power', str(args.pml_power),
            '--pml-mesh-scale', str(args.pml_mesh_scale),
            '--pml-order', str(args.pml_order), '--integration-bonus', '8',
            '--pupil-samples', '64', '--padding', '4', '--archive-mode', 'full',
            '--pupil-cache', str(folder/'pupil_cache')]+mesh_arguments(diameter, args.mesh_scale)
        for name in ('domain-radius-um', 'rough-width-um'):
            value = getattr(args, name.replace('-', '_'))
            if value is not None:
                override.extend(('--'+name, str(value)))
        if args.incidence_deg is not None:
            settings = json.loads((args.source_run/'imaging_report.json').read_text(encoding='utf-8'))['settings']
            power = (settings['irradiance_w_um2']*settings['illumination_length_um']*
                     settings['illumination_width_um']*math.cos(math.radians(args.incidence_deg)))
            override.extend(('--incidence-deg', str(args.incidence_deg), '--illumination-power-w', str(power)))
        command = [args.worker_python, '-m', 'fullwave_prototype.imaging_run']+imaging_args+override
        result = run_guarded(command, args.output/('d{:g}.resources.json'.format(diameter)),
                             maximum=args.max_memory_gib, reserve=3., timeout=1800.,
                             cancel_file=args.output/'cancel')
        records.append(dict(diameter_nm=diameter, result=result))
        (args.output/'runs.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
        if result['status'] != 'completed':
            raise SystemExit(1)


if __name__ == '__main__':
    main()

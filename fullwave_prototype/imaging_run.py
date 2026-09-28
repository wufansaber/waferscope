"""Run the near/far/image pipeline or re-image saved complex volume sources."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from .parameters import Parameters
from .rough_surface import Roughness
from .imaging import calculate_images
from .numerical_settings import numerical_parameters
from .imaging_options import add_options, imaging_settings


def image_saved_sources(source, output, settings, archive_mode='full', cache_directory=None):
    source, output = Path(source), Path(output)
    if archive_mode not in ('compact', 'full'):
        raise ValueError('Unknown archive mode')
    settings.validate()
    scattering = json.loads((source/'report.json').read_text(encoding='utf-8'))
    if scattering.get('source_format_version') != 1:
        raise ValueError('This result lacks reusable volume sources. Run the pipeline again without --source.')
    with np.load(source/'paired_fields.npz', allow_pickle=False) as saved:
        data = {key: saved[key] for key in saved.files}
    p = Parameters(**scattering['parameters'])
    arrays, report = calculate_images(data, p, settings, retain_complex_fields=archive_mode == 'full',
                                      cache_directory=cache_directory)
    report.update(scattering_parameters=asdict(p), roughness=scattering['roughness'],
                  source_directory=str(source.resolve()),
                  status='experimental full pipeline; selected Mie/BV/Airy benchmarks tested separately; '
                         'this rough-surface result is not a certified convergence result')
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'image_fields.npz', **arrays)
    (output/'imaging_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    from .imaging_plotting import plot_pipeline
    plot_pipeline(output, data, arrays, scattering, report)
    if settings.detector.enabled:
        from .detector_plot import plot_detector
        plot_detector(output, arrays)
    print(json.dumps(report, indent=2), flush=True)
    return report


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/fullwave_imaging/default')
    parser.add_argument('--source', help='Reuse a rough_run directory with source_format_version=1')
    add_options(parser)
    parser.add_argument('--beam-waist-um', type=float)
    parser.add_argument('--beam-order', type=int)
    parser.add_argument('--na', type=float, nargs=3, default=(.5, .5, .5))
    parser.add_argument('--tilts', type=float, nargs=3, default=(45, -45, 0))
    parser.add_argument('--magnification', type=float, default=20)
    parser.add_argument('--objective-focal-mm', type=float, default=20)
    parser.add_argument('--irradiance-w-um2', type=float, default=1e-6)
    parser.add_argument('--defocus-um', type=float, default=0)
    parser.add_argument('--transmission', type=float, default=1)
    parser.add_argument('--analyzer', choices=('none', 'x', 'y', '45', '135'), default='none')
    parser.add_argument('--pupil-samples', type=int, default=64)
    parser.add_argument('--padding', type=int, default=8)
    parser.add_argument('--near-diagnostics', action='store_true')
    parser.add_argument('--order', type=int)
    parser.add_argument('--threads', type=int, help='NGSolve worker threads; BLAS remains limited by the resource guard')
    parser.add_argument('--accuracy', choices=('preview', 'balanced', 'refined'),
                        help='preview (CLI default): order 2; balanced: order 3 with tested PML tradeoff; refined: original order 3')
    parser.add_argument('--pml-power', type=int, choices=(1, 2, 3))
    parser.add_argument('--pml-mesh-scale', type=float)
    parser.add_argument('--pml-strength', type=float)
    parser.add_argument('--integration-bonus', type=int)
    parser.add_argument('--pml-order', type=int)
    parser.add_argument('--pml-inner-fraction', type=float)
    parser.add_argument('--local-mesh-scale', type=float)
    parser.add_argument('--direct-solver', choices=('umfpack', 'sparsecholesky'))
    parser.add_argument('--mesh-scale', type=float)
    parser.add_argument('--rms-nm', type=float)
    parser.add_argument('--gap-nm', type=float)
    for name in ('wavelength-nm', 'diameter-nm', 'particle-n', 'particle-k', 'substrate-n',
                 'substrate-k', 'incidence-deg', 'domain-radius-um'):
        parser.add_argument('--'+name, type=float)
    parser.add_argument('--polarization', choices=('s', 'p'))
    parser.add_argument('--roughness-seed', type=int)
    return parser


def scattering_settings(args):
    scale = args.mesh_scale if args.mesh_scale is not None else 1
    local = args.local_mesh_scale if args.local_mesh_scale is not None else 1
    if not np.isfinite(local) or not .4 <= local <= 1:
        raise ValueError('local_mesh_scale must be in [0.4, 1]')
    physical = {key: getattr(args, key) for key in ('wavelength_nm', 'diameter_nm', 'particle_n',
                'particle_k', 'substrate_n', 'substrate_k', 'incidence_deg', 'domain_radius_um', 'polarization',
                'beam_waist_um', 'beam_order')}
    p = numerical_parameters(args.accuracy or 'preview', scale, order=args.order, threads=args.threads,
        gap_nm=args.gap_nm if args.gap_nm is not None else 10,
        pml_power=args.pml_power, pml_mesh_scale=args.pml_mesh_scale,
        pml_strength=args.pml_strength, integration_bonus=args.integration_bonus,
        pml_order=args.pml_order, pml_inner_fraction=args.pml_inner_fraction,
        particle_mesh_um=.015*scale*local,
        direct_solver=args.direct_solver, **physical)
    roughness = Roughness(rms_nm=args.rms_nm if args.rms_nm is not None else 2,
                         seed=args.roughness_seed if args.roughness_seed is not None else 7,
                         mesh_um=.025*scale*local)
    roughness.validate(p)
    from .illumination import illumination_rays
    illumination_rays(p)
    return p, roughness


def main():
    parser = build_parser()
    args = parser.parse_args()
    if args.source and args.near_diagnostics:
        parser.error('Dense near diagnostics require a new FEM solve, not --source.')
    if args.source and any(v is not None for v in (args.order, args.threads, args.mesh_scale, args.rms_nm, args.gap_nm,
                                                  args.accuracy, args.pml_power, args.pml_mesh_scale,
                                                  args.pml_strength, args.integration_bonus,
                                                  args.pml_order, args.pml_inner_fraction, args.local_mesh_scale,
                                                  args.direct_solver, args.wavelength_nm, args.diameter_nm,
                                                  args.particle_n, args.particle_k, args.substrate_n,
                                                  args.substrate_k, args.incidence_deg, args.domain_radius_um,
                                                  args.polarization, args.roughness_seed,
                                                  args.beam_waist_um, args.beam_order)):
        parser.error('FEM parameters cannot change with --source; omit --source to recompute scattering.')
    settings = imaging_settings(args)
    output = Path(args.output)
    source = Path(args.source) if args.source else output/'scattering'
    if not args.source:
        from .rough_run import run
        p, roughness = scattering_settings(args)
        run(p, roughness, source, angular_step=5, quadrature_order=7 if p.order >= 3 else 5,
            near_diagnostics=args.near_diagnostics)
    image_saved_sources(source, output, settings, args.archive_mode, args.pupil_cache)


if __name__ == '__main__':
    main()

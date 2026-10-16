"""Run the near/far/image pipeline or re-image saved complex volume sources."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from .parameters import Parameters
from .rough_surface import Roughness
from .imaging import calculate_images
from .numerical_settings import numerical_parameters
from .imaging_options import add_options, imaging_settings
from .near_planes import saved_planes


def image_saved_sources(source, output, settings, archive_mode='full', cache_directory=None):
    source, output = Path(source), Path(output)
    if archive_mode not in ('compact', 'full'):
        raise ValueError('Unknown archive mode')
    settings.validate()
    if settings.illumination_fit == 'camera':
        raise ValueError('Camera-matched roughness requires a whole-surface hybrid source, not a local FEM archive')
    scattering = json.loads((source/'report.json').read_text(encoding='utf-8'))
    if scattering.get('source_format_version') != 1:
        raise ValueError('This result lacks reusable volume sources. Run the pipeline again without --source.')
    with np.load(source/'paired_fields.npz', allow_pickle=False) as saved:
        data = {key: saved[key] for key in saved.files}
    p = Parameters(**scattering['parameters'])
    p.validate()
    if p.beam_profile == 'gaussian':
        from .finite_illumination import power_area
        if (settings.illumination_model != 'gaussian' or settings.illumination_definition != p.beam_definition or not np.allclose(
                (settings.illumination_length_um, settings.illumination_width_um,
                 settings.illumination_normalization_area_um2),
                (p.beam_length_um, p.beam_width_um, power_area(p)), rtol=1e-10, atol=0)):
            raise ValueError('Saved finite-beam scattering and imaging illumination differ')
    arrays, report = calculate_images(data, p, settings, retain_complex_fields=archive_mode == 'full',
                                      cache_directory=cache_directory)
    report.update(scattering_parameters=asdict(p), roughness=scattering['roughness'],
                  image_channel_figures=[1, 2, 3],
                  near_field_planes=[plane.metadata() for plane in saved_planes(data)],
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
    from .illumination_output import save_illumination
    save_illumination(output, p, settings, Roughness(**scattering['roughness']))
    print(json.dumps(report, indent=2), flush=True)
    return report


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/fullwave_imaging/default')
    parser.add_argument('--source', help='Reuse a rough_run directory with source_format_version=1')
    add_options(parser)
    parser.add_argument('--beam-waist-um', type=float)
    parser.add_argument('--beam-order', type=int)
    parser.add_argument('--beam-profile', choices=('plane', 'gaussian'))
    parser.add_argument('--beam-definition', choices=('wafer_gaussian', 'beam_cross_section',
                        'wafer_rectangle', 'measured_wafer'))
    parser.add_argument('--beam-measurement-path', default='')
    parser.add_argument('--beam-phase-assumption', choices=('flat', 'reject'), default='flat')
    parser.add_argument('--na', type=float, nargs=3, default=(.5, .5, .5))
    parser.add_argument('--tilts', type=float, nargs=3, default=(45, -45, 0))
    parser.add_argument('--magnification', type=float, default=20)
    parser.add_argument('--objective-focal-mm', type=float, default=20)
    parser.add_argument('--irradiance-w-um2', type=float)
    parser.add_argument('--illumination-power-w', type=float)
    parser.add_argument('--illumination-length-um', type=float)
    parser.add_argument('--illumination-width-um', type=float)
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
    parser.add_argument('--particle-mesh-um', type=float,
                        help='Particle-region mesh size in um; 0 or omitted uses the profile value')
    parser.add_argument('--quadrature-order', type=int, choices=(5, 7),
                        help='Volume-source quadrature order; default follows FEM order')
    parser.add_argument('--direct-solver', choices=('umfpack', 'sparsecholesky'))
    parser.add_argument('--mesh-scale', type=float)
    parser.add_argument('--rms-nm', type=float)
    parser.add_argument('--roughness-model', choices=('legacy', 'gaussian', 'abc',
                        'measured_height', 'measured_psd'), default='legacy')
    parser.add_argument('--measurement-path', default='')
    parser.add_argument('--detrend', choices=('none', 'mean', 'plane'), default='plane')
    parser.add_argument('--missing-policy', choices=('reject', 'nearest'), default='reject')
    parser.add_argument('--correlation-length-um', type=float, default=.05)
    parser.add_argument('--psd-exponent', type=float, default=3.)
    parser.add_argument('--rough-width-um', type=float, default=.2)
    parser.add_argument('--rough-grid-points', type=int, default=5)
    parser.add_argument('--gap-nm', type=float)
    for name in ('wavelength-nm', 'diameter-nm', 'particle-n', 'particle-k', 'substrate-n',
                 'substrate-k', 'incidence-deg', 'domain-radius-um'):
        parser.add_argument('--'+name, type=float)
    parser.add_argument('--polarization', choices=('s', 'p'))
    parser.add_argument('--roughness-seed', type=int)
    return parser


def scattering_settings(args, particle_enabled=True):
    scale = args.mesh_scale if args.mesh_scale is not None else 1
    local = args.local_mesh_scale if args.local_mesh_scale is not None else 1
    if not np.isfinite(local) or not .4 <= local <= 1:
        raise ValueError('local_mesh_scale must be in [0.4, 1]')
    particle_mesh = args.particle_mesh_um
    if particle_mesh is not None and (not np.isfinite(particle_mesh) or particle_mesh < 0):
        raise ValueError('particle_mesh_um must be nonnegative (0 = automatic)')
    if not particle_mesh:
        particle_mesh = .015*scale*local
    physical = {key: getattr(args, key) for key in ('wavelength_nm', 'diameter_nm', 'particle_n',
                'particle_k', 'substrate_n', 'substrate_k', 'incidence_deg', 'domain_radius_um', 'polarization',
                'beam_waist_um', 'beam_order')}
    p = numerical_parameters(args.accuracy or 'preview', scale, order=args.order, threads=args.threads,
        gap_nm=args.gap_nm if args.gap_nm is not None else 10,
        pml_power=args.pml_power, pml_mesh_scale=args.pml_mesh_scale,
        pml_strength=args.pml_strength, integration_bonus=args.integration_bonus,
        pml_order=args.pml_order, pml_inner_fraction=args.pml_inner_fraction,
        particle_mesh_um=particle_mesh,
        direct_solver=args.direct_solver, **physical)
    if args.beam_profile == 'gaussian':
        p = replace(p, beam_profile='gaussian', beam_length_um=args.illumination_length_um or 2.,
                    beam_width_um=args.illumination_width_um or 2.,
                    beam_definition=args.beam_definition or 'wafer_gaussian',
                    beam_order=args.beam_order or 31)
        if p.beam_definition == 'measured_wafer':
            from .beam_inputs import fingerprint, measured_envelope
            path = str(Path(args.beam_measurement_path).resolve())
            digest = fingerprint(path)
            _, _, _, info = measured_envelope(path, digest, args.beam_phase_assumption)
            p = replace(p, beam_measurement_path=path, beam_measurement_sha256=digest,
                        beam_phase_assumption=args.beam_phase_assumption,
                        beam_length_um=info['input_extent_y_um'], beam_width_um=info['input_extent_x_um'])
        p.validate()
    if args.illumination_fit == 'camera':
        from detector_model import DetectorParameters
        from .imaging_options import DETECTOR_OPTIONS
        from .field_coverage import fit_camera_beam
        detector = DetectorParameters(enabled=args.detector, **{
            key.replace('-', '_'): getattr(args, 'detector_'+key.replace('-', '_')) for key in DETECTOR_OPTIONS})
        p, _ = fit_camera_beam(p, detector, args.magnification, args.na, args.tilts,
                               args.minimum_field_illumination)
        args.illumination_length_um, args.illumination_width_um = p.beam_length_um, p.beam_width_um
        if args.illumination_power_w is None:
            from .finite_illumination import power_area
            args.illumination_power_w = (1. if args.irradiance_w_um2 is None else
                                         args.irradiance_w_um2*power_area(p))
    roughness = Roughness(rms_nm=args.rms_nm if args.rms_nm is not None else 2,
                         seed=args.roughness_seed if args.roughness_seed is not None else 7,
                         mesh_um=.025*scale*local, model=args.roughness_model,
                         correlation_length_um=args.correlation_length_um,
                         psd_exponent=args.psd_exponent,
                         measurement_path=str(Path(args.measurement_path).resolve()) if args.measurement_path else '',
                         detrend=args.detrend, missing_policy=args.missing_policy,
                         width_um=args.rough_width_um, grid_points=args.rough_grid_points)
    if roughness.model.startswith('measured_'):
        from .experimental_surface import measurement_report
        roughness = replace(roughness, rms_nm=measurement_report(roughness)['rms_nm'])
    roughness.validate(p, particle_enabled=particle_enabled)
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
                                                  args.particle_mesh_um, args.quadrature_order,
                                                  args.direct_solver, args.wavelength_nm, args.diameter_nm,
                                                  args.particle_n, args.particle_k, args.substrate_n,
                                                  args.substrate_k, args.incidence_deg, args.domain_radius_um,
                                                  args.polarization, args.roughness_seed,
                                                  args.beam_waist_um, args.beam_order, args.beam_profile, args.beam_definition)):
        parser.error('FEM parameters cannot change with --source; omit --source to recompute scattering.')
    if args.source and args.beam_measurement_path:
        parser.error('Measured illumination cannot change with saved scattering sources')
    output = Path(args.output)
    source = Path(args.source) if args.source else output/'scattering'
    if args.source:
        source_report = json.loads((source/'report.json').read_text(encoding='utf-8'))
        p = Parameters(**source_report['parameters'])
    else:
        p, roughness = scattering_settings(args)
    settings = imaging_settings(args, incidence_deg=p.incidence_deg, parameters=p)
    if not args.source and settings.illumination_fit == 'camera':
        if roughness.model not in ('gaussian', 'abc'):
            raise ValueError('Camera-matched whole-surface roughness requires Gaussian or ABC PSD')
        if args.near_diagnostics:
            raise ValueError('Dense near diagnostics are local-FEM-only; disable them for camera-matched hybrid imaging')
        from .hybrid_run import prepare_sources, export_images
        inputs = dict(parameters=asdict(p), roughness=asdict(roughness), imaging_settings=asdict(settings),
                      coverage=.99, variance_fraction=.9, taper_width_um=.025,
                      quadrature_order=args.quadrature_order or (7 if p.order >= 3 else 5))
        prepare_sources(inputs, output)
        export_images(output, output, asdict(settings))
        return
    if not args.source:
        from .rough_run import run
        run(p, roughness, source, angular_step=5,
            quadrature_order=args.quadrature_order or (7 if p.order >= 3 else 5),
            near_diagnostics=args.near_diagnostics)
    image_saved_sources(source, output, settings, args.archive_mode, args.pupil_cache)


if __name__ == '__main__':
    main()

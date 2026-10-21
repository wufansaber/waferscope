"""Bounded synthetic patterned-wafer Maxwell verification with explicit inputs."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
import shutil

import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.patterned_wafer import (LinePatch, build_mesh, surface_triangles,
                                               mirror_sources, sample_half_planes)


def case_inputs(args):
    from fullwave_prototype.imaging import ImagingSettings
    from detector_model import DetectorParameters
    from camera_profiles import profile_values
    factor = args.mesh_factor
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError('Positive finite mesh factor required')
    p = Parameters(wavelength_nm=266, diameter_nm=args.diameter_nm, gap_nm=10,
        domain_radius_um=1.5, order=args.order, air_mesh_um=.09*factor,
        substrate_mesh_um=.035*factor, particle_mesh_um=.015*factor,
        pml_thickness_um=.20, pml_mesh_scale=1.8, max_dofs=450000, threads=2,
        direct_solver='sparsecholesky', integration_bonus=args.integration_bonus,
        pml_order=args.pml_order, static_condensation=args.order > 1)
    pattern = LinePatch(pitch_um=args.pitch_nm/1000,
        trench_width_um=(args.pitch_nm-args.line_width_nm)/1000,
        trench_depth_um=args.trench_depth_nm/1000, line_length_um=args.line_length_um,
        half_width_um=args.half_width_um, air_height_um=args.air_height_um)
    detector = DetectorParameters(**dict(profile_values('gemi_area'), pixels_x=8, pixels_y=8,
                                         exposure_ms=args.exposure_ms))
    settings = ImagingSettings(numerical_apertures=(args.na,)*3, magnification=args.magnification,
        pupil_samples=64, padding=8, irradiance_w_um2=1/(4*np.cos(np.deg2rad(75))),
        illumination_power_w=1, illumination_length_um=2, illumination_width_um=2,
        illumination_incidence_deg=75, transmission=1e-6, detector=detector,
        camera_profile='gemi_area', gray_snr_pairs=1000)
    p.validate()
    pattern.validate(p)
    settings.validate()
    return p, pattern, settings


def worker(args):
    output = args.output.resolve()
    half_domain = args.half_domain
    import ngsolve as ng
    from fullwave_prototype.solver import solve
    from fullwave_prototype.volume_sources import polarization_sources
    from fullwave_prototype.near_planes import NearPlane, sample_planes, finish_planes
    from fullwave_prototype.near_plotting import plot_near_planes
    from fullwave_prototype.radiation import radiation_amplitude
    from fullwave_prototype.rough_run import decompose
    from fullwave_prototype.imaging import calculate_images
    from tools.plot_patterned_wafer_case import plot_case
    p, pattern, settings = case_inputs(args)
    output.mkdir(parents=True, exist_ok=True)
    ng.SetNumThreads(p.threads)
    print('Building finite etched lines with a shared particle subdomain...', flush=True)
    mesh = build_mesh(p, pattern, half_domain)
    xy = np.linspace(-.89*pattern.half_width_um, .89*pattern.half_width_um, 161)
    zz = np.linspace(-min(.14, pattern.substrate_depth_um*.9), pattern.air_height_um*.85, 121)
    center_z = (p.gap_nm+p.diameter_nm/2)/1000
    planes = (NearPlane('xy', xy, xy, center_z), NearPlane('xz', xy, zz, 0.), NearPlane('yz', xy, zz, 0.))
    axis = np.linspace(-1, 1, 51 if args.audit_only else 101)
    xx, yy = np.meshgrid(axis, axis)
    mask = xx*xx+yy*yy < 1-1e-12
    directions = np.column_stack((xx[mask], yy[mask], np.sqrt(1-xx[mask]**2-yy[mask]**2)))
    data, metrics = {}, {}
    for enabled, name in ((False, 'without'), (True, 'with')):
        print('Solving patterned wafer: '+name+' particle', flush=True)
        _, field, info = solve(p, mesh, particle_enabled=enabled)
        points, weights, volumes = polarization_sources(mesh, field, p, enabled, order=args.source_order)
        if half_domain:
            points, weights = mirror_sources(points, weights)
        print('Far field from {} volume quadrature sources...'.format(len(points)), flush=True)
        data['source_points_'+name], data['source_polarization_weights_'+name] = points, weights
        data['amplitude_'+name] = radiation_amplitude(points, weights, directions, p)
        check_points, check_weights, _ = polarization_sources(mesh, field, p, enabled, order=args.source_order+2)
        if half_domain:
            check_points, check_weights = mirror_sources(check_points, check_weights)
        selection = np.linspace(0, len(directions)-1, 41).astype(int)
        check = radiation_amplitude(check_points, check_weights, directions[selection], p)
        info['radiation_quadrature_relative_l2'] = float(np.linalg.norm(check-data['amplitude_'+name][selection])/np.linalg.norm(check))
        info['contrast_volumes_um3'] = volumes
        info['source_quadrature_order'] = args.source_order
        info['source_check_quadrature_order'] = args.source_order+2
        info['source_quadrature_points_full_domain'] = len(points)
        info['mesh_volume_fraction'] = .5 if half_domain else 1.
        metrics[name] = info
        (sample_half_planes if half_domain else sample_planes)(mesh, field, p, planes, name, data)
        del field, check_points, check_weights
    finish_planes(data, planes)
    data['surface_triangles_um'] = surface_triangles(pattern)
    for name, values in decompose(data['amplitude_with'], data['amplitude_without']).items():
        disk = np.full(xx.shape, np.nan)
        disk[mask] = values
        data['disk_'+name] = disk
    data['disk_axis'] = axis
    report = dict(parameters=asdict(p), pattern=asdict(pattern), metrics=metrics,
                  half_domain_pmc_symmetry=half_domain,
                  near_field_planes=[plane.metadata() for plane in planes],
                  status='synthetic finite-patch FEM preview; not SEM reconstruction or mesh-converged prediction',
                  scope='Two finite trenches on homogeneous half-space; sphere above central ridge. '
                        'Shared mesh and phase for sample and reference; flat specular excluded. '
                        'No random haze, coating, periodic boundary or infinite-grating approximation. '
                        'Plane-wave illumination; 1 W / 2 x 2 um wafer area is irradiance normalization, not a finite beam.')
    plot_near_planes(output, data, report)
    def provider(directions, state):
        return radiation_amplitude(data['source_points_'+state], data['source_polarization_weights_'+state], directions, p)

    arrays, imaging = calculate_images(data, p, settings, retain_complex_fields=False, amplitude_provider=provider)
    imaging.update(scattering_parameters=asdict(p), pattern=asdict(pattern), status=report['status'])
    np.savez_compressed(output/'scattering_fields.npz', **data)
    np.savez_compressed(output/'image_fields.npz', **arrays)
    (output/'scattering_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    (output/'imaging_report.json').write_text(json.dumps(imaging, indent=2), encoding='utf-8')
    plot_case(output, data, arrays, imaging)
    if args.export_to:
        args.export_to.mkdir(parents=True, exist_ok=True)
        for path in output.iterdir():
            if path.suffix in ('.png', '.json'):
                shutil.copy2(path, args.export_to/path.name)
    print(json.dumps(dict(status='completed', output=str(output), dofs=metrics['with']['dofs'],
        powers=[c['with']['image_power_w'] for c in imaging['channels']],
        snr=[c['detector']['gray_snr'] for c in imaging['channels']]), indent=2), flush=True)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('build/patterned_wafer_20261019/preview'))
    parser.add_argument('--mesh-factor', type=float, default=1.)
    parser.add_argument('--order', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--pml-order', type=int, choices=(0, 1, 2, 3), default=0)
    parser.add_argument('--integration-bonus', type=int, choices=range(13), default=10)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--half-domain', action='store_true')
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--diameter-nm', type=float, default=100)
    parser.add_argument('--pitch-nm', type=float, default=1000)
    parser.add_argument('--line-width-nm', type=float, default=500)
    parser.add_argument('--trench-depth-nm', type=float, default=100)
    parser.add_argument('--line-length-um', type=float, default=1.4)
    parser.add_argument('--half-width-um', type=float, default=.9)
    parser.add_argument('--air-height-um', type=float, default=.30)
    parser.add_argument('--magnification', type=float, default=20)
    parser.add_argument('--na', type=float, default=.5)
    parser.add_argument('--exposure-ms', type=float, default=.01)
    parser.add_argument('--source-order', type=int, choices=(5, 7, 9, 11, 13), default=11)
    parser.add_argument('--export-to', type=Path)
    return parser


def main():
    args = build_parser().parse_args()
    case_inputs(args)
    if args.worker:
        worker(args)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        command = [sys.executable, '-m', 'tools.patterned_wafer_case', '--worker']+sys.argv[1:]
        report = run_guarded(command,
            args.output.parent/(args.output.name+'_resources.json'), maximum=9., reserve=3., timeout=900.)
        if report['status'] != 'completed' or report['returncode']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()

"""Paired coherent Maxwell solves on one finite rough-surface realization."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np

from .rough_surface import Roughness, build_rough_mesh, height_mesh
from .solver import solve, evaluate_fields
from .volume_sources import polarization_sources
from .far_field import angular_grid, radiation_amplitude
from .run import hemisphere_integral
from .numerical_settings import numerical_parameters


def decompose(amplitude_with, amplitude_without):
    delta = amplitude_with-amplitude_without
    return dict(total=np.sum(abs(amplitude_with)**2, axis=-1),
                background=np.sum(abs(amplitude_without)**2, axis=-1),
                change_field=np.sum(abs(delta)**2, axis=-1),
                interference=2*np.real(np.sum(delta*np.conj(amplitude_without), axis=-1)))


def run(p, roughness, output, angular_step=10, quadrature_order=5, near_diagnostics=False):
    import ngsolve as ng
    p.validate()
    roughness.validate(p)
    if angular_step not in (2, 5, 10, 15):
        raise ValueError('angular_step must be 2, 5, 10 or 15 degrees')
    if quadrature_order not in range(3, 10):
        raise ValueError('quadrature_order must be 3..9')
    start = time.perf_counter()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    ng.SetNumThreads(p.threads)
    print('Generating one shared rough-surface mesh for the paired solves...', flush=True)
    mesh = build_rough_mesh(p, roughness)
    theta, phi, directions = angular_grid(angular_step)
    lateral = np.linspace(-.6*p.domain_radius_um, .6*p.domain_radius_um, 101)
    height = np.linspace(-.299*p.domain_radius_um, .75*p.domain_radius_um, 101)
    yy, zz = np.meshgrid(lateral, height)
    near_points = np.column_stack((np.zeros(yy.size), yy.ravel(), zz.ravel()))
    data = dict(theta=theta, phi=phi, near_y_um=lateral, near_z_um=height)
    metrics = {}
    for enabled, name in ((False, 'without'), (True, 'with')):
        print('Solving {} particle...'.format(name), flush=True)
        _, field, info = solve(p, mesh=mesh, particle_enabled=enabled)
        points, sources, volumes = polarization_sources(mesh, field, p, enabled, quadrature_order)
        amplitude = radiation_amplitude(points, sources, directions.reshape(-1, 3), p)
        cp, cs, _ = polarization_sources(mesh, field, p, enabled, quadrature_order+2)
        sampled = directions.reshape(-1, 3)[::12]
        check = radiation_amplitude(cp, cs, sampled, p)
        info['radiation_quadrature_amplitude_relative_l2'] = float(
            np.linalg.norm(check-amplitude[::12])/max(np.linalg.norm(check), 1e-300))
        info['contrast_volumes_um3'] = volumes
        info['source_quadrature_points'] = len(points)
        info['source_quadrature_order'] = quadrature_order
        e, h, _ = evaluate_fields(mesh, field, near_points, p)
        data['amplitude_'+name] = amplitude.reshape(directions.shape[:2]+(2,))
        data['near_e_'+name] = e.reshape(yy.shape+(3,))
        data['near_hz0_'+name] = h.reshape(yy.shape+(3,))
        # Retain phase-resolved volume sources so objectives can change without another FEM solve.
        data['source_points_'+name] = points
        data['source_polarization_weights_'+name] = sources
        metrics[name] = info
        if enabled and near_diagnostics:
            from .near_audit import export_diagnostics
            local_scale = (p.particle_mesh_um/.015)/(p.air_mesh_um/.08)
            export_diagnostics(mesh, field, p, roughness, info, output/'near_audit', local_scale)
        del field, points, sources, cp, cs
    components = decompose(data['amplitude_with'], data['amplitude_without'])
    for name, value in components.items():
        data['dsc_'+name] = value
    data['amplitude_change'] = data['amplitude_with']-data['amplitude_without']
    data['near_e_change'] = data['near_e_with']-data['near_e_without']
    vertices, triangles = height_mesh(roughness)
    data['surface_vertices_um'], data['surface_triangles_um'] = vertices, triangles
    if not all(np.isfinite(value).all() for value in data.values()):
        raise RuntimeError('Non-finite paired result; export rejected')
    integrals = {name: hemisphere_integral(value, theta) for name, value in components.items()}
    integrals['intensity_difference'] = integrals['total']-integrals['background']
    metrics['hemisphere_integrals_um2'] = integrals
    metrics['decomposition_absolute_error_um2_sr'] = float(np.max(abs(
        components['total']-components['background']-components['change_field']-components['interference'])))
    metrics['elapsed_seconds'] = time.perf_counter()-start
    payload = dict(parameters=asdict(p), roughness=asdict(roughness), metrics=metrics,
                   source_format_version=1,
                   angular_step_deg=angular_step, status='experimental; roughness not physically converged',
                   scope='Finite triangulated rough patch on half-space, sphere entirely in air; same mesh in both solves.',
                   units='E/E_inc; H*Z0/E_inc; amplitudes um; DSC um^2/sr; exp(-i omega t).',
                   surface_mean_height_nm=float(triangles[..., 2].mean()*1000),
                   surface_peak_height_nm=float(vertices[..., 2].max()*1000),
                   surface_min_height_nm=float(vertices[..., 2].min()*1000),
                   note='Change-field power is NOT camera background subtraction. Total-background includes signed interference. No finite beam or ensemble average.')
    np.savez_compressed(output/'paired_fields.npz', **data)
    (output/'report.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    from .rough_plotting import plot_pair
    plot_pair(output, data, payload)
    print(json.dumps(payload, indent=2), flush=True)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/fem_roughness/default')
    parser.add_argument('--rms-nm', type=float, default=2)
    parser.add_argument('--width-um', type=float, default=.2)
    parser.add_argument('--grid-points', type=int, default=5)
    parser.add_argument('--smoothing-cells', type=float, default=.6)
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--gap-nm', type=float, default=10)
    parser.add_argument('--order', type=int)
    parser.add_argument('--accuracy', choices=('preview', 'refined'), default='preview')
    parser.add_argument('--pml-power', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--pml-mesh-scale', type=float, default=1.)
    parser.add_argument('--mesh-scale', type=float, default=1)
    parser.add_argument('--angular-step', type=int, choices=(2, 5, 10, 15), default=10)
    parser.add_argument('--pml-thickness', type=float, default=.16)
    parser.add_argument('--polarization', choices=('s', 'p'), default='p')
    args = parser.parse_args()
    p = numerical_parameters(args.accuracy, args.mesh_scale, order=args.order,
                   gap_nm=args.gap_nm, polarization=args.polarization,
                   pml_thickness_um=args.pml_thickness, pml_power=args.pml_power,
                   pml_mesh_scale=args.pml_mesh_scale)
    r = Roughness(width_um=args.width_um, rms_nm=args.rms_nm, grid_points=args.grid_points,
                  smoothing_cells=args.smoothing_cells, seed=args.seed, mesh_um=.025*args.mesh_scale)
    run(p, r, args.output, angular_step=args.angular_step, quadrature_order=7 if p.order >= 3 else 5)


if __name__ == '__main__':
    main()

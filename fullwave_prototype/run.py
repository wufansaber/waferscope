"""CLI benchmark: native FEM fields, half-space radiation and optional SCATMECH comparison."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from .parameters import Parameters
from .solver import solve, evaluate_fields
from .far_field import sphere_quadrature, radiation_amplitude, angular_grid


def run(p, output, reference_python=None):
    if p.defect_shape != 'sphere':
        raise ValueError('This benchmark is sphere-only; use imaging_run for finite defect patterns')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    mesh, field, metrics = solve(p)
    points, weights = sphere_quadrature(p)
    electric, _, _ = evaluate_fields(mesh, field, points, p)
    theta, phi, directions = angular_grid(5)
    polarization = (p.eps_particle-1)*electric*weights[:, None]
    amplitudes = radiation_amplitude(points, polarization, directions.reshape(-1, 3), p)
    dsc = np.sum(abs(amplitudes)**2, axis=-1).reshape(directions.shape[:2])
    metrics['hemisphere_dsc_um2'] = hemisphere_integral(dsc, theta)
    check_points, check_weights = sphere_quadrature(p, radial_order=14)
    check_e, _, _ = evaluate_fields(mesh, field, check_points, p)
    check_amplitude = radiation_amplitude(check_points, (p.eps_particle-1)*check_e*check_weights[:, None],
                                         directions.reshape(-1, 3)[::24], p)
    check_dsc = np.sum(abs(check_amplitude)**2, axis=-1)
    metrics['radiation_quadrature_relative_l2'] = float(np.linalg.norm(check_dsc-dsc.ravel()[::24])
                                                       /max(np.linalg.norm(check_dsc), 1e-300))
    # Sample between interfaces, not at discontinuous dielectric boundaries.
    lateral = np.linspace(-.6*p.domain_radius_um, .6*p.domain_radius_um, 101)
    height = np.linspace(-.299*p.domain_radius_um, .75*p.domain_radius_um, 101)
    yy, zz = np.meshgrid(lateral, height)
    near_points = np.column_stack((np.zeros(yy.size), yy.ravel(), zz.ravel()))
    e, h, esc = evaluate_fields(mesh, field, near_points, p)
    if not all(np.isfinite(a).all() for a in (amplitudes, e, h, esc)):
        raise RuntimeError('Non-finite FEM fields; result rejected')
    np.savez_compressed(output/'fields.npz', theta=theta, phi=phi, dsc_um2_sr=dsc,
                        amplitude_sp_um=amplitudes.reshape(directions.shape[:2]+(2,)),
                        near_y_um=lateral, near_z_um=height,
                        near_e=e.reshape(yy.shape+(3,)), near_hz0=h.reshape(yy.shape+(3,)),
                        near_scattered_e=esc.reshape(yy.shape+(3,)))
    payload = dict(parameters=asdict(p), metrics=metrics,
                   status='experimental, not production validated',
                   scope='Sphere in air above smooth half-space only; no films or roughness.',
                   units='Lengths um; E/E_inc; H*Z0/E_inc; DSC um^2/sr. exp(-i omega t).',
                   near_plane='YZ (incidence plane), x=0; interface samples are not one-sided limits.')
    config = output/'report.json'
    config.write_text(json.dumps(payload, indent=2), encoding='utf-8')
    reference = None
    if reference_python:
        subprocess.run([reference_python, '-m', 'fullwave_prototype.reference',
                        '--config', str(config.resolve()), '--output', str((output/'scatmech.npz').resolve())],
                       cwd=str(Path(__file__).resolve().parents[1]), check=True)
        with np.load(output/'scatmech.npz') as data:
            reference = data['dsc_um2_sr']
        ref_total = hemisphere_integral(reference, theta)
        metrics['reference_hemisphere_dsc_um2'] = ref_total
        metrics['hemisphere_relative_error'] = abs(metrics['hemisphere_dsc_um2']/ref_total-1)
        w = np.sin(theta)[:, None]
        metrics['angular_weighted_relative_l2'] = float(np.sqrt(
            np.sum(w*(dsc-reference)**2)/np.sum(w*reference**2)))
    metrics['elapsed_seconds'] = time.perf_counter()-start
    config.write_text(json.dumps(payload, indent=2), encoding='utf-8')
    from .plotting import plot_result
    plot_result(output, p, metrics, reference)
    print(json.dumps(payload, indent=2), flush=True)
    return payload


def hemisphere_integral(dsc, theta):
    weighted = np.mean(dsc, axis=1)*np.sin(theta)
    return float(np.pi*np.sum((weighted[:-1]+weighted[1:])*np.diff(theta)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/fem_validation/default')
    parser.add_argument('--reference-python')
    parser.add_argument('--order', type=int, default=2)
    parser.add_argument('--gap-nm', type=float, default=0)
    parser.add_argument('--mesh-scale', type=float, default=1)
    parser.add_argument('--radius', type=float, default=.22)
    parser.add_argument('--pml-thickness', type=float, default=.16)
    parser.add_argument('--polarization', choices=['s', 'p'], default='p')
    parser.add_argument('--substrate-n', type=float, default=5.57)
    parser.add_argument('--substrate-k', type=float, default=.387)
    parser.add_argument('--wavelength-nm', type=float, default=355)
    parser.add_argument('--diameter-nm', type=float, default=100)
    parser.add_argument('--incidence-deg', type=float, default=75)
    args = parser.parse_args()
    p = Parameters(order=args.order, gap_nm=args.gap_nm,
                   air_mesh_um=.08*args.mesh_scale, substrate_mesh_um=.04*args.mesh_scale,
                   particle_mesh_um=.015*args.mesh_scale, domain_radius_um=args.radius,
                   pml_thickness_um=args.pml_thickness, polarization=args.polarization,
                   substrate_n=args.substrate_n, substrate_k=args.substrate_k,
                   wavelength_nm=args.wavelength_nm, diameter_nm=args.diameter_nm,
                   incidence_deg=args.incidence_deg)
    run(p, args.output, args.reference_python)


if __name__ == '__main__':
    main()

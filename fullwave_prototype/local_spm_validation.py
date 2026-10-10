"""Independent haze-only low-height and FEM-order limits on one stored surface."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

import numpy as np

from .background_sampling import electric_background
from .far_field import radiation_amplitude
from .hybrid_radiation import triangle_source, triangle_amplitude
from .hybrid_run import load_sources
from .parameters import Parameters
from .rough_surface import Roughness, height_mesh
from .shared_surface import local_crop, load_mother


def probe_directions():
    theta = np.deg2rad([10, 30, 50, 70, 85])
    phi = np.deg2rad(np.arange(0, 360, 30))
    t, f = np.meshgrid(theta, phi, indexing='ij')
    return np.stack((np.sin(t)*np.cos(f), np.sin(t)*np.sin(f), np.cos(t)), axis=-1).reshape(-1, 3)


def comparison(fem, spm):
    reference = max(np.linalg.norm(spm), 1e-300)
    error = float(np.linalg.norm(fem-spm)/reference)
    overlap = np.vdot(spm, fem)
    gain = overlap/max(float(np.vdot(spm, spm).real), 1e-300)
    return dict(complex_relative_l2=error, best_scalar_gain=[float(gain.real), float(gain.imag)],
        scalar_corrected_relative_l2=float(np.linalg.norm(fem-gain*spm)/reference),
        sampled_intensity_ratio=float(np.sum(abs(fem)**2)/max(np.sum(abs(spm)**2), 1e-300)),
        scope='Selected direction samples, NOT hemisphere-integrated energy; no fitted gain is applied to results')


def run_case(p, rough, axis, height, scale, order, output, power_checks=False):
    import ngsolve as ng
    from .rough_geometry import build_rough_mesh
    from .solver import solve, background_coefficient
    from .volume_sources import polarization_sources
    output.mkdir(parents=True, exist_ok=False)
    local = local_crop(output, axis, height*scale, rough)
    p = replace(p, order=order, pml_order=1, static_condensation=order > 1)
    p.validate()
    local.validate(p, particle_enabled=False)
    started = time.perf_counter()
    print('Haze-only validation: height scale {}, FEM order {}...'.format(scale, order), flush=True)
    mesh = build_rough_mesh(p, local, particle_enabled=False)
    _, field, solver = solve(p, mesh=mesh, particle_enabled=False)
    points, sources, volumes = polarization_sources(mesh, field, p, False, 5, vectorized_background=True)
    samples = points[::max(1, len(points)//32)][:32]
    locations = mesh(samples[:, 0], samples[:, 1], samples[:, 2])
    original = np.asarray(background_coefficient(p)(locations)).reshape(-1, 3)
    batch = electric_background(samples, p)
    electric_error = float(np.linalg.norm(batch-original)/max(np.linalg.norm(original), 1e-300))
    if electric_error > 1e-11:
        raise RuntimeError('Batched background differs from the native source reference')
    q = probe_directions()
    fem = radiation_amplitude(points, sources, q, p)
    fine_points, fine_sources, _ = polarization_sources(mesh, field, p, False, 7, vectorized_background=True)
    quadrature_error = float(np.linalg.norm(fem-radiation_amplitude(fine_points, fine_sources, q, p))/
                             max(np.linalg.norm(fem), 1e-300))
    _, triangles = height_mesh(local)
    sp, sw = triangle_source(triangles, p, order=8)
    spm = triangle_amplitude(sp, sw, q, p)
    result = dict(height_scale=scale, input_stationary_rms_nm=rough.rms_nm*scale,
        fem_order=order, parameters=asdict(p), roughness=asdict(local),
        solver=solver, contrast_volumes_um3=volumes, comparison=comparison(fem, spm),
        source_quadrature_relative_l2_change=quadrature_error,
        batched_vs_native_background_relative_l2=electric_error,
        elapsed_seconds=time.perf_counter()-started,
        scope='Haze only; NO particle, no near/image export, same mother realization scaled in height. '
              'Each physical height requires its own geometry/mesh. Material constants stay fixed.')
    if power_checks:
        from .energy_balance import check_balance
        from .finite_illumination import power_area
        from .surface_observables import integrals
        from .radiation import radiation_amplitude as blocked_amplitude

        class VolumeRadiation:
            def amplitude(self, directions):
                shape = np.asarray(directions).shape[:-1]
                return blocked_amplitude(points, sources, np.asarray(directions).reshape(-1, 3), p).reshape(shape+(2,))

        class TriangleRadiation:
            def amplitude(self, directions):
                shape = np.asarray(directions).shape[:-1]
                return triangle_amplitude(sp, sw, np.asarray(directions).reshape(-1, 3), p).reshape(shape+(2,))

        irradiance = 1/power_area(p)
        a, tilts = (.5, .5, .5), (45, -45, 0)
        fem_power = integrals(VolumeRadiation(), a, tilts, 18)*irradiance
        fem_coarse = integrals(VolumeRadiation(), a, tilts, 12)*irradiance
        spm_power = integrals(TriangleRadiation(), a, tilts, 18)*irradiance
        result['absolute_power_check'] = dict(laser_power_w=1., hemisphere_and_channels_fem_w=fem_power.tolist(),
            hemisphere_and_channels_spm_w=spm_power.tolist(),
            fem_to_spm_ratio=(fem_power/np.maximum(spm_power, 1e-300)).tolist(),
            angular_orders=[12, 18], fem_angular_relative_changes=(abs(fem_power-fem_coarse)/
                np.maximum(fem_power, 1e-300)).tolist(), scope='Scattered optical powers, no specular term. Fixed test apertures.')
        result['energy_balance'] = check_balance(mesh, field, p, irradiance)
    result['elapsed_seconds'] = time.perf_counter()-started
    np.savez_compressed(output/'fields.npz', directions=q, fem=fem, spm=spm, triangles_um=triangles)
    (output/'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(dict(scale=scale, order=order, **result['comparison']), indent=2), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--orders', type=int, nargs='+', choices=(1, 2), default=[1, 2])
    parser.add_argument('--direct-solver', choices=('umfpack', 'sparsecholesky'))
    parser.add_argument('--power-checks', action='store_true')
    parser.add_argument('--height-scales', type=float, nargs='+', default=[1., .25, .1])
    args = parser.parse_args()
    _, metadata = load_sources(args.source)
    axis, height, mother = load_mother(args.source)
    p = Parameters(**metadata['parameters'])
    if args.direct_solver:
        p = replace(p, direct_solver=args.direct_solver)
    # Retain the original local physical taper, not a newly randomized local patch.
    rough = Roughness(**metadata['stationary_roughness'])
    rough = replace(rough, taper_fraction=metadata['roughness']['taper_fraction'])
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    report = dict(model='shared_surface_low_height_fem_spm_validation', status='running',
        source=str(Path(args.source).resolve()), mother_sha256=mother['sha256'], cases=[],
        input_stationary_rms_nm=rough.rms_nm,
        full_model_converged=False, scope='Independent deterministic limit test, not ensemble or full-model certification')
    path = output/'validation_report.json'
    if any(not np.isfinite(s) or not 0 < s <= 1 for s in args.height_scales):
        raise ValueError('Height scales must be in (0,1]')
    cases = [(s, order) for order in args.orders for s in (args.height_scales if order == 1 else [min(args.height_scales)])]
    for scale, order in cases:
        name = 'scale_{}_order_{}'.format(str(scale).replace('.', 'p'), order)
        case = run_case(p, rough, axis, height, scale, order, output/name, args.power_checks)
        report['cases'].append(case)
        path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    report['status'] = 'complete_with_physical_limits'
    path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    from .validation_plotting import plot_local
    plot_local(output, report)


if __name__ == '__main__':
    main()

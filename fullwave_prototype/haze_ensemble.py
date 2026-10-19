"""Independent finite-surface FEM haze realizations and intensity statistics."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import subprocess
import sys
import hashlib
import time

import numpy as np

from .far_field import angular_grid, radiation_amplitude
from .imaging_run import build_parser, scattering_settings
from .parameters import Parameters
from .rough_surface import (Roughness, build_rough_mesh, height_mesh, surface_statistics,
                            taper_effective_area)
from .run import hemisphere_integral
from .solver import solve
from .volume_sources import polarization_sources


def simulate_sample(parameters, roughness, output, angular_step=5, quadrature_order=5):
    """One haze-only solve; do not retain or average its complex amplitude."""
    parameters.validate()
    start = time.perf_counter()
    if roughness.model not in ('gaussian', 'abc', 'measured_psd', 'measured_height') or roughness.rms_nm <= 0:
        raise ValueError('Haze requires nonzero statistical or measured roughness')
    if parameters.beam_waist_um != 0:
        raise ValueError('Use plane or Gaussian footprint mode for haze, not legacy beam waist')
    roughness.validate(parameters, particle_enabled=False)
    mesh = build_rough_mesh(parameters, roughness, particle_enabled=False)
    _, field, solver_info = solve(parameters, mesh=mesh, particle_enabled=False)
    points, sources, volumes = polarization_sources(
        mesh, field, parameters, particle_enabled=False, order=quadrature_order)
    theta, phi, directions = angular_grid(angular_step)
    amplitude = radiation_amplitude(points, sources, directions.reshape(-1, 3), parameters)
    dsc = np.sum(abs(amplitude.reshape(directions.shape[:2]+(2,)))**2, axis=-1)
    if not np.isfinite(dsc).all() or np.any(dsc < 0):
        raise RuntimeError('Non-finite haze DSC')
    vertices, triangles = height_mesh(roughness)
    from .coupled_surface import realization_signature
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    report = dict(seed=roughness.seed, parameters=asdict(parameters), roughness=asdict(roughness),
                  surface_statistics=surface_statistics(roughness, vertices, triangles),
                  solver=solver_info, source_volumes_um3=volumes,
                  hemisphere_dsc_um2=hemisphere_integral(dsc, theta),
                  surface_sha256=hashlib.sha256(vertices.tobytes()).hexdigest(),
                  realization_sha256=realization_signature(roughness),
                  elapsed_seconds=time.perf_counter()-start,
                  scope='One finite tapered surface; haze only; no specular reflection or particle.')
    from .illumination_output import illumination_metadata
    report['illumination'] = illumination_metadata(parameters, roughness)
    np.savez_compressed(output/'haze_dsc.npz', theta=theta, phi=phi, dsc=dsc)
    np.savez_compressed(output/'haze_sources.npz', points=points, weights=sources)
    (output/'sample_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def summarize(samples, reference, theta, effective_area_um2, comparison_area_um2,
              roughness, finite_reference=None, finite_beam=False):
    """Average independent DSC intensities on a common solid-angle grid."""
    values = np.stack(samples)
    if values.ndim != 3 or values.shape[1] != len(theta) or not np.isfinite(values).all():
        raise ValueError('Sample maps must be finite and share one angular grid')
    if np.any(values < 0):
        raise ValueError('Sample DSC must be nonnegative')
    n = len(values)
    if effective_area_um2 <= 0 or comparison_area_um2 <= 0:
        raise ValueError('Comparison and effective roughness areas must be positive')
    scale = (comparison_area_um2/effective_area_um2
             if roughness.model != 'measured_height' and not finite_beam else 1.)
    mean = values.mean(axis=0)*scale
    standard_error = (values.std(axis=0, ddof=1)*scale/np.sqrt(n)
                      if n > 1 else np.full(mean.shape, np.nan))
    hemispheres = np.array([hemisphere_integral(value, theta) for value in values])*scale
    from .rough_statistics import nyquist_variance_fraction
    spacing = roughness.width_um/(roughness.grid_points-1)
    metrics = dict(samples=n, comparison_area_um2=comparison_area_um2,
                   patch_geometric_area_um2=roughness.width_um**2,
                   patch_effective_area_um2=effective_area_um2,
                   patch_taper_area_fraction=effective_area_um2/roughness.width_um**2,
                   fem_mean_hemisphere_dsc_um2=float(hemispheres.mean()),
                   fem_hemisphere_sample_std_um2=(float(hemispheres.std(ddof=1)) if n > 1 else None),
                   fem_hemisphere_standard_error_um2=(float(hemispheres.std(ddof=1)/np.sqrt(n))
                                                      if n > 1 else None),
                   patch_width_to_psd_length=roughness.width_um/roughness.correlation_length_um,
                   psd_model=roughness.model,
                   psd_length_um=roughness.correlation_length_um,
                   psd_exponent=roughness.psd_exponent if roughness.model == 'abc' else None,
                   nyquist_variance_fraction=(nyquist_variance_fraction(
                       roughness.model, spacing, roughness.correlation_length_um,
                       roughness.psd_exponent) if roughness.model in ('gaussian', 'abc') else None),
                   edge_taper_width_um=roughness.width_um*roughness.taper_fraction,
                   physical_convergence='unverified')
    metrics['area_extrapolation_applied'] = not finite_beam and roughness.model != 'measured_height'
    if roughness.model == 'gaussian':
        metrics['patch_width_to_correlation_length'] = metrics['patch_width_to_psd_length']
    if reference is not None:
        if reference.shape != mean.shape or not np.isfinite(reference).all() or np.any(reference < 0):
            raise ValueError('NIST reference must be finite, nonnegative and on the FEM grid')
        weight = np.sin(theta)[:, None]
        denominator = np.sqrt(np.sum(weight*mean*mean)*np.sum(weight*reference*reference))
        reference_total = hemisphere_integral(reference, theta)
        metrics.update(nist_hemisphere_dsc_um2=reference_total,
                       fem_to_nist_hemisphere_ratio=(float(hemispheres.mean()/reference_total)
                                                     if reference_total > 0 else None),
                       angular_shape_cosine=float(np.sum(weight*mean*reference)/denominator)
                            if denominator > 0 else None)
    if finite_reference is not None:
        if (finite_reference.shape != mean.shape or not np.isfinite(finite_reference).all() or
                np.any(finite_reference < 0)):
            raise ValueError('Finite-patch reference must be finite, nonnegative and on the FEM grid')
        weight = np.sin(theta)[:, None]
        denominator = np.sqrt(np.sum(weight*mean*mean)*
                              np.sum(weight*finite_reference*finite_reference))
        finite_total = hemisphere_integral(finite_reference, theta)
        metrics.update(finite_patch_hemisphere_dsc_um2=finite_total,
                       fem_to_finite_patch_hemisphere_ratio=(float(hemispheres.mean()/finite_total)
                                                            if finite_total > 0 else None),
                       finite_patch_angular_shape_cosine=float(
                           np.sum(weight*mean*finite_reference)/denominator)
                           if denominator > 0 else None)
        if reference is not None:
            denominator = np.sqrt(np.sum(weight*reference*reference)*
                                  np.sum(weight*finite_reference*finite_reference))
            metrics.update(finite_patch_to_nist_hemisphere_ratio=(float(
                               finite_total/metrics['nist_hemisphere_dsc_um2'])
                               if metrics['nist_hemisphere_dsc_um2'] > 0 else None),
                           finite_patch_to_nist_angular_shape_cosine=float(
                               np.sum(weight*reference*finite_reference)/denominator)
                               if denominator > 0 else None)
    return mean, standard_error, metrics


def run_ensemble(parameters, roughness, output, samples, comparison_area_um2=.04,
                 reference_path=None, angular_step=5, quadrature_order=5,
                 apertures=(.5, .5, .5), tilts=(45, -45, 0), target_relative_ci=.05):
    if not 1 <= samples <= 128:
        raise ValueError('Sample count must be 1..128')
    if roughness.seed+samples-1 > 1000000:
        raise ValueError('Last roughness seed exceeds 1000000')
    if comparison_area_um2 <= 0 or not np.isfinite(comparison_area_um2):
        raise ValueError('Comparison area must be positive and finite')
    parameters.validate()
    roughness.validate(parameters, particle_enabled=False)
    if roughness.model not in ('gaussian', 'abc', 'measured_psd') or roughness.rms_nm <= 0 or parameters.beam_waist_um != 0:
        raise ValueError('Use statistical roughness and plane/Gaussian footprint mode; measured height is deterministic')
    theta, phi, _ = angular_grid(angular_step)
    from .haze_statistics import statistical_report
    statistical_report([np.zeros((len(theta), len(phi)))], theta, phi, 1.,
                       apertures, tilts, target_relative_ci)
    reference = None
    finite_reference = None
    finite_beam = parameters.beam_profile == 'gaussian'
    if finite_beam and reference_path is not None:
        raise ValueError('A plane-wave NIST PSD reference is not a finite coherent-beam reference')
    if reference_path is not None:
        with np.load(reference_path, allow_pickle=False) as saved:
            if not np.array_equal(theta, saved['theta']) or not np.array_equal(phi, saved['phi']):
                raise ValueError('NIST reference angular grid differs from FEM grid')
            if 'psd_model' in saved.files and str(saved['psd_model']) != roughness.model:
                raise ValueError('NIST reference PSD model differs from FEM roughness model')
            if ('comparison_area_um2' in saved.files and
                    not np.isclose(float(saved['comparison_area_um2']), comparison_area_um2)):
                raise ValueError('NIST reference comparison area differs from FEM area')
            for key, actual in (('parameters_json', asdict(parameters)), ('roughness_json', asdict(roughness))):
                if key in saved.files and json.loads(str(saved[key])) != actual:
                    raise ValueError('NIST reference input parameters differ from FEM: '+key)
            if 'measurement_sha256' in saved.files and roughness.model == 'measured_psd':
                from .experimental_surface import input_fingerprint
                if str(saved['measurement_sha256']) != input_fingerprint(roughness.measurement_path):
                    raise ValueError('Measurement changed after the NIST reference was computed')
            reference = saved['dsc']
            if 'finite_patch_dsc' in saved.files:
                finite_reference = saved['finite_patch_dsc']
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    input_path = output/'haze_ensemble_inputs.json'
    measurement_sha = None
    if roughness.model == 'measured_psd':
        from .experimental_surface import input_fingerprint
        measurement_sha = input_fingerprint(roughness.measurement_path)
    input_path.write_text(json.dumps(dict(parameters=asdict(parameters), roughness=asdict(roughness),
        samples=samples, comparison_area_um2=comparison_area_um2,
        measurement_sha256=measurement_sha,
        angular_step_deg=angular_step, quadrature_order=quadrature_order), indent=2), encoding='utf-8')
    maps, sample_reports = [], []
    for offset in range(samples):
        seed = roughness.seed+offset
        folder = output/'samples'/('seed_'+str(seed))
        print('Haze FEM sample {}/{} (seed {})'.format(offset+1, samples, seed), flush=True)
        subprocess.run([sys.executable, '-m', 'fullwave_prototype.haze_ensemble', 'sample',
                        '--input', str(input_path), '--seed', str(seed), '--output', str(folder)],
                       check=True)
        report = json.loads((folder/'sample_report.json').read_text(encoding='utf-8'))
        with np.load(folder/'haze_dsc.npz', allow_pickle=False) as saved:
            if not np.array_equal(theta, saved['theta']) or not np.array_equal(phi, saved['phi']):
                raise ValueError('Sample angular grid mismatch')
            maps.append(saved['dsc'])
        sample_reports.append(report)
    effective_area = taper_effective_area(roughness)
    mean, standard_error, metrics = summarize(
        maps, reference, theta, effective_area, comparison_area_um2, roughness,
        finite_reference, finite_beam=finite_beam)
    rms_values = [item['surface_statistics']['triangle_rms_nm'] for item in sample_reports]
    from .haze_statistics import statistical_report, plot_statistics
    statistics = statistical_report(maps, theta, phi, 1. if finite_beam else comparison_area_um2/effective_area,
                                     apertures, tilts, target_relative_ci)
    metrics.update(statistics)
    metrics['realized_patch_rms_mean_nm'] = float(np.mean(rms_values))
    metrics['input_stationary_rms_nm'] = roughness.rms_nm
    warnings = []
    if roughness.model == 'gaussian' and roughness.width_um/roughness.correlation_length_um < 10:
        warnings.append('Patch width is below 10 Gaussian correlation lengths; finite-window PSD bias may remain.')
    if roughness.model == 'abc' and metrics['nyquist_variance_fraction'] < .9:
        warnings.append('Below 90% of the infinite ABC PSD variance lies within the inscribed height-grid Nyquist circle; high-frequency aliasing/discretization may dominate.')
    if roughness.model == 'abc' and roughness.psd_exponent <= 4:
        warnings.append('Infinite ABC PSD with C <= 4 has no finite RMS slope; FEM slopes depend on grid spacing.')
    if roughness.rms_nm/parameters.wavelength_nm >= .05:
        warnings.append('RMS/wavelength >= 0.05: the NIST first-order smooth-surface reference may be unreliable.')
    if abs(np.mean(rms_values)/roughness.rms_nm-1) > .1:
        warnings.append('Mean realized patch RMS differs from stationary input RMS by more than 10%.')
    if samples < 16:
        warnings.append('Fewer than 16 surfaces is an exploratory sample, not statistical convergence.')
    if roughness.model == 'measured_psd':
        from .experimental_surface import measurement_report
        metrics['measurement'] = measurement_report(roughness)
        warnings.append('Measured PSD is band-limited; synthesis assumes Gaussian heights. Unmeasured frequencies are zero.')
    warnings.append('Mesh, PML and patch-width convergence have not been established.')
    report = dict(status='complete', parameters=asdict(parameters), roughness=asdict(roughness),
                  seeds=[item['seed'] for item in sample_reports], metrics=metrics,
                  warnings=warnings, sample_hemisphere_dsc_um2=[
                      item['hemisphere_dsc_um2']*(1. if finite_beam else comparison_area_um2/effective_area)
                      for item in sample_reports],
                  nist_reference=reference_path is not None,
                  finite_patch_reference=finite_reference is not None,
                  definition='Mean of independent DSC intensities, never complex amplitudes. '
                             'Raw sample DSC is retained; FEM ensemble is divided by squared-window '
                             'effective roughness area, then scaled to the comparison area. '
                             'Finite-patch reference is the analytic nodal-covariance PSD of '
                             'the tapered piecewise-linear surface, passed through NIST Unit PSD '
                             'first-order response; infinite NIST statistical PSD is separate.',
                  figures=(['haze_ensemble.png', 'haze_ensemble_fem.png',
                            'haze_ensemble_stderr.png', 'haze_sample_convergence.png']+
                           (['haze_ensemble_nist.png'] if reference is not None else [])+
                           (['haze_ensemble_finite.png'] if finite_reference is not None else [])))
    from .illumination_output import illumination_metadata
    report['illumination'] = illumination_metadata(parameters, roughness)
    if finite_beam:
        from .finite_illumination import power_area
        report['definition'] = ('Mean of independently solved coherent finite-beam DSC intensities. '
            'Only the modeled rough patch scatters; no effective-area extrapolation or plane-wave NIST reference.')
        warnings += report['illumination']['warnings']
        metrics['collected_power_fraction_per_unit_laser_power'] = {
            name: item['mean_um2']/power_area(parameters) for name, item in statistics['observables'].items()}
    maps_to_save = dict(theta=theta, phi=phi, fem_mean=mean,
                        fem_standard_error=standard_error,
                        nist_psd=(reference if reference is not None else np.empty((0, 0))),
                        finite_patch=(finite_reference if finite_reference is not None
                                      else np.empty((0, 0))))
    if roughness.model == 'gaussian':
        maps_to_save['nist_gaussian'] = maps_to_save['nist_psd']
    np.savez_compressed(output/'haze_ensemble_maps.npz', **maps_to_save)
    from .haze_ensemble_plotting import plot_ensemble
    plot_ensemble(output/'haze_ensemble.png', theta, phi, mean, standard_error,
                  reference, report, finite_reference)
    plot_statistics(output/'haze_sample_convergence.png', statistics)
    (output/'haze_ensemble_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)
    return report


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'sample':
        parser = argparse.ArgumentParser(description='Solve one haze-only FEM surface')
        parser.add_argument('--input', required=True)
        parser.add_argument('--seed', type=int, required=True)
        parser.add_argument('--output', required=True)
        args = parser.parse_args(sys.argv[2:])
        inputs = json.loads(Path(args.input).read_text(encoding='utf-8'))
        p = Parameters(**inputs['parameters'])
        r = replace(Roughness(**inputs['roughness']), seed=args.seed)
        if inputs.get('measurement_sha256'):
            from .experimental_surface import input_fingerprint
            if input_fingerprint(r.measurement_path) != inputs['measurement_sha256']:
                raise ValueError('Measurement changed after the job was configured; recompute its references')
        simulate_sample(p, r, args.output, inputs['angular_step_deg'],
                        inputs['quadrature_order'])
        return
    parser = build_parser()
    parser.add_argument('--samples', type=int, default=16)
    parser.add_argument('--comparison-area-um2', type=float, default=.04)
    parser.add_argument('--taper-fraction', type=float, default=.125)
    parser.add_argument('--reference')
    parser.add_argument('--target-relative-ci', type=float, default=.05)
    args = parser.parse_args()
    p, r = scattering_settings(args, particle_enabled=False)
    r = replace(r, taper_fraction=args.taper_fraction)
    from .imaging_options import imaging_settings
    from .illumination_output import save_illumination
    settings = imaging_settings(args, incidence_deg=p.incidence_deg, parameters=p)
    report = run_ensemble(p, r, args.output, args.samples, args.comparison_area_um2,
                 args.reference, quadrature_order=args.quadrature_order or (7 if p.order >= 3 else 5),
                 apertures=args.na, tilts=args.tilts, target_relative_ci=args.target_relative_ci)
    report['illumination'] = save_illumination(args.output, p, settings, r)
    report['metrics']['mean_power_w'] = {name: item['mean_um2']*settings.irradiance_w_um2
        for name, item in report['metrics']['observables'].items()}
    report['illumination_settings'] = asdict(settings)
    (Path(args.output)/'haze_ensemble_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()

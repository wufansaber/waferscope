"""Bounded release experiments: run this module inside run_guarded."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time
import numpy as np


def beam_case(output):
    from .numerical_settings import numerical_parameters
    from .solver import solve, evaluate_fields
    from .theory_reference import mie_fields
    p = numerical_parameters('balanced', substrate_n=1., substrate_k=0., gap_nm=10.,
                             beam_waist_um=2., beam_order=5, threads=4)
    mesh, field, metrics = solve(p)
    # Non-interface, non-pole samples inside and outside the sphere.
    points = np.array(((.01, .012, .06), (.02, -.01, .08), (.09, .04, .05),
                       (-.08, .07, .11), (.03, -.08, .14), (.06, .06, .16)))
    e, h, _ = evaluate_fields(mesh, field, points, p)
    em, hm = mie_fields(points, p)
    er, hr = mie_fields(points, replace(p, beam_order=7))
    from .physical_audit import sphere_rule
    from .volume_sources import polarization_sources
    from .far_field import radiation_amplitude
    directions, weights = sphere_rule(8)
    source_points, source_weights, _ = polarization_sources(mesh, field, p, order=7)
    amplitude = radiation_amplitude(source_points, source_weights, directions, p)
    s = np.column_stack((-directions[:, 1], directions[:, 0], np.zeros(len(directions))))
    s /= np.linalg.norm(s, axis=1)[:, None]
    pol_p = np.cross(directions, s)
    references = []
    for radius in (200., 400.):
        scattered, _ = mie_fields(radius*directions, p, include_incident=False)
        far = radius*np.exp(-1j*p.k0*radius)*scattered
        references.append(np.column_stack((np.sum(far*s, axis=1), np.sum(far*pol_p, axis=1))))
    report = dict(parameters=asdict(p), solver=metrics,
        e_relative_l2=float(np.linalg.norm(e-em)/np.linalg.norm(em)),
        h_relative_l2=float(np.linalg.norm(h-hm)/np.linalg.norm(hm)),
        mie_quadrature_e_relative_l2=float(np.linalg.norm(em-er)/np.linalg.norm(er)),
        mie_quadrature_h_relative_l2=float(np.linalg.norm(hm-hr)/np.linalg.norm(hr)))
    report.update(far_amplitude_relative_l2=float(np.linalg.norm(amplitude-references[-1])/np.linalg.norm(references[-1])),
        far_reference_radius_change=float(np.linalg.norm(references[0]-references[1])/np.linalg.norm(references[1])),
        radiated_power_relative_error=float(abs(np.sum(weights*np.sum(abs(amplitude)**2, axis=1)) /
            np.sum(weights*np.sum(abs(references[-1])**2, axis=1))-1)))
    report['passed'] = max(report['e_relative_l2'], report['h_relative_l2'],
                           report['far_amplitude_relative_l2'], report['radiated_power_relative_error']) < .02
    np.savez_compressed(output/'fields.npz', points=points, e=e, h=h, mie_e=em, mie_h=hm,
                        directions=directions, amplitude=amplitude, mie_amplitude=references[-1])
    return report


def reimage_case(output, source):
    from .parameters import Parameters
    from .imaging import ImagingSettings
    from .imaging_run import image_saved_sources
    from .far_field import radiation_amplitude_scalar
    from .radiation import radiation_amplitude
    source = Path(source)
    report = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    scattering = Path(report['source_directory'])
    p = Parameters(**report['scattering_parameters'])
    with np.load(scattering/'paired_fields.npz') as data:
        points, weights = data['source_points_with'], data['source_polarization_weights_with']
    phi = np.linspace(0, 2*np.pi, 100, endpoint=False)
    directions = np.column_stack((.5*np.cos(phi), .5*np.sin(phi), np.full(100, np.sqrt(.75))))
    timings, values = {}, {}
    for name, function in (('scalar', radiation_amplitude_scalar), ('blocked', radiation_amplitude)):
        trials = []
        for _ in range(3):
            start = time.perf_counter()
            values[name] = function(points, weights, directions, p)
            trials.append(time.perf_counter()-start)
        timings[name] = trials
    settings = dict(report['settings'])
    settings.pop('detector', None)
    start = time.perf_counter()
    image_saved_sources(scattering, output, ImagingSettings(**settings), archive_mode='compact',
                        cache_directory=output/'pupil_cache')
    elapsed = time.perf_counter()-start
    start = time.perf_counter()
    warm = image_saved_sources(scattering, output/'warm', ImagingSettings(**settings), archive_mode='compact',
                               cache_directory=output/'pupil_cache')
    warm_elapsed = time.perf_counter()-start
    errors = {}
    with np.load(source/'image_fields.npz') as old, np.load(output/'image_fields.npz') as new:
        for key in new.files:
            if key in old and np.issubdtype(new[key].dtype, np.number):
                errors[key] = float(np.linalg.norm(new[key]-old[key])/max(np.linalg.norm(old[key]), 1e-300))
    result = dict(radiation_timings_seconds=timings,
        radiation_relative_l2=float(np.linalg.norm(values['blocked']-values['scalar'])/np.linalg.norm(values['scalar'])),
        compact_reimage_seconds=elapsed, original_bytes=(source/'image_fields.npz').stat().st_size,
        cached_reimage_seconds=warm_elapsed, cache_hits=[c['pupil_cache_hit'] for c in warm['channels']],
        compact_bytes=(output/'image_fields.npz').stat().st_size, saved_array_relative_l2=errors)
    result['passed'] = max(errors.values()) < 1e-11 and result['radiation_relative_l2'] < 1e-12
    return result


def sensitivity_case(output, variant, threads=4):
    from .numerical_settings import numerical_parameters
    from .rough_surface import Roughness, build_rough_mesh
    from .solver import solve
    from .volume_sources import polarization_sources
    from .imaging import ImagingSettings, calculate_channel
    p = numerical_parameters('balanced', gap_nm=10, threads=threads, pml_order=2, pml_inner_fraction=.5)
    roughness = Roughness()
    if variant == 'domain':
        p = replace(p, domain_radius_um=.24)
    elif variant == 'thickness':
        p = replace(p, pml_thickness_um=.20)
    elif variant.startswith('seed'):
        roughness = replace(roughness, seed=int(variant[4:]))
    mesh = build_rough_mesh(p, roughness)
    data, solvers = {}, {}
    for name, enabled in (('without', False), ('with', True)):
        _, field, info = solve(p, mesh, enabled)
        points, weights, _ = polarization_sources(mesh, field, p, order=7, particle_enabled=enabled)
        data['source_points_'+name], data['source_polarization_weights_'+name] = points, weights
        solvers[name] = info
        del field
    settings = ImagingSettings(pupil_samples=32, padding=2)
    channels = [calculate_channel(data, p, settings, i)[1] for i in range(3)]
    return dict(parameters=asdict(p), roughness=asdict(roughness), channels=channels, solvers=solvers,
                note='32-sample pupil sensitivity screen; not a mesh or statistical convergence certificate.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case', choices=('beam', 'reimage', 'sensitivity'))
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--source', default='build/thread_audit/pipeline8')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--variant', default='baseline', choices=('baseline', 'domain', 'thickness', 'seed8', 'seed9'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.case == 'beam':
        report = beam_case(args.output)
    elif args.case == 'reimage':
        report = reimage_case(args.output, args.source)
    else:
        report = sensitivity_case(args.output, args.variant, args.threads)
    (args.output/'audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)
    if report.get('passed') is False:
        raise SystemExit(1)


if __name__ == '__main__':
    main()

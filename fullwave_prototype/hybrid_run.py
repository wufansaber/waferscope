"""Whole-beam SPM with a shared-surface, phase-aligned local FEM correction."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters
from .far_field import angular_grid
from .finite_illumination import power_area
from .hybrid_radiation import HybridRadiation, triangle_source, triangle_amplitude
from .imaging import ImagingSettings, calculate_images
from .near_planes import saved_planes
from .parameters import Parameters
from .rough_run import decompose
from .rough_surface import Roughness
from .shared_surface import make_mother, save_mother, sample_mother, local_crop, load_mother
from .whole_surface import plan_surface, window


MODEL = 'whole_spm_local_fem_hybrid'
SCOPE = ('Whole-surface first-order SPM plus local FEM minus the same local triangular SPM. '
         'One shared synthetic surface and global complex phase reference. '
         'Not a rigorous full-wafer Maxwell solution; outer/local multiple scattering is omitted.')


def source_signature(metadata):
    values = {key: metadata[key] for key in ('parameters', 'roughness', 'stationary_roughness',
                                            'plan', 'source_sha256')}
    values['mother_sha256'] = metadata['mother_surface']['sha256']
    return hashlib.sha256(json.dumps(values, sort_keys=True).encode('utf-8')).hexdigest()


def optical_settings(values, p):
    values = dict(values)
    detector = DetectorParameters(**values.pop('detector', {}))
    settings = ImagingSettings(**values, detector=detector)
    settings.validate()
    if settings.illumination_model != 'gaussian' or settings.illumination_definition != p.beam_definition or not np.allclose(
            (settings.illumination_length_um, settings.illumination_width_um, settings.illumination_incidence_deg),
            (p.beam_length_um, p.beam_width_um, p.incidence_deg), rtol=1e-10, atol=0):
        raise ValueError('Hybrid scattering and imaging must use the same finite coherent beam')
    return replace(settings, illumination_normalization_area_um2=power_area(p),
                   irradiance_w_um2=settings.illumination_power_w/power_area(p))


def prepare_sources(inputs, output):
    from .rough_run import run
    p, rough = Parameters(**inputs['parameters']), Roughness(**inputs['roughness'])
    p.validate()
    optical_settings(inputs['imaging_settings'], p)
    p, plan = plan_surface(p, rough, inputs.get('coverage', .99), inputs.get('variance_fraction', .9),
                           inputs.get('taper_width_um', .025))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    mother_axis, mother_height = make_mother(rough, plan, inputs.get('mother_range_factor', 2))
    mother = save_mother(output, mother_axis, mother_height, rough)
    local = local_crop(output, mother_axis, mother_height, rough)
    local.validate(p)
    local_report = run(p, local, output/'local_scattering', angular_step=5,
                       quadrature_order=inputs.get('quadrature_order', 7 if p.order >= 3 else 5))
    with np.load(output/'local_scattering'/'paired_fields.npz', allow_pickle=False) as saved:
        data = {key: saved[key] for key in saved.files}
    n, step = plan['grid_points'], plan['spacing_um']
    axis = (np.arange(n)-n//2)*step
    height = sample_mother(mother_axis, mother_height, axis, axis)
    taper = window(axis, axis, plan['width_um'], plan['taper_width_um'])
    from .surface_beam import surface_fields_order
    background, ei, hi = surface_fields_order(axis, axis, p, plan.get('surface_beam_order', p.beam_order))
    flux = -np.real(np.cross(ei, hi.conj()))[..., 2]
    coverage = float(np.sum(flux*taper**2)*step**2/power_area(p))
    if not plan['target_power_coverage'] <= coverage <= 1.01:
        raise ValueError('Shared whole-surface illumination coverage failed on the final grid')
    plan['checked_taper_weighted_power_coverage'] = coverage
    positions, weights = triangle_source(data['surface_triangles_um'], p)
    data.update(whole_axis_um=axis, whole_height_um=height, whole_height_window=taper,
        whole_weighted_fields=height[..., None]*taper[..., None]*background,
        whole_incident_flux_unit=flux, local_spm_positions_um=positions, local_spm_weights=weights)
    # Check the local subtraction on a higher independent triangle quadrature.
    q = np.array([[.3, .4, np.sqrt(.75)], [-.4, .1, np.sqrt(.83)], [0., 0., 1.]])
    check_positions, check_weights = triangle_source(data['surface_triangles_um'], p, order=8)
    a = triangle_amplitude(positions, weights, q, p)
    b = triangle_amplitude(check_positions, check_weights, q, p)
    error = float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-300))
    if error > .01:
        raise ValueError('Local SPM subtraction quadrature did not reach the 1% target')
    np.savez_compressed(output/'hybrid_fields.npz', **data)
    metadata = dict(model=MODEL, source_format_version=2, scope=SCOPE, parameters=asdict(p),
        roughness=asdict(local), stationary_roughness=asdict(rough), plan=plan, mother_surface=mother,
        local_scattering_report=local_report, local_spm_triangle_relative_l2_change=error,
        source_sha256=hashlib.sha256((output/'hybrid_fields.npz').read_bytes()).hexdigest(),
        near_field_scope='Local paired FEM only; whole-surface outer near fields are not reconstructed',
        realization='One shared synthetic mother surface; NOT a statistically averaged wafer')
    metadata['physical_signature'] = source_signature(metadata)
    (output/'hybrid_source_report.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    (output/'hybrid_input.json').write_text(json.dumps(inputs, indent=2), encoding='utf-8')
    return data, metadata


def load_sources(source):
    source = Path(source)
    metadata = json.loads((source/'hybrid_source_report.json').read_text(encoding='utf-8'))
    Parameters(**metadata['parameters']).validate()
    if metadata.get('model') != MODEL or metadata.get('source_format_version') != 2:
        raise ValueError('Select a complete phase-aligned hybrid source')
    if source_signature(metadata) != metadata.get('physical_signature'):
        raise ValueError('Hybrid physical metadata changed; recompute scattering')
    _, _, mother = load_mother(source)
    if mother['sha256'] != metadata['mother_surface']['sha256']:
        raise ValueError('Shared surface no longer matches the paired FEM source')
    if hashlib.sha256((source/'hybrid_fields.npz').read_bytes()).hexdigest() != metadata['source_sha256']:
        raise ValueError('Hybrid complex source changed; do not reuse it for imaging')
    with np.load(source/'hybrid_fields.npz', allow_pickle=False) as saved:
        data = {key: saved[key] for key in saved.files}
    return data, metadata


def export_images(source, output, settings_values):
    data, metadata = load_sources(source)
    p = Parameters(**metadata['parameters'])
    settings = optical_settings(settings_values, p)
    radiation = HybridRadiation(data, p)
    theta, phi, q = angular_grid(5)
    terms = radiation.components(q)
    data.update(theta=theta, phi=phi, amplitude_with=terms['with_particle'], amplitude_without=terms['without'])
    components = decompose(data['amplitude_with'], data['amplitude_without'])
    data.update({'dsc_'+key: value for key, value in components.items()})
    change = data['amplitude_with']-data['amplitude_without']
    local_change = terms['local_fem_with']-terms['local_fem_without']
    change_error = float(np.linalg.norm(change-local_change)/max(np.linalg.norm(local_change), 1e-300))
    spm_check = float(np.linalg.norm(terms['local_fem_without']-terms['local_spm'])/
                      max(np.linalg.norm(terms['local_spm']), 1e-300))
    arrays, report = calculate_images(data, p, settings, retain_complex_fields=False,
                                      amplitude_provider=radiation.amplitude)
    from .hybrid_surface_audit import audit_surface
    print('Checking outer SPM range and spatial quadrature on the same mother surface...', flush=True)
    surface_audit = audit_surface(source, data, metadata, settings)
    report.update(model=MODEL, scope=SCOPE, near_field_scope=metadata['near_field_scope'],
        scattering_parameters=asdict(p), roughness=metadata['roughness'],
        stationary_roughness=metadata['stationary_roughness'], whole_surface_plan=metadata['plan'],
        mother_surface_sha256=metadata['mother_surface']['sha256'],
        source_physical_signature=metadata['physical_signature'], surface_convergence=surface_audit,
        image_channel_figures=[1, 2, 3], near_field_planes=[plane.metadata() for plane in saved_planes(data)],
        source_directory=str(Path(source).resolve()), local_change_relative_l2_error=change_error,
        local_fem_vs_spm_complex_relative_l2=spm_check,
        subtraction='A_hybrid = A_whole_SPM + A_local_FEM - A_same_local_SPM; '
                    'amplitudes share coordinates/polarization/incident flux; NOT intensity subtraction',
        status='experimental coherent hybrid; higher-order outer/local coupling and convergence unverified')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'image_fields.npz', **arrays)
    np.savez_compressed(output/'hybrid_far_fields.npz', theta=theta, phi=phi,
                        **terms, **{'dsc_'+key: value for key, value in components.items()})
    (output/'imaging_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    (output/'hybrid_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    (output/'hybrid_surface_audit.json').write_text(json.dumps(surface_audit, indent=2), encoding='utf-8')
    scattering = dict(metadata['local_scattering_report'], near_field_scope=metadata['near_field_scope'])
    from .imaging_plotting import plot_pipeline
    plot_pipeline(output, data, arrays, scattering, report)
    from .hybrid_plotting import plot_surface
    plot_surface(output, data, metadata, settings, surface_audit)
    if settings.detector.enabled:
        from .detector_plot import plot_detector
        plot_detector(output, arrays)
    print(json.dumps(dict(model=MODEL, image_power_w=[c['with']['image_power_w'] for c in report['channels']],
                          local_change_relative_l2_error=change_error,
                          local_fem_vs_spm_complex_relative_l2=spm_check), indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', help='Full physical input for a new paired FEM solve')
    parser.add_argument('--source', help='Reuse a verified hybrid source without a new FEM solve')
    parser.add_argument('--imaging-settings-json')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if bool(args.input) == bool(args.source):
        parser.error('Choose exactly one of --input or --source')
    if args.source:
        source = Path(args.source)
        settings = (json.loads(Path(args.imaging_settings_json).read_text(encoding='utf-8'))
                    if args.imaging_settings_json else
                    json.loads((source/'hybrid_input.json').read_text(encoding='utf-8'))['imaging_settings'])
    else:
        if args.imaging_settings_json:
            parser.error('New solves get their illumination and optics from --input')
        inputs = json.loads(Path(args.input).read_text(encoding='utf-8'))
        source = Path(args.output)
        prepare_sources(inputs, source)
        settings = inputs['imaging_settings']
    export_images(source, args.output, settings)


if __name__ == '__main__':
    main()

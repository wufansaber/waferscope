"""Compare independent particle-only and haze-only FEM experiments."""
import argparse
from dataclasses import asdict, replace
from datetime import datetime
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters
from .imaging import ImagingSettings, calculate_images
from .parameters import Parameters
from .rough_run import run as run_scattering
from .rough_surface import Roughness


def _read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def _source_for_result(result, image_report):
    candidates = [result/'scattering']
    recorded = image_report.get('source_directory')
    if recorded:
        source = Path(recorded)
        candidates.append(source if source.is_absolute() else result/source)
    for source in candidates:
        if (source/'report.json').is_file() and (source/'paired_fields.npz').is_file():
            report = _read_json(source/'report.json')
            if report.get('source_format_version') == 1:
                return source.resolve(), report
    raise ValueError('A reusable rough-surface scattering source is required.')


def _settings(saved):
    values = dict(saved)
    values['detector'] = DetectorParameters(**values['detector'])
    for name in ('numerical_apertures', 'tilt_degrees', 'aberrations_nm'):
        values[name] = tuple(values[name])
    settings = ImagingSettings(**values)
    settings.validate()
    return replace(settings, detector=replace(settings.detector, enabled=False))


def _clean_source(result, source, parameters, roughness, angular_step, quadrature_order):
    candidates = [source.parent/'particle_only_scattering', result/'particle_only_scattering']
    for candidate in candidates:
        try:
            report = _read_json(candidate/'report.json')
            if (report.get('source_format_version') == 1 and
                    report['parameters'] == asdict(parameters) and
                    report['roughness']['rms_nm'] == 0 and
                    report['angular_step_deg'] == angular_step and
                    report['metrics']['with']['source_quadrature_order'] == quadrature_order and
                    (candidate/'paired_fields.npz').is_file()):
                return candidate, report
        except (OSError, ValueError, KeyError, TypeError):
            pass
    destination = result/'particle_only_scattering'
    report = run_scattering(parameters, replace(roughness, rms_nm=0.), destination,
                            angular_step=angular_step, quadrature_order=quadrature_order)
    return destination, report


def compare(result, cache_directory=None, imaging_settings_path=None):
    result = Path(result).resolve()
    image_report = _read_json(result/'imaging_report.json')
    source, rough_report = _source_for_result(result, image_report)
    parameters = Parameters(**rough_report['parameters'])
    if image_report['scattering_parameters'] != asdict(parameters):
        raise ValueError('Scattering source and image settings describe different physical parameters.')
    roughness = Roughness(**rough_report['roughness'])
    requested_settings = (_read_json(imaging_settings_path) if imaging_settings_path else
                          image_report['settings'])
    settings = _settings(requested_settings)
    angular_step = rough_report['angular_step_deg']
    quadrature_order = rough_report['metrics']['with']['source_quadrature_order']
    clean_source, clean_report = _clean_source(
        result, source, parameters, roughness, angular_step, quadrature_order)
    with np.load(source/'paired_fields.npz', allow_pickle=False) as saved:
        rough_data = {name: saved[name] for name in saved.files}
    rough_images, haze_report = calculate_images(
        rough_data, parameters, settings, retain_complex_fields=False,
        cache_directory=cache_directory)
    haze_fields = {key: rough_images[key] for index in (1, 2, 3)
                   for key in (f'ch{index}_image_axis_um', f'ch{index}_image_without_w_um2')}
    del rough_images, rough_data
    with np.load(clean_source/'paired_fields.npz', allow_pickle=False) as saved:
        clean_data = {name: saved[name] for name in saved.files}
    particle_images, particle_report = calculate_images(
        clean_data, parameters, settings, retain_complex_fields=False,
        cache_directory=cache_directory)

    output = result/'component_comparison'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    output.mkdir(parents=True, exist_ok=False)
    fields = {}
    channels = []
    exposure_s = settings.detector.exposure_ms*1e-3
    for index in (1, 2, 3):
        axis = haze_fields[f'ch{index}_image_axis_um']
        clean_axis = particle_images[f'ch{index}_image_axis_um']
        if axis.shape != clean_axis.shape or not np.allclose(axis, clean_axis, rtol=1e-10, atol=1e-10):
            raise ValueError('Particle and haze images do not share the same image coordinates.')
        particle = particle_images[f'ch{index}_image_with_w_um2']
        haze = haze_fields[f'ch{index}_image_without_w_um2']
        if not np.isfinite(particle).all() or not np.isfinite(haze).all():
            raise ValueError('Non-finite component image.')
        area = (axis[1]-axis[0])**2
        particle_power = particle_report['channels'][index-1]['with']['image_power_w']
        haze_power = haze_report['channels'][index-1]['without']['image_power_w']
        for name, value, power in (('particle', particle, particle_power),
                                   ('haze', haze, haze_power)):
            if not np.isclose(value.sum()*area, power, rtol=1e-8, atol=1e-25):
                raise ValueError(name+' image power fails integral check.')
            fields[f'ch{index}_{name}_w_um2'] = value
        fields[f'ch{index}_axis_um'] = axis
        channels.append(dict(channel=index, na=settings.numerical_apertures[index-1],
                             tilt_deg=settings.tilt_degrees[index-1],
                             particle_power_w=float(particle_power), haze_power_w=float(haze_power),
                             particle_energy_j=float(particle_power*exposure_s),
                             haze_energy_j=float(haze_power*exposure_s)))
    np.savez_compressed(output/'image_components.npz', **fields)
    with np.load(clean_source/'paired_fields.npz', allow_pickle=False) as particle_far, \
            np.load(source/'paired_fields.npz', allow_pickle=False) as haze_far:
        theta, phi = particle_far['theta'], particle_far['phi']
        if not np.array_equal(theta, haze_far['theta']) or not np.array_equal(phi, haze_far['phi']):
            raise ValueError('Particle and haze far fields do not share the same angular grid.')
        far = dict(theta=theta, phi=phi, particle_dsc=particle_far['dsc_total'],
                   haze_dsc=haze_far['dsc_background'])
    np.savez_compressed(output/'far_components.npz', **far)
    particle_integral = clean_report['metrics']['hemisphere_integrals_um2']['total']
    haze_integral = rough_report['metrics']['hemisphere_integrals_um2']['background']
    from .run import hemisphere_integral
    for name, integral in (('particle', particle_integral), ('haze', haze_integral)):
        if not np.isclose(hemisphere_integral(far[name+'_dsc'], theta), integral,
                          rtol=1e-8, atol=1e-20):
            raise ValueError(name+' far-field integral check failed.')
    comparison = dict(status='complete', definition=dict(
        particle='Sphere above a smooth half-space; no roughness.',
        haze='The saved rough surface without a particle.',
        warning='Independent experiments; their powers are not additive and do not reconstruct the sphere-on-rough-surface result.'),
        parameters=asdict(parameters), roughness=asdict(roughness),
        imaging_settings=asdict(settings),
        source_directories=dict(particle=str(clean_source), haze=str(source)),
        far_hemisphere=dict(particle_dsc_um2=float(particle_integral),
                            haze_dsc_um2=float(haze_integral),
                            particle_equivalent_power_w=float(particle_integral*settings.irradiance_w_um2),
                            haze_equivalent_power_w=float(haze_integral*settings.irradiance_w_um2)),
        exposure_ms=settings.detector.exposure_ms, channels=channels,
        figures=dict(far='far_components.png', far_log='far_components_log.png',
                     image='image_components.png', image_log='image_components_log.png'),
        scope='Only scattered light; no specular leakage. FEM/PML/mesh convergence not guaranteed.')
    from .component_compare_plotting import plot_comparison
    plot_comparison(output, far, fields, comparison)
    (output/'comparison.json').write_text(json.dumps(comparison, indent=2), encoding='utf-8')
    manifest = dict(directory=str(output.relative_to(result)), report='comparison.json',
                    totals='totals.txt', figures=comparison['figures'])
    temporary = result/'component_comparison_manifest.tmp'
    temporary.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    temporary.replace(result/'component_comparison_manifest.json')
    print(json.dumps(comparison, indent=2), flush=True)
    return comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', required=True)
    parser.add_argument('--pupil-cache')
    parser.add_argument('--imaging-settings-json')
    args = parser.parse_args()
    compare(args.result, args.pupil_cache, args.imaging_settings_json)


if __name__ == '__main__':
    main()

"""Verify archived 100 nm images and camera SNR without changing illumination."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

import numpy as np

from detector_model import DetectorParameters
from fullwave_prototype.imaging_audit import check_saved
from tools.particle100_camera import camera_case, independent_budget, monte_carlo, direct_pixels


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2), encoding='utf-8')


def load_image(source):
    report = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    if report['scattering_parameters']['diameter_nm'] != 100:
        raise ValueError('Only the 100 nm benchmark is accepted')
    with np.load(source/'image_fields.npz', allow_pickle=False) as saved:
        arrays = {k: saved[k] for k in saved.files if '_detector_' not in k}
    return arrays, report


def export_case(output, arrays, report, parameters, transmission=1):
    from fullwave_prototype.detector_plot import plot_detector
    output.mkdir(parents=True, exist_ok=True)
    summary, detectors = camera_case(arrays, report, parameters, transmission)
    camera_arrays = {'ch{}_detector_{}'.format(i+1, k): v for i, d in enumerate(detectors)
                     for k, v in d.items() if isinstance(v, np.ndarray)}
    plot_detector(output, camera_arrays)
    np.savez_compressed(output/'camera_fields.npz', **camera_arrays)
    save_json(output/'camera_report.json', summary)
    return summary, detectors


def refine_hybrid(source, output, samples):
    from fullwave_prototype.hybrid_run import load_sources, optical_settings, SCOPE
    from fullwave_prototype.hybrid_radiation import HybridRadiation
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.imaging import calculate_images
    _, original = load_image(source)
    data, metadata = load_sources(source)
    p = Parameters(**metadata['parameters'])
    settings = optical_settings(original['settings'], p)
    settings = replace(settings, pupil_samples=samples, padding=8, detector=DetectorParameters())
    start = time.perf_counter()
    print('Refining 100 nm hybrid pupil to {} samples...'.format(samples), flush=True)
    arrays, report = calculate_images(data, p, settings, retain_complex_fields=False,
                                     amplitude_provider=HybridRadiation(data, p).amplitude)
    report.update(scattering_parameters=asdict(p), roughness=metadata['roughness'],
        stationary_roughness=metadata['stationary_roughness'], model=metadata['model'], scope=SCOPE,
        optical_source=str(source.resolve()), elapsed_seconds=time.perf_counter()-start,
        near_field_scope='Reused archived local FEM sources; no new FEM solve',
        source_physical_signature=metadata['physical_signature'])
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'image_fields.npz', **arrays)
    save_json(output/'imaging_report.json', report)
    print('Pupil {} complete in {:.1f}s'.format(samples, report['elapsed_seconds']), flush=True)
    return arrays, report


def plane_check(source, output):
    arrays, report = load_image(source)
    propagation, _ = check_saved(source)
    parameters = replace(DetectorParameters(**report['settings']['detector']), exposure_ms=.01,
                         reference_frames=1)
    summary, detectors = export_case(output/'camera_10us', arrays, report, parameters)
    checks = [dict(channel=i+1, independent_budget=independent_budget(d, parameters,
        report['scattering_parameters']['wavelength_nm']), monte_carlo=monte_carlo(d, parameters, seed=781+i),
        pixel_integration=direct_pixels(arrays, report, d, i, parameters)) for i, d in enumerate(detectors)]
    saturated, _ = export_case(output/'camera_1ms', arrays, report, replace(parameters, exposure_ms=1))
    payload = dict(source=str(source.resolve()), conditions=report['settings'],
        scattering=report['scattering_parameters'], propagation=propagation,
        camera_10us=summary, camera_1ms=saturated, checks=checks)
    save_json(output/'report.json', payload)
    return payload


def hybrid_check(source, output):
    levels = [refine_hybrid(source, output/('pupil'+str(n)), n) for n in (64, 128)]
    parameters = DetectorParameters(enabled=True, pixels_x=16, pixels_y=16, exposure_ms=.01,
                                    prnu_fraction=.01, dsnu_e=.5, reference_frames=1)
    camera_levels = [camera_case(a, r, parameters)[1] for a, r in levels]
    finest, report = levels[-1]
    original, _ = export_case(output/'camera_original_1W_10us', finest, report, parameters)
    peak = max(float(max(d['mean_total_e'].max(),
                        (d['haze_signal_e']+d['mean_background_e']).max()))
               for d in camera_levels[-1])
    # Illustrative uniform optical attenuation, never applied to the original result.
    attenuation = min(1., .5*parameters.full_well_e/peak)
    attenuated, detectors = export_case(output/'camera_attenuated_10us', finest, report, parameters, attenuation)
    coarse_detectors = camera_case(*levels[0], parameters, attenuation)[1]
    convergence = []
    for i in range(3):
        coarse, fine = camera_levels[0][i], camera_levels[1][i]
        delta_coarse, delta_fine = coarse['defect_contrast_e'], fine['defect_contrast_e']
        convergence.append(dict(channel=i+1,
            total_pixel_power_relative_l2=float(np.linalg.norm(coarse['pixel_power_w']-fine['pixel_power_w'])/
                                               np.linalg.norm(fine['pixel_power_w'])),
            contrast_relative_l2=float(np.linalg.norm(delta_coarse-delta_fine)/np.linalg.norm(delta_fine)),
            template_snr_relative_change=abs(coarse_detectors[i]['detection_template_snr']/
                                              detectors[i]['detection_template_snr']-1),
            collected_power_relative_change=abs(levels[0][1]['channels'][i]['with']['image_power_w']/
                                                 report['channels'][i]['with']['image_power_w']-1)))
    propagation, _ = check_saved(output/'pupil128')
    checks = [dict(channel=i+1, independent_budget=independent_budget(d, parameters,
        report['scattering_parameters']['wavelength_nm']), monte_carlo=monte_carlo(d, parameters, seed=901+i),
        pixel_integration=direct_pixels(finest, report, camera_levels[-1][i], i, parameters))
        for i, d in enumerate(detectors)]
    payload = dict(source=str(source.resolve()), scattering=report['scattering_parameters'],
        conditions=report['settings'], convergence=convergence, propagation=propagation,
        original_1W_10us=original, extra_attenuation_example=attenuated,
        extra_transmission=attenuation, checks=checks,
        warning='Gaussian illumination unchanged; archived first-order local FEM + outer SPM, single surface. '
            'Material n,k and camera QE are assumed inputs. Not a certified full-wafer model or detection probability.')
    save_json(output/'report.json', payload)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='build/particle100_validation_20261016')
    parser.add_argument('--mode', choices=('plane', 'hybrid', 'both'), default='both')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    result = {}
    if args.mode in ('plane', 'both'):
        print('Checking 355 nm / 100 nm local FEM baseline...', flush=True)
        result['plane'] = plane_check(Path('build/camera_validation_20261012/fem_rough'), output/'plane355')
    if args.mode in ('hybrid', 'both'):
        result['hybrid'] = hybrid_check(Path('build/model_campaign_20261010/large_hybrid'), output/'hybrid266')
    from tools.particle100_plot import plot_summary
    plot_summary(output, result)
    save_json(output/'report.json', result)
    failed = [c for section in result.values() for c in section['checks']
              if not c['monte_carlo']['passed'] or c['independent_budget']['variance_relative_max_error'] > 1e-12]
    failed += [c for section in result.values() for c in section['propagation']['channels'] if not c['passed']]
    if failed:
        raise RuntimeError('Propagation or camera validation failed; inspect report.json')
    print('100 nm validation complete: '+str(output.resolve()), flush=True)


if __name__ == '__main__':
    main()

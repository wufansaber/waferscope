"""Attach paired camera/SNR results to a saved FEM or hybrid optical image."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import shutil

import numpy as np

from defect_detection import calculate_detection
from detector_model import DetectorParameters
from gray_acquisition import acquire_gray, validate_acquisition
from camera_profiles import describe_profile
from .detector_plot import plot_detector


def run(source, output, exposure_ms=None, reference_frames=1, haze_mismatch_fraction=0, camera_settings=None):
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or source in output.parents or output in source.parents:
        raise ValueError('Replay output must be separate from the source tree')
    report = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    parameters = replace(DetectorParameters(**report['settings'].get('detector', {})), enabled=True,
                         reference_frames=reference_frames, haze_mismatch_fraction=haze_mismatch_fraction)
    camera_settings = dict(camera_settings or {})
    if 'detector' in camera_settings:
        parameters = DetectorParameters(**camera_settings['detector'])
        parameters = replace(parameters, enabled=True)
    if exposure_ms is not None:
        parameters = replace(parameters, exposure_ms=exposure_ms)
    parameters.validate()
    pairs = camera_settings.get('gray_snr_pairs', 0)
    seed = camera_settings.get('gray_snr_seed', 20261016)
    profile = camera_settings.get('camera_profile', 'custom')
    describe_profile(profile)
    validate_acquisition(parameters, pairs, seed)
    if profile == 'gemi_area' and report['scattering_parameters']['wavelength_nm'] != 266:
        raise ValueError('Gemi preset QE is specified at 266 nm; use a custom calibrated QE at other wavelengths')
    from .parameters import Parameters
    from .field_coverage import illumination_coverage
    physical = Parameters(**report['scattering_parameters'])
    settings = report['settings']
    coverage = illumination_coverage(physical, parameters, settings['magnification'],
        settings['numerical_apertures'], settings['tilt_degrees'], settings.get('minimum_field_illumination', .5))
    if settings.get('illumination_fit') == 'camera' and not coverage['all_channels_covered']:
        raise ValueError('The saved beam does not cover the new sensor; calculate new scattering')
    seeds = np.random.SeedSequence(parameters.seed).spawn(3)
    with np.load(source/'image_fields.npz', allow_pickle=False) as saved:
        arrays = {k: saved[k] for k in saved.files if '_detector_' not in k}
    for i, channel in enumerate(report['channels']):
        prefix = 'ch{}_'.format(i+1)
        def optical(name):
            return dict(image_axis_um=arrays[prefix+'image_axis_um'],
                irradiance_w_um2=arrays[prefix+'image_'+name+'_w_um2'],
                power_w=channel[name]['image_power_w'])
        detector = calculate_detection(optical('with'), optical('without'),
            report['scattering_parameters']['wavelength_nm'], parameters, np.random.default_rng(seeds[i]), i)
        detector.update(acquire_gray(detector, parameters, pairs, seed, i))
        detector['camera_profile_note'] = describe_profile(profile)
        channel['detector'] = {k: v for k, v in detector.items() if not isinstance(v, np.ndarray)}
        arrays.update({prefix+'detector_'+k: v for k, v in detector.items() if isinstance(v, np.ndarray)})
    report['settings']['detector'] = asdict(parameters)
    report['settings'].update(camera_profile=profile, gray_snr_pairs=pairs, gray_snr_seed=seed)
    report['field_coverage'] = coverage
    report.update(camera_detection_model_version=1, camera_optical_source=str(source),
                  camera_replay_scope='Original paired optical irradiance unchanged; new camera/reference frames only.')
    output.mkdir(parents=True, exist_ok=True)
    # Retain the scientific figures required by the existing UI result loader.
    for path in source.glob('*.png'):
        if not path.name.startswith(('detector', 'gray_')):
            shutil.copy2(path, output/path.name)
    for filename in ('illumination_report.json',):
        if (source/filename).is_file():
            shutil.copy2(source/filename, output/filename)
    plot_detector(output, arrays)
    np.savez_compressed(output/'image_fields.npz', **arrays)
    (output/'imaging_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--exposure-ms', type=float)
    parser.add_argument('--reference-frames', type=int, default=1)
    parser.add_argument('--haze-mismatch-fraction', type=float, default=0)
    parser.add_argument('--camera-settings-json')
    args = parser.parse_args()
    camera = json.loads(Path(args.camera_settings_json).read_text(encoding='utf-8')) if args.camera_settings_json else None
    report = run(args.source, args.output, args.exposure_ms, args.reference_frames, args.haze_mismatch_fraction, camera)
    print(json.dumps([dict(channel=i+1, **{k: c['detector'][k] for k in
        ('defect_contrast_roi_e', 'haze_roi_e', 'detection_roi_snr', 'detection_template_snr',
         'detection_invalid_pixels')}) for i, c in enumerate(report['channels'])], indent=2))


if __name__ == '__main__':
    main()

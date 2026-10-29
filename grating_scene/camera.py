"""Paired area-camera observations of a saved periodic S4 scene, no TDI."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil

import numpy as np

from camera_readout import export_raw_camera
from defect_detection import contrast_metrics, reference_average
from detector_model import DetectorParameters, calculate_pixel_detector
from gray_acquisition import acquire_gray, validate_acquisition
from .imaging import mode_coefficients
from .pipeline import digest, is_within, load_scene, write_json
from .pixels import integrate_modes


def run_camera(source, output, settings):
    source, output = Path(source).resolve(), Path(output).resolve()
    if is_within(output, source) or is_within(source, output):
        raise ValueError('Camera output must be separate from its optical source tree')
    config, optical, arrays = load_scene(source)
    p = DetectorParameters(**settings['detector'])
    p.validate()
    pairs, seed = settings.get('gray_pairs', 1000), settings.get('gray_seed', 20261027)
    if not p.enabled or p.haze_mismatch_fraction:
        raise ValueError('Enable the camera; this deterministic grating scene has no haze mismatch model')
    validate_acquisition(p, pairs, seed)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Use an empty camera folder')
    output.mkdir(parents=True, exist_ok=True)
    snapshot = output/'optical_scene'
    shutil.copytree(source, snapshot)
    write_json(output/'camera_settings.json', settings)
    channels, raw = [], {}
    streams = np.random.SeedSequence(p.seed).spawn(3)
    data = {}
    for label in ('reference', 'sample'):
        with np.load(source/label/'modes.npz', allow_pickle=False) as saved:
            data[label] = dict(saved)
    for i, (mask, tilt) in enumerate(zip(arrays['far']['channel_masks'], config.tilts_deg)):
        print('CAMERA_CHANNEL '+str(i+1), flush=True)
        observations = {}
        rng = np.random.default_rng(streams[i])
        for label in ('sample', 'reference'):
            orders, coefficients = mode_coefficients(data[label], mask, config, tilt)
            x, y, power = integrate_modes(orders, coefficients, config, p)
            observations[label] = calculate_pixel_detector(x, y, power, config.scene.wavelength_um*1000, p, rng, i)
        sample, reference = observations['sample'], observations['reference']
        average, valid = reference_average(reference, p, rng)
        sample.update(contrast_metrics(sample, reference, p, valid))
        sample.update(reference_mean_dn=average, reference_digital_dn=reference['digital_dn'],
            difference_dn=sample['digital_dn'].astype(float)-average)
        sample.update(acquire_gray(sample, p, pairs, seed, i))
        # The shared readout calls its reference "haze"; here it is the deterministic bare grating.
        names = {'haze_signal_e': 'reference_signal_e', 'haze_roi_e': 'reference_roi_e',
            'reference_haze_shot_variance_e2': 'reference_shot_variance_e2',
            'reference_haze_roi_variance_e2': 'reference_shot_roi_variance_e2'}
        sample = {names.get(k, k): v for k, v in sample.items() if not k.startswith('haze_mismatch')
                  and k != 'detection_scope'}
        channel = {k: v for k, v in sample.items() if not isinstance(v, np.ndarray)}
        channel.update(channel=i+1, enabled=config.numerical_apertures[i] > 0,
            reference_roi_power_w=reference['roi_power_w'],
            ideal_relay_periods_in_roi=(p.pixels_x*p.pixel_pitch_um/(config.magnification*config.scene.lx_um),
                p.pixels_y*p.pixel_pitch_um/(config.magnification*config.scene.ly_um)))
        raw.update({'ch{}_{}'.format(i+1, k): v for k, v in sample.items() if isinstance(v, np.ndarray)})
        export_raw_camera(output/('ch'+str(i+1)), [sample, reference])
        (output/('ch'+str(i+1))/'camera_ch1.tiff').rename(output/('ch'+str(i+1))/'sample_raw.tiff')
        (output/('ch'+str(i+1))/'camera_ch2.tiff').rename(output/('ch'+str(i+1))/'reference_first_raw.tiff')
        channels.append(channel)
    np.savez_compressed(output/'camera_data.npz', **raw)
    report = dict(format='waferscope-periodic-camera', version=1, status='completed',
        settings=dict(detector=asdict(p), gray_pairs=pairs, gray_seed=seed), channels=channels,
        optical_report_sha256=digest(source/'scene_report.json'), optical_model_validated=False,
        scope='Simulated area camera of an ideal periodic wafer-conjugate relay; no TDI or actual camera measurement.',
        limitations=optical['limitations']+[
            'The ROI can span repeated cells. Those are physical copies of this periodic model, not isolated defects.',
            'Reference is the deterministic bare grating, not haze. No roughness, motion or registration noise.',
            'Raw DN temporal ROI SNR is conditional on this optical model and sensor inputs, not optical validation.'])
    from .camera_display import brightness_rows
    from .camera_plot import camera_figure
    from matplotlib.figure import Figure
    report['brightness'] = brightness_rows(raw, report)
    figure = Figure(figsize=(12, 11))
    camera_figure(figure, raw, report)
    figure.savefig(output/'camera.png', dpi=150)
    report['files_sha256'] = {str(path.relative_to(output)).replace('\\', '/'): digest(path)
        for path in sorted(output.rglob('*')) if path.is_file()}
    write_json(output/'camera_report.json', report)
    return report


def load_camera(folder):
    folder = Path(folder).resolve()
    report = json.loads((folder/'camera_report.json').read_text(encoding='utf-8'))
    if (report.get('format') != 'waferscope-periodic-camera' or report.get('version') != 1
            or report.get('status') != 'completed'):
        raise ValueError('Select a completed periodic camera result')
    required = {'camera_data.npz', 'camera_settings.json', 'optical_scene/scene_report.json'}
    if not required.issubset(report['files_sha256']):
        raise ValueError('Camera bundle lacks original data')
    for name, expected in report['files_sha256'].items():
        path = (folder/name).resolve()
        if not is_within(path, folder) or digest(path) != expected:
            raise ValueError('Camera file path/hash check failed: '+name)
    load_scene(folder/'optical_scene')
    if digest(folder/'optical_scene/scene_report.json') != report['optical_report_sha256']:
        raise ValueError('Camera optical provenance mismatch')
    with np.load(folder/'camera_data.npz', allow_pickle=False) as saved:
        arrays = dict(saved)
    return report, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--settings', type=Path, required=True)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--cancel-file', type=Path)
    parser.add_argument('--timeout', type=float, default=600.)
    parser.add_argument('--max-memory-gib', type=float, default=8.)
    args = parser.parse_args()
    settings = json.loads(args.settings.read_text(encoding='utf-8'))
    if args.worker:
        run_camera(args.source, args.output, settings)
        return
    from fullwave_prototype.guarded_run import run_guarded
    import sys
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = run_guarded([sys.executable, '-m', 'grating_scene.camera', '--worker']+sys.argv[1:],
        args.output.with_suffix('.resources.json'), maximum=args.max_memory_gib, reserve=3.,
        timeout=args.timeout, cancel_file=args.cancel_file)
    if result['status'] != 'completed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()

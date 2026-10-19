"""Save simulated raw gray stacks, reload them and measure empirical SNR/NSR."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from camera_readout import digitize
from detector_model import DetectorParameters
from tools.gray_frame_snr import measure_paired_frames


def sample_stack(detector, parameters, frames, rng, reference=False):
    mean = (detector['haze_signal_e'] if reference else detector['mean_signal_e'])
    total = mean+detector['mean_background_e']
    shape = (frames,)+mean.shape
    charge = rng.poisson(np.broadcast_to(total, shape))
    clipped = int(np.count_nonzero(charge >= parameters.full_well_e))
    measured = np.minimum(charge, parameters.full_well_e).astype(float)
    measured += rng.normal(0, parameters.read_noise_e, shape)
    measured += detector['dsnu_offset_e']
    raw = digitize(measured, total, parameters, detector['dsnu_offset_e'])
    return raw['digital_dn'], dict(full_well_clipped_values=clipped,
        adc_clipped_low_values=raw['adc_clipped_low_pixels'],
        adc_clipped_high_values=raw['adc_clipped_high_pixels'])


def camera_states(source, parameters=None):
    report = json.loads((source/'imaging_report.json').read_text(encoding='utf-8'))
    with np.load(source/'image_fields.npz', allow_pickle=False) as saved:
        if parameters is None:
            states = []
            for i in (1, 2, 3):
                prefix = 'ch{}_detector_'.format(i)
                states.append({name[len(prefix):]: saved[name] for name in saved.files if name.startswith(prefix)})
        else:
            from tools.particle100_camera import camera_case
            arrays = {k: saved[k] for k in saved.files if '_detector_' not in k}
            _, states = camera_case(arrays, report, parameters)
    p = parameters or DetectorParameters(**report['settings']['detector'])
    report['settings']['detector'] = asdict(p)
    return report, p, states


def generate_stacks(source, output, frames, seed, parameters=None):
    report, p, states = camera_states(source, parameters)
    if p.reference_frames != 1 or p.haze_mismatch_fraction != 0:
        raise ValueError('This acquisition example needs one same-surface reference per sample frame')
    seeds = np.random.SeedSequence(seed).spawn(6)
    arrays, status = {}, []
    for i, d in enumerate(states, 1):
        if not d['detection_valid'].all():
            raise ValueError('Camera output includes invalid pixels')
        for j, name in enumerate(('sample', 'reference')):
            raw, check = sample_stack(d, p, frames, np.random.default_rng(seeds[2*(i-1)+j]), j == 1)
            raw[0] = d['digital_dn' if j == 0 else 'reference_digital_dn']
            arrays['ch{}_{}_dn'.format(i, name)] = raw
            status.append(dict(channel=i, state=name, **check))
    if any(c[k] for c in status for k in ('full_well_clipped_values', 'adc_clipped_low_values', 'adc_clipped_high_values')):
        raise ValueError('Simulated acquisition clipped; adjust exposure/transmission and acquire again')
    path = output/'raw_gray_frames.npz'
    np.savez_compressed(path, **arrays)
    return report, status, path


def plot_frames(output, saved, metrics, offset, bit_depth=16):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from PIL import Image
    maximum = max(float(saved['ch{}_{}_dn'.format(i, state)][0].max())
                  for i in (1, 2, 3) for state in ('sample', 'reference'))
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    for i in (1, 2, 3):
        for row, state in enumerate(('sample', 'reference')):
            raw = saved['ch{}_{}_dn'.format(i, state)][0]
            Image.fromarray(raw).save(output/'ch{}_{}_raw16.tiff'.format(i, state))
            im = axes[row, i-1].imshow(raw, origin='lower', cmap='gray', vmin=offset, vmax=maximum,
                                      interpolation='nearest')
            axes[row, i-1].set(title='Ch{} | {}'.format(i, state), xlabel='Pixel X', ylabel='Pixel Y')
    fig.colorbar(im, ax=list(axes.ravel()), label='Raw DN | shared display window; no change to stored data')
    fig.suptitle('Simulated {}-bit ADC gray frames | uint16 TIFF storage | first pair\n'.format(bit_depth)+
                 '{}-pair empirical ROI SNR: '.format(metrics[0]['full_roi']['frames'])+
                 ', '.join('{:.2f}'.format(m['full_roi']['snr']) for m in metrics))
    fig.savefig(output/'raw_gray_frames.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', default='build/field_coverage_20261016/result')
    parser.add_argument('--output', default='build/field_coverage_20261016/gray_snr')
    parser.add_argument('--frames', type=int, default=5000)
    parser.add_argument('--seed', type=int, default=20261016)
    parser.add_argument('--camera-config', help='JSON containing detector parameters and their provenance')
    args = parser.parse_args()
    if not 100 <= args.frames <= 20000:
        raise ValueError('Use 100..20000 frame pairs')
    source, output = Path(args.source), Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(args.camera_config).read_text(encoding='utf-8')) if args.camera_config else None
    parameters = DetectorParameters(**config['detector']) if config else None
    if parameters:
        parameters.validate()
    report, acquisition, raw_path = generate_stacks(source, output, args.frames, args.seed, parameters)
    channels = []
    # Estimation starts from reloaded integer DN; optical/electron expectation arrays are not used.
    with np.load(raw_path, allow_pickle=False) as saved:
        for i in (1, 2, 3):
            sample, reference = (saved['ch{}_{}_dn'.format(i, name)] for name in ('sample', 'reference'))
            full, series = measure_paired_frames(sample, reference,
                bit_depth=report['settings']['detector']['bit_depth'])
            np.savetxt(output/'ch{}_roi_gray_series.csv'.format(i), series, delimiter=',',
                       header='sample_roi_dn,reference_roi_dn,difference_roi_dn', comments='')
            channels.append(dict(channel=i, full_roi=full))
        plot_frames(output, saved, channels, report['settings']['detector']['offset_dn'],
                    report['settings']['detector']['bit_depth'])
    result = dict(provenance='SIMULATED camera raw DN frames, NOT hardware camera measurements',
        source=str(source.resolve()), raw_frames=str(raw_path.resolve()), frame_pairs=args.frames,
        acquisition_seed=args.seed, settings=report['settings'], acquisition_clipping=acquisition,
        channels=channels, camera_profile=config,
        limitations=['Same static synthetic haze surface in each pair; spatial haze mismatch is not sampled.',
                     'Temporal empirical ROI SNR, not single-image spatial contrast or detection probability.',
                     'Source full-wave/outer-surface physical convergence is not certified.',
                     'Frame zero uses the original readout if no override, otherwise the new camera configuration.'])
    (output/'gray_snr_report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result['channels'], indent=2))


if __name__ == '__main__':
    main()

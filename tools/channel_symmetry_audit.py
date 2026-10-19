"""Audit the current symmetric channels: flat wafer, mirror remesh, independent haze."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters, integrate_pixels
from fullwave_prototype.hybrid_run import load_sources, optical_settings
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.surface_beam import surface_fields_order
from .channel_symmetry import relative_error, reflect_centered_grid
from .channel_symmetry_checks import smooth_check, mirror_check
from .channel_symmetry_ensemble import run_ensemble
from .gemi_area_profile import profile


def save_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False), encoding='utf-8')


def check_padding(source, current, camera):
    comparisons = []
    with np.load(Path(source)/'image_fields.npz', allow_pickle=False) as saved:
        for i, row in enumerate(current):
            prefix = 'ch{}_'.format(i+1)
            values = {}
            for state, image_state in (('with_particle', 'with'), ('without', 'without')):
                _, _, pixels = integrate_pixels(saved[prefix+'image_axis_um'],
                    saved[prefix+'image_'+image_state+'_w_um2'], camera)
                reference = float(pixels.sum())
                values[state] = abs(row[state]['roi_power_w']/reference-1)
            comparisons.append(values)
    return dict(roi_power_relative_changes_padding4_vs_saved_padding8=comparisons,
                scope='Same 128-sample pupil and camera ROI; checks FFT image sampling only.')


def plot_pairs(path, pairs, titles):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(pairs), 2, figsize=(9, 3.7*len(pairs)), squeeze=False,
                             constrained_layout=True)
    for row, (images, title) in enumerate(zip(pairs, titles)):
        maximum = max(float(image['flux'].max()) for image in images)
        for col, image in enumerate(images):
            axis = image['axis_um']
            selected = abs(axis) <= 50
            im = axes[row, col].pcolormesh(axis[selected], axis[selected],
                image['flux'][np.ix_(selected, selected)], shading='auto', cmap='inferno',
                vmin=0, vmax=maximum)
            axes[row, col].set(title='{} | CH{}'.format(title, col+1), aspect='equal',
                               xlabel='Image x (um)', ylabel='Image y (um)')
        fig.colorbar(im, ax=list(axes[row]), label='W/um^2 (shared scale within each row)')
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_ensemble(output, summary, arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import t
    powers = arrays['powers_w']
    n = np.arange(1, len(powers)+1)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for i, name in enumerate(('aperture', 'camera_roi')):
        delta = powers[:, 0, i]-powers[:, 1, i]
        mean = np.cumsum(delta)/n
        scale = np.cumsum(powers[:, :, i].mean(axis=1))/n
        variance = np.maximum(0, (np.cumsum(delta**2)-n*mean**2)/np.maximum(n-1, 1))
        half = t.ppf(.975, np.maximum(n-1, 1))*np.sqrt(variance/n)
        axes[i].plot(n, 100*mean/scale, color='black')
        axes[i].fill_between(n[1:], 100*(mean-half)[1:]/scale[1:],
                             100*(mean+half)[1:]/scale[1:], alpha=.25, color='teal')
        axes[i].axhline(0, color='gray', linestyle='--')
        axes[i].set(title=name+' | independent surfaces', xlabel='Surface count',
                    ylabel='(CH1 - CH2) / pooled mean (%)')
    fig.suptitle('Fixed {} surfaces | shaded: pointwise paired 95% CI'.format(len(powers)))
    fig.savefig(output/'ensemble_confidence.png', dpi=170)
    plt.close(fig)
    axis = arrays['image_axis_um']
    plot_pairs(output/'haze_ensemble.png',
        [[dict(axis_um=axis, flux=a) for a in arrays[key]] for key in ('first_images', 'mean_images')],
        ['One independent surface', '{}-surface intensity mean'.format(len(powers))])


def run(source, output, samples, seed):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    data, metadata = load_sources(source)
    native = json.loads((output/'native_report.json').read_text(encoding='utf-8'))
    if native['physical_signature'] != metadata['physical_signature']:
        raise ValueError('Independent FEM sources do not match the current source')
    p = Parameters(**metadata['parameters'])
    config = json.loads((Path(source)/'hybrid_input.json').read_text(encoding='utf-8'))
    settings = replace(optical_settings(config['imaging_settings'], p), padding=4)
    if (tuple(settings.numerical_apertures[:2]) != (.5, .5)
            or tuple(settings.tilt_degrees[:2]) != (45, -45) or settings.analyzer != 'none'
            or any(settings.aberrations_nm) or p.polarization != 'p'):
        raise ValueError('This benchmark requires the symmetric .5 NA, +/-45 deg, p-polarized case')
    camera = DetectorParameters(**profile()['detector'])
    print('Checking smooth-wafer symmetry...', flush=True)
    with np.load(output/'flat_sources.npz', allow_pickle=False) as saved:
        flat, flat_images = smooth_check(saved, p, settings, camera)
    save_json(output/'flat_report.json', flat)
    plot_pairs(output/'flat_symmetry.png', list(flat_images.values()),
               ['Plane-wave BV (separate reference)', 'Current Gaussian FEM'])
    del flat_images
    axis = data['whole_axis_um']
    background = surface_fields_order(axis, axis, p,
        metadata['plan'].get('surface_beam_order', p.beam_order))[0]
    background_error = relative_error(background,
        reflect_centered_grid(background)*[-1, 1, 1])
    with np.load(output/'mirror_sources.npz', allow_pickle=False) as fresh:
        mirror, original_images, mirror_images = mirror_check(data, fresh, background, p, settings, camera)
    save_json(output/'mirror_report.json', mirror)
    padding_check = check_padding(source, mirror['original_channels'], camera)
    save_json(output/'padding_report.json', padding_check)
    plot_pairs(output/'rough_mirror.png',
        [[a['with_particle'] for a in original_images], [a['with_particle'] for a in mirror_images]],
        ['Original surface + particle', 'Mirrored surface + fresh local FEM'])
    del original_images, mirror_images
    print('Starting predeclared independent-surface ensemble...', flush=True)
    ensemble, arrays = run_ensemble(data, metadata, background, p, settings, camera, samples, seed)
    np.savez_compressed(output/'ensemble.npz', **arrays)
    save_json(output/'ensemble_report.json', ensemble)
    np.savetxt(output/'ensemble_powers.csv', arrays['powers_w'].reshape(samples, 4), delimiter=',',
        header='ch1_aperture_W,ch1_roi_W,ch2_aperture_W,ch2_roi_W', comments='')
    plot_ensemble(output, ensemble, arrays)
    report = dict(source=str(Path(source).resolve()), physical_signature=metadata['physical_signature'],
        parameters=asdict(p), imaging_settings=asdict(settings), camera_roi=asdict(camera),
        flat=flat, mirror=mirror, ensemble=ensemble,
        illumination_x_reflection_relative_l2=background_error,
        native_report=native,
        padding_check=padding_check,
        sampling_scope='Current pupil sampling, FFT padding 4 rather than source padding 8; '
                       '100 x 100 um camera ROI. Powers before QE, PRNU, dark offset or noise. '
                       'No production simulation or UI parameters were modified.')
    save_json(output/'report.json', report)
    print('Audit completed: '+str((output/'report.json').resolve()), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', default='build/field_coverage_20261016/result')
    parser.add_argument('--output', default='build/channel_symmetry_20261016')
    parser.add_argument('--samples', type=int, default=128)
    parser.add_argument('--seed', type=int, default=202610160)
    args = parser.parse_args()
    run(args.source, args.output, args.samples, args.seed)

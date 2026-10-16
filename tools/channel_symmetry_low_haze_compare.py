"""Compare a freshly solved low-RMS control at unchanged illumination and geometry."""
import argparse
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters, integrate_pixels
from fullwave_prototype.hybrid_run import load_sources
from .channel_symmetry import relative_error
from .channel_symmetry_local import brightness, comparison
from .gemi_area_profile import profile


def measure(source, camera):
    channels, pictures = [], []
    with np.load(source/'image_fields.npz') as saved:
        for i in (1, 2):
            prefix = 'ch{}_'.format(i)
            axis = saved[prefix+'image_axis_um']
            row = {}
            for state in ('with', 'without', 'change'):
                flux = saved[prefix+'image_'+state+'_w_um2']
                _, _, pixels = integrate_pixels(axis, flux, camera)
                image = dict(axis_um=axis, flux=flux, pixels=pixels)
                row[state] = brightness(image, camera)
                if state == 'with':
                    pictures.append(image)
            center = np.argmin(abs(axis))
            row['signed_interference_at_center_w_um2'] = float(
                saved[prefix+'image_interference_w_um2'][center, center])
            channels.append(row)
    return channels, pictures


def run(source, low_source, output):
    data, metadata = load_sources(source)
    low, low_metadata = load_sources(low_source)
    if metadata['parameters'] != low_metadata['parameters']:
        raise ValueError('Low-RMS contrast must keep all illumination and FEM parameters fixed')
    before, after = metadata['stationary_roughness'], low_metadata['stationary_roughness']
    if {k: v for k, v in before.items() if k != 'rms_nm'} != {
            k: v for k, v in after.items() if k != 'rms_nm'}:
        raise ValueError('Only roughness RMS may change in this control')
    np.testing.assert_allclose(data['whole_axis_um'], low['whole_axis_um'], rtol=1e-12, atol=1e-12)
    scale = after['rms_nm']/before['rms_nm']
    height_error = relative_error(low['whole_height_um'], scale*data['whole_height_um'])
    if height_error > 1e-10:
        raise ValueError('The low-RMS surface is not the same scaled realization')
    del data, low
    camera = DetectorParameters(**profile()['detector'])
    original, images = measure(source, camera)
    control, control_images = measure(low_source, camera)
    report = dict(original_rms_nm=before['rms_nm'], low_rms_nm=after['rms_nm'],
        height_scale=scale, same_surface_height_relative_l2=height_error,
        original_channels=original, low_rms_channels=control,
        low_rms_ch1_vs_ch2={state: {key: comparison(control[0][state][key], control[1][state][key])
            for key in control[0][state]} for state in ('with', 'without', 'change')},
        haze_roi_power_ratio_low_vs_original=[control[i]['without']['camera_roi_power_w']/
            original[i]['without']['camera_roi_power_w'] for i in (0, 1)],
        scope='Fresh paired local FEM and whole-surface SPM at identical physical parameters '
              'except RMS; same synthetic morphology, not a measured wafer or forced symmetry. '
              'No per-channel intensity normalization; optical irradiance before camera noise.')
    (output/'low_haze_comparison.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    all_images = images+control_images
    maximum = max(row['with']['central_peak_w_um2'] for row in original+control)
    for i, (ax, image) in enumerate(zip(axes.flat, all_images)):
        axis = image['axis_um']
        selected = abs(axis) <= 12
        im = ax.pcolormesh(axis[selected], axis[selected], image['flux'][np.ix_(selected, selected)],
                          shading='auto', cmap='inferno', vmin=0, vmax=maximum)
        ax.set(title='RMS {} nm | CH{}'.format(before['rms_nm'] if i < 2 else after['rms_nm'], i % 2+1),
               aspect='equal', xlabel='Image x (um)', ylabel='Image y (um)')
    fig.colorbar(im, ax=list(axes.flat), label='W/um^2 | one absolute scale for all panels')
    fig.suptitle('Fresh local FEM | same surface shape, illumination, particle and optics')
    fig.savefig(output/'low_haze_comparison.png', dpi=160)
    plt.close(fig)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('build/field_coverage_20261016/result'))
    parser.add_argument('--output', type=Path, default=Path('build/channel_symmetry_20261016'))
    parser.add_argument('--low-source', type=Path, default=Path('build/channel_symmetry_20261016/low_haze_02nm'))
    args = parser.parse_args()
    run(args.source, args.low_source, args.output)

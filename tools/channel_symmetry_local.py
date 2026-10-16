"""Compare absolute central brightness separately from total collected power."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from detector_model import DetectorParameters, integrate_pixels
from fullwave_prototype.far_field import radiation_amplitude
from fullwave_prototype.hybrid_radiation import triangle_amplitude
from fullwave_prototype.hybrid_run import load_sources, optical_settings
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.surface_scattering import SurfaceRadiation
from optical_propagation import pupil_geometry
from .channel_symmetry import channel_image, reflect_sources
from .channel_symmetry_audit import save_json
from .gemi_area_profile import profile


def brightness(image, camera):
    axis, flux = image['axis_um'], image['flux']
    selected = abs(axis) <= 10
    center = np.argmin(abs(axis))
    box = replace(camera, pixels_x=1, pixels_y=1, pixel_pitch_um=20, fill_factor=1)
    _, _, pixels = integrate_pixels(axis, flux, box)
    return dict(center_w_um2=float(flux[center, center]),
        central_peak_w_um2=float(flux[np.ix_(selected, selected)].max()),
        central_20um_box_power_w=float(pixels.sum()),
        camera_roi_power_w=float(image['pixels'].sum()))


def comparison(a, b):
    if a <= 0 or b <= 0:
        raise ValueError('Positive intensities or powers required')
    return dict(ratio=a/b, relative_to_second=(a-b)/b,
                symmetric_relative_difference=2*(a-b)/(a+b))


def run(source, output):
    data, metadata = load_sources(source)
    p = Parameters(**metadata['parameters'])
    config = json.loads((Path(source)/'hybrid_input.json').read_text(encoding='utf-8'))
    settings = replace(optical_settings(config['imaging_settings'], p), padding=8)
    camera = DetectorParameters(**profile()['detector'])
    images, metrics = [], []
    geometries = [pupil_geometry(.5, tilt, settings.pupil_samples) for tilt in (45, -45)]
    with np.load(Path(source)/'image_fields.npz') as saved:
        for i, geometry in enumerate(geometries):
            field = saved['ch{}_amplitude_with'.format(i+1)]
            _, image = channel_image(field[geometry[3]], geometry, p, settings, camera)
            images.append(image)
            metrics.append(brightness(image, camera))
    print('Computing central-brightness swap for the independent mirrored solve...', flush=True)
    reflected = reflect_sources(data)
    surface = SurfaceRadiation(data['whole_axis_um'], data['whole_axis_um'],
                               reflected['whole_weighted_fields'], p)
    with np.load(output/'mirror_sources.npz') as fresh:
        for geometry in geometries:
            q = geometry[1][geometry[3]]
            amplitude = surface.amplitude(q)-triangle_amplitude(reflected['local_spm_positions_um'],
                reflected['local_spm_weights'], q, p)+radiation_amplitude(
                    fresh['source_points_with'], fresh['source_polarization_weights_with'], q, p)
            _, image = channel_image(amplitude, geometry, p, settings, camera)
            images.append(image)
            metrics.append(brightness(image, camera))
    report = dict(original_channels=metrics[:2], independent_mirror_channels=metrics[2:],
        original_ch1_vs_ch2={key: comparison(metrics[0][key], metrics[1][key]) for key in metrics[0]},
        mirror_swap=[{key: comparison(metrics[i][key], metrics[3-i][key]) for key in metrics[i]}
                     for i in (0, 1)],
        scope='Same absolute irradiance scales; no per-channel normalization or fitted alignment. '
              'Central peak and fixed 20x20 um image box are distinct from the whole sensor ROI. '
              '128 pupil samples, padding 8; peak intensities remain image-sampling dependent.')
    save_json(output/'central_brightness_report.json', report)
    np.savez_compressed(output/'central_brightness_images.npz',
        axis_um=images[0]['axis_um'], irradiance_w_um2=np.array([i['flux'] for i in images]))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    maximum = max(m['central_peak_w_um2'] for m in metrics)
    for i, (ax, image) in enumerate(zip(axes.flat, images)):
        axis = image['axis_um']
        selected = abs(axis) <= 12
        im = ax.pcolormesh(axis[selected], axis[selected], image['flux'][np.ix_(selected, selected)],
                          shading='auto', cmap='inferno', vmin=0, vmax=maximum)
        ax.set(title='{} | CH{}'.format('Original' if i < 2 else 'Fresh mirrored', i % 2+1),
               xlabel='Image x (um)', ylabel='Image y (um)', aspect='equal')
    fig.colorbar(im, ax=list(axes.flat), label='W/um^2 | one absolute scale for all panels')
    fig.suptitle('Central brightness | mirror pairs: top-left vs bottom-right, top-right vs bottom-left')
    fig.savefig(output/'central_brightness.png', dpi=160)
    plt.close(fig)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('build/field_coverage_20261016/result'))
    parser.add_argument('--output', type=Path, default=Path('build/channel_symmetry_20261016'))
    args = parser.parse_args()
    run(args.source, args.output)

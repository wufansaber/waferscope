"""Absolute optical/camera comparisons of the saved 200 nm and 2 um FEM cases."""
import argparse
import json
from pathlib import Path

import numpy as np


def load(folder):
    report = json.loads((folder/'imaging_report.json').read_text(encoding='utf-8'))
    keys = ['ch3_'+key for key in ('image_axis_um', 'image_with_w_um2', 'detector_pixel_power_w',
                                 'detector_pixel_x_um', 'detector_pixel_y_um', 'detector_gray_first_sample_dn')]
    with np.load(folder/'image_fields.npz') as saved:
        arrays = {key: saved[key] for key in keys}
    return report, arrays


def plot(small, large, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Polygon
    from fullwave_prototype.defects import footprint_nm
    output.mkdir(parents=True, exist_ok=True)
    cases = [load(small), load(large)]
    fig, axes = plt.subplots(2, 4, figsize=(18, 9), constrained_layout=True)
    rows = []
    for row, (report, data) in enumerate(cases):
        p, s = report['scattering_parameters'], report['settings']
        size = p['defect_length_nm']/1000
        outline = footprint_nm(p)/1000
        axes[row, 0].add_patch(Polygon(outline, facecolor='#d78ba6', edgecolor='#812d59'))
        limit = size*.65
        axes[row, 0].set(title='Object | L {:g}, W {:g}, H {:g} um'.format(
            size, p['defect_width_nm']/1000, p['defect_height_nm']/1000),
            xlim=(-limit, limit), ylim=(-limit, limit), xlabel='X (um)', ylabel='Y (um)', aspect='equal')
        axis, image = data['ch3_image_axis_um'], data['ch3_image_with_w_um2']
        limit = max(12., size*s['magnification']*.8)
        selected = abs(axis) <= limit+2*(axis[1]-axis[0])
        im = axes[row, 1].pcolormesh(axis[selected], axis[selected], image[np.ix_(selected, selected)],
            shading='auto', cmap='magma', norm=PowerNorm(.5, vmin=0, vmax=max(float(image.max()), 1e-30)))
        axes[row, 1].set(title='CH3 optical image | 20x, NA 0.5', xlim=(-limit, limit), ylim=(-limit, limit),
                        xlabel='Image x (um)', ylabel='Image y (um)', aspect='equal')
        fig.colorbar(im, ax=axes[row, 1], label='W / um^2 (sqrt colors)', shrink=.75)
        px, py = data['ch3_detector_pixel_x_um'], data['ch3_detector_pixel_y_um']
        power = data['ch3_detector_pixel_power_w']
        im = axes[row, 2].pcolormesh(px, py, power, shading='auto', cmap='gray',
            norm=PowerNorm(.5, vmin=0, vmax=max(float(power.max()), 1e-30)))
        axes[row, 2].set(title='Expected 5 um-pixel power', xlim=(-50, 50), ylim=(-50, 50),
                        xlabel='Sensor x (um)', ylabel='Sensor y (um)', aspect='equal')
        fig.colorbar(im, ax=axes[row, 2], label='W / pixel (sqrt colors)', shrink=.75)
        raw = data['ch3_detector_gray_first_sample_dn']
        im = axes[row, 3].pcolormesh(px, py, raw, shading='auto', cmap='gray', vmin=60,
            vmax=max(65, float(raw.max())))
        axes[row, 3].set(title='One simulated area-camera frame', xlabel='Sensor x (um)',
                        ylabel='Sensor y (um)', aspect='equal')
        fig.colorbar(im, ax=axes[row, 3], label='Raw DN; labeled display window', shrink=.75)
        rows.append(dict(length_um=size, width_um=p['defect_width_nm']/1000,
            image_length_um=size*s['magnification'],
            image_width_um=p['defect_width_nm']/1000*s['magnification'],
            length_pixels=size*s['magnification']/s['detector']['pixel_pitch_um'],
            width_pixels=p['defect_width_nm']/1000*s['magnification']/s['detector']['pixel_pitch_um'],
            power_w=[ch['with']['image_power_w'] for ch in report['channels']],
            gray_snr=[ch['detector']['gray_snr'] for ch in report['channels']],
            gray_valid=[ch['detector']['gray_valid'] for ch in report['channels']],
            preview_mesh_um=dict(air=p['air_mesh_um'], substrate=p['substrate_mesh_um'], particle=p['particle_mesh_um'])))
    fig.suptitle('ACTUAL FEM COMPARISON | 266 nm, p, incidence 75 deg from -Y | M20, NA0.5\n'
        '1 W plane-wave normalization over 2 x 2 um, transmission 1e-6 | '
        'absolute scales shown separately; no energy normalization\n'
        'Large case uses a coarser preview mesh; not a certified convergence result')
    fig.savefig(output/'cross_size_comparison.png', dpi=140)
    plt.close(fig)
    # Keep the full three-channel large-object optical view for tilted-channel inspection.
    result = dict(baseline=str(small.resolve()), large=str(large.resolve()), cases=rows,
                  scope='Actual independent FEM solves, smooth wafer, preview meshes. '
                        'No filled geometric cross was substituted for the diffraction image.')
    (output/'comparison_summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--small', type=Path, default=Path('build/fullwave_ui/20261017_164144_979822'))
    parser.add_argument('--large', type=Path, default=Path('build/cross_2um_20261017/preview_scale3'))
    parser.add_argument('--output', type=Path, default=Path('build/cross_2um_20261017'))
    args = parser.parse_args()
    plot(args.small, args.large, args.output)


if __name__ == '__main__':
    main()

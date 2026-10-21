"""Export signed optical and single-raw-frame sample/reference differences."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def difference_arrays(arrays):
    result = {}
    for i in (1, 2, 3):
        prefix = 'ch{}_'.format(i)
        optical = arrays[prefix+'image_with_w_um2']-arrays[prefix+'image_without_w_um2']
        if not np.isfinite(optical).all():
            raise ValueError('Non-finite optical difference')
        np.testing.assert_allclose(optical, arrays[prefix+'image_difference_w_um2'], atol=1e-25)
        sample = arrays[prefix+'detector_gray_first_sample_dn']
        reference = arrays[prefix+'detector_gray_first_reference_dn']
        if sample.shape != reference.shape or sample.ndim != 2:
            raise ValueError('Raw sample/reference frame shapes differ')
        result[prefix+'optical_difference_w_um2'] = optical
        result[prefix+'image_axis_um'] = arrays[prefix+'image_axis_um']
        # Subtract after promotion: uint16 subtraction would wrap negative changes.
        result[prefix+'raw_difference_dn'] = sample.astype(np.int32)-reference.astype(np.int32)
    return result


def plot_differences(output, arrays, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    data = difference_arrays(arrays)
    settings = report['settings']
    detector = settings['detector']
    lx, ly = (detector['pixel_pitch_um']*detector['pixels_x']/2,
              detector['pixel_pitch_um']*detector['pixels_y']/2)
    diameter = report['scattering_parameters']['diameter_nm']
    metadata = '{:g} nm particle | {:g}x | NA {}'.format(diameter, settings['magnification'],
        '/'.join('{:g}'.format(v) for v in settings['numerical_apertures']))
    for kind, units, filename in (('optical', 'W/um^2', 'optical_sample_minus_reference.png'),
                                  ('raw', 'DN', 'camera_sample_minus_reference.png')):
        views = []
        for i in (1, 2, 3):
            prefix = 'ch{}_'.format(i)
            if kind == 'optical':
                axis = data[prefix+'image_axis_um']
                x, y = abs(axis) <= lx, abs(axis) <= ly
                views.append(data[prefix+'optical_difference_w_um2'][np.ix_(y, x)])
            else:
                views.append(data[prefix+'raw_difference_dn'])
        limit = max(float(np.max(abs(view))) for view in views) or 1.
        norm = Normalize(-limit, limit)
        fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
        fig.subplots_adjust(left=.06, right=.85, top=.80, bottom=.23, wspace=.32)
        for i, (ax, value) in enumerate(zip(axes, views), 1):
            if kind == 'optical':
                axis = data['ch{}_image_axis_um'.format(i)]
                artist = ax.pcolormesh(axis[abs(axis) <= lx], axis[abs(axis) <= ly], value,
                                      cmap='RdBu_r', norm=norm, shading='auto')
            else:
                artist = ax.imshow(value, origin='lower', extent=(-lx, lx, -ly, ly),
                    cmap='RdBu_r', norm=norm, interpolation='nearest')
            ax.set(title='CH{} | Sample - Reference'.format(i), aspect='equal',
                xlabel='{} x (um)'.format('Image' if kind == 'optical' else 'Sensor'),
                ylabel='{} y (um)'.format('Image' if kind == 'optical' else 'Sensor'))
        fig.colorbar(artist, cax=fig.add_axes([.89, .26, .014, .49]), label=units+' (shared signed scale)')
        fig.suptitle(metadata+'\n'+('Optical intensity difference' if kind == 'optical'
                                   else 'Single raw camera frame-pair difference'), fontsize=13)
        note = ('Red: sample brighter | Blue: sample darker | No per-channel normalization\n'
                'Unvalidated FEM preview; mesh/PML convergence and physical detectability not established.')
        if kind == 'raw':
            note += ' Includes camera noise; not frame averaging.'
        fig.text(.06, .055, note, fontsize=9, color='#982626', va='bottom')
        fig.savefig(output/filename, dpi=180)
        plt.close(fig)
    np.savez_compressed(output/'sample_minus_reference.npz', **data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    report = json.loads((args.source/'imaging_report.json').read_text(encoding='utf-8'))
    with np.load(args.source/'image_fields.npz', allow_pickle=False) as saved:
        arrays = {key: saved[key] for key in saved.files}
    plot_differences(args.output, arrays, report)
    metadata = dict(operation='sample - reference', source=str(args.source.resolve()),
        report_sha256=hashlib.sha256((args.source/'imaging_report.json').read_bytes()).hexdigest(),
        raw_pair='First saved sample and first saved reference; no averaging or negative-value clipping',
        optical_units='W/um^2', camera_units='signed DN', per_channel_normalization=False,
        physical_validation='Unvalidated FEM preview, not a calibrated measurement')
    (args.output/'difference_export.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()

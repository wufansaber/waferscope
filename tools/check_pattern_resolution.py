"""Audit a real saved FEM pattern without changing geometry or forcing its image."""
import argparse
import json
from pathlib import Path

import numpy as np

from fullwave_prototype.imaging import ImagingSettings, image_from_amplitude
from fullwave_prototype.imaging_audit import check_saved
from fullwave_prototype.imaging_convergence_audit import dipole_amplitude
from optical_propagation import pupil_geometry


def translation_check():
    local, k, rotation, mask, du = pupil_geometry(.5, 0, 64)
    reference = dipole_amplitude(k[mask], (1, 0, 0))
    reports, images = [], []
    offset = .532
    for m in (20, 40):
        amplitude = np.zeros(mask.shape+(2,), complex)
        amplitude[mask] = reference*np.exp(-2j*np.pi/.266*k[mask, 0]*offset)[:, None]
        settings = ImagingSettings(magnification=m, pupil_samples=64, padding=8)
        result = image_from_amplitude(amplitude, local, k, rotation, du, .266, settings)
        peak = np.unravel_index(np.argmax(result['flux']), result['flux'].shape)
        measured = result['axis_um'][peak[1]]
        if not np.isclose(measured, -m*offset, atol=1e-10):
            raise RuntimeError('Image translation does not match magnification')
        reports.append(dict(magnification=m, source_x_um=offset,
                            predicted_image_x_um=-m*offset, image_peak_x_um=float(measured)))
        images.append(result['flux']*m*m)
    error = float(np.linalg.norm(images[0]-images[1])/np.linalg.norm(images[0]))
    if error > .001:
        raise RuntimeError('Magnification changes the object-space profile unexpectedly')
    return dict(points=reports, object_space_profile_relative_l2=error,
                scope='Translated synthetic dipole tests propagation and scale, not FEM defect accuracy')


def dense_image(saved, settings, wavelength_um, index, xy):
    from fullwave_prototype.direct_imaging import direct_image
    prefix = 'ch{}_'.format(index)
    axis = saved[prefix+'pupil_u_axis']
    u, v = np.meshgrid(axis, axis)
    mask = saved[prefix+'pupil_mask']
    uv = np.column_stack((u[mask], v[mask]))
    weights = np.full(len(uv), (axis[1]-axis[0])**2)
    if any(settings['aberrations_nm']):
        raise ValueError('This pattern audit expects a zero-aberration saved pupil')
    _, _, flux = direct_image(saved[prefix+'amplitude_with'][mask], uv, weights, xy,
        wavelength_um, settings['tilt_degrees'][index-1], settings['magnification'],
        settings['image_defocus_um'], settings['transmission'], settings['analyzer'])
    return flux*settings['irradiance_w_um2']


def finer_pupil(source, settings, output):
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.imaging import calculate_channel
    from fullwave_prototype.far_field import radiation_amplitude
    p = Parameters(**json.loads((source/'report.json').read_text(encoding='utf-8'))['parameters'])
    with np.load(source/'paired_fields.npz') as saved:
        data = {key: saved[key] for key in ('source_points_with', 'source_points_without',
            'source_polarization_weights_with', 'source_polarization_weights_without')}
    config = dict(settings)
    config.pop('detector')
    config.update(pupil_samples=64, padding=8, gray_snr_pairs=0)
    config = ImagingSettings(**config)
    def amplitude(directions, name):
        weight = data['source_polarization_weights_'+name]
        if not np.any(weight):
            return np.zeros((len(directions), 2), complex)
        return radiation_amplitude(data['source_points_'+name], weight, directions, p)
    channel, metrics = calculate_channel(data, p, config, 2, amplitude_provider=amplitude)
    np.savez_compressed(output/'finer_ch3.npz', axis_um=channel['image_axis_um'],
        irradiance_w_um2=channel['image_with_w_um2'], pupil=channel['amplitude_with'])
    return channel, metrics


def plot_comparison(output, report, saved, axis, dense, fine):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from matplotlib.patches import Polygon
    from fullwave_prototype.defects import footprint_nm
    p = report['scattering_parameters']
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6), constrained_layout=True)
    axes[0].add_patch(Polygon(footprint_nm(p), facecolor='#d78ba6', edgecolor='#812d59'))
    edge = max(p['defect_length_nm']*.65, 100)
    axes[0].set(title='Actual defect | wafer XY', xlim=(-edge, edge), ylim=(-edge, edge),
                xlabel='X (nm)', ylabel='Y (nm)', aspect='equal')
    peak = max(dense.max(), fine['image_with_w_um2'].max())
    norm = PowerNorm(.5, vmin=0, vmax=peak)
    axes[1].pcolormesh(axis, axis, dense, shading='auto', cmap='magma', norm=norm)
    axes[1].set(title='CH3 | direct vector integral\n32 pupil, 0.2 um image sampling',
                xlim=(-12, 12), ylim=(-12, 12), xlabel='Image x (um)', ylabel='Image y (um)', aspect='equal')
    im = axes[2].pcolormesh(fine['image_axis_um'], fine['image_axis_um'], fine['image_with_w_um2'],
                           shading='auto', cmap='magma', norm=norm)
    axes[2].set(title='CH3 | same FEM source\n64 pupil, padding 8', xlim=(-12, 12), ylim=(-12, 12),
                xlabel='Image x (um)', ylabel='Image y (um)', aspect='equal')
    fig.colorbar(im, ax=list(axes[1:3]), label='W / um^2 | shared sqrt colors', shrink=.7)
    pixels = saved['ch3_detector_pixel_power_w']
    px = saved['ch3_detector_pixel_x_um']
    py = saved['ch3_detector_pixel_y_um']
    im = axes[3].pcolormesh(px, py, pixels, shading='auto', cmap='gray',
        norm=PowerNorm(.5, vmin=0, vmax=max(float(pixels.max()), 1e-30)))
    axes[3].set(title='CH3 | original 5 um pixels\nExpected power, before camera noise',
                xlim=(-15, 15), ylim=(-15, 15), xlabel='Sensor x (um)', ylabel='Sensor y (um)', aspect='equal')
    fig.colorbar(im, ax=axes[3], label='W / pixel', shrink=.7)
    fig.suptitle('REAL SAVED FEM | L {} nm, W {} nm | {} nm, NA {}, M {}\n'
        'Geometry and physics unchanged; finer numerical sampling does not recover sub-resolution detail'.format(
        p['defect_length_nm'], p['defect_width_nm'], p['wavelength_nm'],
        report['settings']['numerical_apertures'][2], report['settings']['magnification']))
    fig.savefig(output/'pattern_resolution.png', dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    report = json.loads((args.source/'imaging_report.json').read_text(encoding='utf-8'))
    settings, p = report['settings'], report['scattering_parameters']
    wavelength, m = p['wavelength_nm']/1000, settings['magnification']
    pitch = settings['detector']['pixel_pitch_um']
    audit, _ = check_saved(args.source)
    if not all(ch['passed'] for ch in audit['channels']):
        raise RuntimeError('Saved FFT disagrees with independent direct diffraction')
    axis = np.linspace(-12, 12, 121)
    x, y = np.meshgrid(axis, axis)
    with np.load(args.source/'image_fields.npz') as loaded:
        saved = dict(loaded)
    dense = dense_image(saved, settings, wavelength, 3, np.column_stack((x.ravel(), y.ravel()))).reshape(x.shape)
    fine, metric = finer_pupil(Path(report['source_directory']), settings, args.output)
    old = saved['ch3_image_axis_um']
    selected = np.flatnonzero(abs(old) <= 12)
    new = fine['image_axis_um']
    ix = np.rint((old[selected]-new[0])/(new[1]-new[0])).astype(int)
    a = saved['ch3_image_with_w_um2'][np.ix_(selected, selected)]
    b = fine['image_with_w_um2'][np.ix_(ix, ix)]
    result = dict(status='passed', source=str(args.source.resolve()), same_pupil_direct_audit=audit,
        translated_point=translation_check(),
        object_length_nm=p['defect_length_nm'], object_width_nm=p['defect_width_nm'],
        image_geometric_length_um=p['defect_length_nm']/1000*m,
        image_geometric_width_um=p['defect_width_nm']/1000*m,
        length_in_camera_pixels=p['defect_length_nm']/1000*m/pitch,
        width_in_camera_pixels=p['defect_width_nm']/1000*m/pitch,
        object_rayleigh_reference_nm=[.61*p['wavelength_nm']/na if na else None
                                     for na in settings['numerical_apertures']],
        object_pixel_nm=pitch/m*1000,
        pupil_32_vs_64_ch3_profile_relative_l2=float(np.linalg.norm(a-b)/np.linalg.norm(b)),
        pupil_32_vs_64_ch3_power_relative_change=float(metric['with']['image_power_w']/
            report['channels'][2]['with']['image_power_w']-1),
        scope='Same saved FEM source; no new FEM solve, no fitted shape, no guarantee of Maxwell mesh convergence. '
              'Rayleigh value is an incoherent-point scale reference, not an exact coherent dark-field shape cutoff.')
    plot_comparison(args.output, report, saved, axis, dense, fine)
    (args.output/'resolution_report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()

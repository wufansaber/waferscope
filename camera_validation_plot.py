"""Inspectable synthetic camera validation and absolute, unnormalized images."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_flat(output, arrays, report):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].plot(arrays['photons'], arrays['signal_dn'], 'o', ms=3)
    axes[0].set(xlabel='Incident photons / pixel', ylabel='Dark-subtracted DN', title='Sensitivity')
    axes[1].plot(arrays['signal_dn'], arrays['variance_dn2'], 'o', ms=3, label='Pair differences')
    axes[1].plot(arrays['signal_dn'], arrays['fitted_variance_dn2'], label='Linear fit')
    axes[1].set(xlabel='Dark-subtracted DN', ylabel='Temporal variance (DN^2)',
                title='Photon transfer: K={:.5g} DN/e-'.format(report['recovered_gain_dn_e']))
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=.25)
    fig.savefig(output/'photon_transfer.png', dpi=180)
    plt.close(fig)


def example_figures(result):
    optical, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    peak = max(float(c['irradiance_w_um2'].max()) for c in result['channels'])
    for i, channel in enumerate(result['channels']):
        axis = channel['image_axis_um']
        im = axes[i].imshow(channel['irradiance_w_um2'], origin='lower',
            extent=(axis[0], axis[-1], axis[0], axis[-1]), vmin=0, vmax=peak, cmap='inferno')
        axes[i].set(title='CH{} | tilt={} deg'.format(i+1, channel['tilt_deg']),
                    xlabel='Image x (um)', ylabel='Image y (um)', xlim=(-20, 20), ylim=(-20, 20))
    optical.colorbar(im, ax=list(axes), label='Irradiance (W/um^2)')
    optical.suptitle('10 nm particle | 266 nm | NA=0.5 | M=20 | absolute optical image')
    camera, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for row, key in enumerate(('mean_signal_e', 'digital_dn')):
        high = max(1, max(float(c['detector'][key].max()) for c in result['channels']))
        for i, channel in enumerate(result['channels']):
            d = channel['detector']
            im = axes[row, i].imshow(d[key], origin='lower', vmin=0, vmax=high,
                cmap='inferno' if row == 0 else 'gray', interpolation='nearest')
            label = 'Expected photoelectrons' if row == 0 else 'Raw CMOS, total SNR={:.3g}'.format(d['emva_total_roi_snr'])
            axes[row, i].set(title='CH{} | {}'.format(i+1, label), xlabel='Pixel x', ylabel='Pixel y')
        camera.colorbar(im, ax=list(axes[row]), label='e- / exposure' if row == 0 else 'DN / exposure')
    p = result['imaging_parameters'].detector
    camera.suptitle('10 nm particle | {:.4g} ms | {:.3g} um pixels | {} bit | {:.3g} e-/DN | black {:.3g} DN'.format(
        p.exposure_ms, p.pixel_pitch_um, p.bit_depth, p.gain_e_dn, p.offset_dn))
    return optical, camera

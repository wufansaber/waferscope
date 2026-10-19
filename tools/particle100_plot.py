"""Three-channel validation: absolute optics and Monte Carlo SNR, no peak scaling."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_summary(output, result):
    for name, info in result.items():
        source = (output/'hybrid266/pupil128' if name == 'hybrid' else
                  Path(info['source']))
        with np.load(source/'image_fields.npz') as saved:
            fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
            for row, field in enumerate(('with', 'without')):
                maximum = max(float(saved['ch{}_image_{}_w_um2'.format(i, field)].max()) for i in (1, 2, 3))
                for i in range(1, 4):
                    axis = saved['ch{}_image_axis_um'.format(i)]
                    selected = np.flatnonzero(abs(axis) <= 20)
                    image = saved['ch{}_image_{}_w_um2'.format(i, field)][np.ix_(selected, selected)]
                    im = axes[row, i-1].pcolormesh(axis[selected], axis[selected], image,
                        vmin=0, vmax=maximum, cmap='inferno', shading='auto')
                    axes[row, i-1].set(title='Ch{} | {}'.format(i, 'Particle + haze' if row == 0 else 'Haze only'),
                                       xlabel='Image x (um)', ylabel='Image y (um)', aspect='equal')
                fig.colorbar(im, ax=list(axes[row]), label='W/um^2 (absolute, separate row scales)')
            fig.suptitle('100 nm particle | {} nm | NA=0.5 | M=20'.format(info['scattering']['wavelength_nm']))
            fig.savefig(output/(name+'_optical_images.png'), dpi=180)
            plt.close(fig)
    fig, axes = plt.subplots(1, len(result), figsize=(6*len(result), 4), squeeze=False, constrained_layout=True)
    for ax, (name, info) in zip(axes[0], result.items()):
        checks = info['checks']
        x = np.arange(3)
        ax.bar(x-.18, [c['monte_carlo']['theoretical_template_snr'] for c in checks], .36, label='Theory')
        ax.bar(x+.18, [c['monte_carlo']['empirical_template_snr'] for c in checks], .36, label='3000-frame MC')
        ax.set(xticks=x, xticklabels=('+45 deg', '-45 deg', '0 deg'), ylabel='Known-template detection SNR',
               title='355 nm plane, 10 us' if name == 'plane' else '266 nm Gaussian, 10 us\nExtra transmission {:.3g}'.format(info['extra_transmission']))
        ax.legend()
    fig.savefig(output/'snr_validation.png', dpi=180)
    plt.close(fig)

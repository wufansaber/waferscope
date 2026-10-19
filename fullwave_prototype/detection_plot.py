"""Absolute paired defect/haze/noise and signed pixel detection SNR."""
import numpy as np
import matplotlib.pyplot as plt


def plot_detection(output, arrays):
    if 'ch1_detector_defect_contrast_e' not in arrays:
        return
    keys = ('defect_contrast_e', 'haze_signal_e', 'detection_noise_rms_e')
    labels = ('Defect contrast (e-)', 'Haze photons converted (e-)', 'Detection noise RMS (e-)')
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), constrained_layout=True)
    for row, (key, label) in enumerate(zip(keys, labels)):
        values = [arrays['ch{}_detector_{}'.format(i, key)] for i in (1, 2, 3)]
        high = max(1e-12, max(float(abs(value).max()) for value in values))
        for col, value in enumerate(values):
            im = axes[row, col].imshow(value, origin='lower', interpolation='nearest',
                vmin=-high if row == 0 else 0, vmax=high, cmap='RdBu_r' if row == 0 else 'viridis')
            axes[row, col].set(title='Ch{} | {}'.format(col+1, label), xlabel='Pixel x', ylabel='Pixel y')
        fig.colorbar(im, ax=list(axes[row]), label=label)
    fig.savefig(output/'detector_components.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    values = [np.ma.array(arrays['ch{}_detector_detection_snr'.format(i)],
                         mask=~arrays['ch{}_detector_detection_valid'.format(i)]) for i in (1, 2, 3)]
    high = max(1e-12, max(float(np.max(np.abs(v.compressed()))) if v.count() else 0 for v in values))
    for i, value in enumerate(values):
        cmap = plt.get_cmap('RdBu_r').copy()
        cmap.set_bad('#808080')
        im = axes[i].imshow(value, origin='lower', interpolation='nearest', vmin=-high, vmax=high, cmap=cmap)
        axes[i].set(title='Ch{} | signed pixel SNR\n{} invalid px'.format(i+1, int(value.mask.sum())),
                    xlabel='Pixel x', ylabel='Pixel y')
    fig.colorbar(im, ax=list(axes), label='Signed detection contrast / noise RMS')
    fig.savefig(output/'detector_snr.png', dpi=160)
    plt.close(fig)

"""Absolute electron counts, with no per-channel normalization."""
import matplotlib.pyplot as plt


def plot_detector(output, arrays):
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for row, key in enumerate(('mean_signal_e', 'measured_e')):
        values = [arrays['ch{}_detector_{}'.format(i, key)] for i in (1, 2, 3)]
        low = min(0, min(float(a.min()) for a in values))
        high = max(1, max(float(a.max()) for a in values))
        for i, value in enumerate(values):
            # Pixel indices remain well-defined for a one-pixel ROI.
            im = axes[row, i].imshow(value, origin='lower', vmin=low, vmax=high, cmap='viridis')
            axes[row, i].set(title='Ch{} {}'.format(i+1, key), xlabel='Pixel x', ylabel='Pixel y')
        fig.colorbar(im, ax=list(axes[row]), label='Electrons / exposure')
    fig.savefig(output/'detector.png', dpi=160)
    plt.close(fig)

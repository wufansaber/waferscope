"""Absolute-scale audit plots, with no image-by-image normalization."""
import numpy as np


def plot(output, arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for column in range(3):
        prefix = 'ch{}_'.format(column+1)
        axis = arrays[prefix+'axis_um']
        direct = arrays[prefix+'with_direct']*1e12
        fft = arrays[prefix+'with_fft']*1e12
        half = (axis[1]-axis[0])/2
        extent = (axis[0]-half, axis[-1]+half, axis[0]-half, axis[-1]+half)
        image = axes[0, column].imshow(direct, origin='lower', extent=extent,
            cmap='magma', vmin=0, vmax=max(float(direct.max()), 1e-30))
        fig.colorbar(image, ax=axes[0, column], label='pW / um2')
        axes[0, column].set_title('Ch{}: direct vector integral'.format(column+1))
        mid = len(axis)//2
        axes[1, column].plot(axis, direct[:, mid], label='Direct integral')
        axes[1, column].plot(axis, fft[:, mid], 'o', ms=3, fillstyle='none', label='Production FFT')
        axes[1, column].set(xlabel='Image y (um), x=0', ylabel='pW / um2')
        axes[1, column].legend(fontsize=8)
        axes[0, column].set(xlabel='Image x (um)', ylabel='Image y (um)')
    fig.savefig(output/'rough_independent.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    axis = arrays['flat_fine_axis_um']
    half = (axis[1]-axis[0])/2
    extent = (axis[0]-half, axis[-1]+half, axis[0]-half, axis[-1]+half)
    for row, pol in enumerate(('p', 's')):
        for column in range(3):
            prefix = '{}_ch{}'.format(pol, column+1)
            reference = arrays[prefix+'_fine_direct']*1e12
            image = axes[row, column].imshow(reference, origin='lower', cmap='magma', vmin=0,
                extent=extent)
            fig.colorbar(image, ax=axes[row, column], label='pW / um2')
            axes[row, column].set(title='BV + direct: {} / Ch{}'.format(pol, column+1),
                                  xlabel='Image x (um)', ylabel='Image y (um)')
    fig.savefig(output/'flat_bv_independent.png', dpi=180)
    plt.close(fig)

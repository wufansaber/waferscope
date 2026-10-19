"""Raw DN images and frame statistics; display window never alters saved values."""
import json
import numpy as np


def export_gray_acquisition(output, arrays):
    from PIL import Image
    import matplotlib.pyplot as plt
    if 'ch1_detector_gray_roi_series_dn' not in arrays:
        return
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    pictures = [[arrays['ch{}_detector_gray_first_{}_dn'.format(i, state)] for i in (1, 2, 3)]
                for state in ('sample', 'reference_mean')]
    low = min(float(a.min()) for row in pictures for a in row)
    high = max(low+1, max(float(a.max()) for row in pictures for a in row))
    raw_arrays = {}
    for i in (1, 2, 3):
        prefix = 'ch{}_detector_'.format(i)
        for row, (values, state) in enumerate(zip(pictures, ('Sample', 'Mean reference'))):
            im = axes[row, i-1].imshow(values[i-1], origin='lower', cmap='gray', vmin=low, vmax=high,
                                        interpolation='nearest')
            axes[row, i-1].set(title='CH{} {}'.format(i, state), xlabel='Pixel x', ylabel='Pixel y')
        for key in ('sample', 'reference'):
            raw = arrays[prefix+'gray_first_'+key+'_dn']
            Image.fromarray(raw).save(output/'gray_ch{}_{}.tiff'.format(i, key))
        series = arrays[prefix+'gray_roi_series_dn']
        np.savetxt(output/'gray_ch{}_series.csv'.format(i), series, delimiter=',',
                   header='sample_roi_DN,mean_reference_roi_DN,difference_roi_DN', comments='')
        for key in ('gray_sample_raw_dn', 'gray_reference_raw_dn'):
            if prefix+key in arrays:
                raw_arrays['ch{}_{}'.format(i, key)] = arrays[prefix+key]
    fig.colorbar(im, ax=list(axes.flat), label='Raw DN | one shared display window')
    fig.suptitle('SNR acquisition: first sample / reference pair | simulated area camera, no TDI')
    fig.savefig(output/'gray_acquisition.png', dpi=160)
    plt.close(fig)
    if raw_arrays:
        np.savez_compressed(output/'gray_raw_stacks.npz', **raw_arrays)

"""Absolute electron counts, with no per-channel normalization."""
import matplotlib.pyplot as plt
from camera_readout import export_raw_camera


def plot_detector(output, arrays):
    paired = 'ch1_detector_reference_mean_dn' in arrays
    keys = ('digital_dn', 'reference_mean_dn', 'difference_dn') if paired else ('mean_signal_e', 'digital_dn')
    fig, axes = plt.subplots(len(keys), 3, figsize=(12, 3.5*len(keys)), constrained_layout=True)
    for row, key in enumerate(keys):
        values = [arrays['ch{}_detector_{}'.format(i, key)] for i in (1, 2, 3)]
        signed = key == 'difference_dn'
        high = max(1, max(float(abs(a).max()) for a in values))
        low = -high if signed else min(0, min(float(a.min()) for a in values))
        for i, value in enumerate(values):
            # Pixel indices remain well-defined for a one-pixel ROI.
            im = axes[row, i].imshow(value, origin='lower', vmin=low, vmax=high,
                                     cmap='RdBu_r' if signed else 'viridis' if key == 'mean_signal_e' else 'gray',
                                     interpolation='nearest')
            label = {'digital_dn': 'Sample DN', 'reference_mean_dn': 'Haze reference DN',
                     'difference_dn': 'Sample - reference DN', 'mean_signal_e': 'Expected e-'}[key]
            title = 'Ch{} {}'.format(i+1, label)
            if key == 'digital_dn':
                well = int(arrays['ch{}_detector_saturated'.format(i+1)].sum())
                adc = int(arrays['ch{}_detector_adc_saturated'.format(i+1)].sum())
                title += '\nWell / ADC saturated: {} / {} px'.format(well, adc)
            axes[row, i].set(title=title, xlabel='Pixel x', ylabel='Pixel y')
            if key == 'digital_dn' and (well or adc):
                axes[row, i].title.set_color('firebrick')
        fig.colorbar(im, ax=list(axes[row]), label='Electrons / exposure' if key == 'mean_signal_e' else 'DN / exposure')
    fig.savefig(output/'detector.png', dpi=160)
    plt.close(fig)
    export_raw_camera(output, [dict(digital_dn=arrays['ch{}_detector_digital_dn'.format(i)]) for i in (1, 2, 3)])
    from .detection_plot import plot_detection
    plot_detection(output, arrays)
    from .gray_acquisition_plot import export_gray_acquisition
    export_gray_acquisition(output, arrays)

"""Raw ADC counts and signed paired differences, never peak normalized."""
from matplotlib.colors import Normalize, TwoSlopeNorm
from .camera_display import adc_maximum, brightness_rows, count_text, pixel_text


def camera_figure(figure, arrays, report=None, full_adc=False):
    figure.clear()
    rows = brightness_rows(arrays, report)
    maximum = adc_maximum(report)
    values = [row for row in rows if row['key'] != 'difference_dn']
    scale = Normalize(0 if maximum is not None and full_adc else min(row['min_dn'] for row in values),
        maximum if maximum is not None and full_adc else max(row['max_dn'] for row in values)+1e-9)
    limit = max([1.]+[max(abs(row['min_dn']), abs(row['max_dn']))
        for row in rows if row['key'] == 'difference_dn'])
    difference_scale = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    frames = {}
    for index, row in enumerate(rows, 1):
        key = row['key']
        data = arrays['ch{}_{}'.format(row['channel'], key)]
        difference = key == 'difference_dn'
        ax = figure.add_subplot(3, 3, index)
        image = ax.imshow(data, origin='lower', interpolation='nearest',
            cmap='RdBu_r' if difference else 'gray', norm=difference_scale if difference else scale)
        title = 'CH{} {}\nmin / max: {:.6g} / {:.6g} DN | mean: {:.6g} DN'.format(
            row['channel'], row['title'], row['min_dn'], row['max_dn'], row['mean_dn'])
        ax.set_title(title, fontsize=9)
        ax.set(xlabel='Pixel x', ylabel='Pixel y')
        ax.tick_params(labelsize=8)
        warnings = []
        if row['at_adc_max_pixels']:
            warnings.append('ADC MAX '+count_text(row['at_adc_max_pixels'], row['pixels']))
        if row['full_well_pixels']:
            warnings.append('FULL WELL '+count_text(row['full_well_pixels'], row['pixels']))
        if warnings:
            ax.text(.03, .97, '\n'.join(warnings), transform=ax.transAxes, va='top',
                color='#b00020', fontsize=8, bbox=dict(facecolor='white', alpha=.9, edgecolor='none'))
        ax.format_coord = lambda x, y, frame=data: pixel_text(frame, x, y)
        figure.colorbar(image, ax=ax, label='DN', shrink=.65)
        frames[ax] = (row, data)
    figure.suptitle('Periodic area camera | raw DN | no TDI', fontsize=11)
    figure.tight_layout(rect=(0, 0, 1, .97))
    return frames

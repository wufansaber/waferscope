"""Display statistics and pixel lookup on original, unscaled ADC data."""
import numpy as np


FRAMES = (('digital_dn', 'Sample'), ('reference_mean_dn', 'Mean reference'),
          ('difference_dn', 'Sample - reference'))


def adc_maximum(report=None):
    detector = (report or {}).get('settings', {}).get('detector', {})
    bits = detector.get('bit_depth')
    return (1 << int(bits))-1 if bits is not None else None


def brightness_rows(arrays, report=None):
    maximum = adc_maximum(report)
    rows = []
    for channel in (1, 2, 3):
        prefix = 'ch{}_'.format(channel)
        for key, title in FRAMES:
            data = np.asarray(arrays[prefix+key])
            if data.ndim != 2 or not data.size or not np.isfinite(data).all():
                raise ValueError('Camera display requires finite, nonempty 2D frames')
            full_well = arrays.get(prefix+'saturated') if key == 'digital_dn' else None
            if full_well is not None and np.shape(full_well) != data.shape:
                raise ValueError('Full-well mask must match the camera frame')
            rows.append(dict(channel=channel, key=key, title=title, pixels=int(data.size),
                min_dn=float(data.min()), max_dn=float(data.max()), mean_dn=float(data.mean()),
                sum_dn=float(data.sum(dtype=float)), adc_max_dn=maximum,
                at_adc_max_pixels=(int(np.count_nonzero(data == maximum))
                    if maximum is not None and key != 'difference_dn' else None),
                full_well_pixels=(int(np.count_nonzero(full_well)) if full_well is not None else None)))
    return rows


def pixel_value(data, x, y):
    if x is None or y is None or not np.isfinite([x, y]).all():
        return None
    height, width = data.shape
    if not (-.5 <= x < width-.5 and -.5 <= y < height-.5):
        return None
    col, row = int(np.floor(x+.5)), int(np.floor(y+.5))
    return col, row, float(data[row, col])


def pixel_text(data, x, y):
    value = pixel_value(data, x, y)
    return '' if value is None else 'x={}  y={}  DN={:.10g}'.format(*value)


def count_text(count, total):
    return 'N/A' if count is None else '{}/{} ({:.3g}%)'.format(count, total, 100*count/total)

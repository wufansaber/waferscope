import copy

import numpy as np
import pytest

from grating_scene.camera_display import adc_maximum, brightness_rows, pixel_value


@pytest.fixture
def frames():
    arrays = {}
    for channel in (1, 2, 3):
        prefix = 'ch{}_'.format(channel)
        arrays[prefix+'digital_dn'] = np.full((2, 3), 1023, dtype=np.uint16)
        arrays[prefix+'reference_mean_dn'] = np.array([[63.5, 64., 65.], [66., 67., 68.]])
        arrays[prefix+'difference_dn'] = arrays[prefix+'digital_dn'].astype(float)-arrays[prefix+'reference_mean_dn']
        arrays[prefix+'saturated'] = np.ones((2, 3), dtype=bool)
    return arrays


def test_statistics_use_actual_adc_bit_depth_not_storage_dtype(frames):
    report = dict(settings=dict(detector=dict(bit_depth=10)))
    assert adc_maximum(report) == 1023
    rows = brightness_rows(frames, report)
    assert rows[0]['min_dn'] == rows[0]['max_dn'] == rows[0]['mean_dn'] == 1023
    assert rows[0]['sum_dn'] == 6138
    assert rows[0]['at_adc_max_pixels'] == rows[0]['full_well_pixels'] == 6
    assert rows[1]['min_dn'] == 63.5
    assert rows[1]['full_well_pixels'] is None
    assert rows[2]['at_adc_max_pixels'] is rows[2]['full_well_pixels'] is None
    assert rows[2]['sum_dn'] == pytest.approx(6*1023-frames['ch1_reference_mean_dn'].sum())
    assert adc_maximum() is None
    assert brightness_rows(frames)[0]['at_adc_max_pixels'] is None


def test_signed_difference_and_reference_mean_do_not_invent_clipping(frames):
    frames['ch1_difference_dn'][0, 0] = -2.5
    frames['ch1_reference_mean_dn'][:] = 1022.75
    row = brightness_rows(frames, dict(settings=dict(detector=dict(bit_depth=10))))
    assert row[2]['min_dn'] == -2.5
    assert row[1]['at_adc_max_pixels'] == 0
    assert row[1]['full_well_pixels'] is None


def test_exact_pixel_lookup_not_transposed_or_outside_frame():
    data = np.array([[1., 2., 3.], [4., 5., 6.25]])
    assert pixel_value(data, 2., 1.) == (2, 1, 6.25)
    assert pixel_value(data, -.5, -.5) == (0, 0, 1.)
    assert pixel_value(data, .5, .5) == (1, 1, 5.)
    for x, y in ((2.5, 0), (0, 1.5), (-.51, 0), (np.nan, 0), (None, 0)):
        assert pixel_value(data, x, y) is None


def test_plot_displays_numeric_saturation_and_preserves_raw_values(frames):
    from matplotlib.figure import Figure
    from grating_scene.camera_plot import camera_figure
    before = copy.deepcopy(frames)
    figure = Figure(figsize=(12, 11))
    panels = camera_figure(figure, frames, dict(settings=dict(detector=dict(bit_depth=10))), full_adc=True)
    assert len(panels) == 9
    ax = next(iter(panels))
    assert ax.images[0].get_clim() == (0, 1023)
    assert '1023 / 1023 DN' in ax.get_title()
    assert 'ADC MAX 6/6 (100%)' in ax.texts[0].get_text()
    assert 'FULL WELL' in ax.texts[0].get_text()
    assert ax.format_coord(0., 0.) == 'x=0  y=0  DN=1023'
    limits = [ax.images[0].get_clim() for ax, (row, _) in panels.items() if row['key'] == 'difference_dn']
    assert limits[0] == limits[1] == limits[2]
    for key in frames:
        np.testing.assert_array_equal(frames[key], before[key])
    panels = camera_figure(figure, frames, dict(settings=dict(detector=dict(bit_depth=10))))
    assert next(iter(panels)).images[0].get_clim()[0] == 63.5


def test_invalid_display_data_rejected(frames):
    frames['ch1_digital_dn'] = np.array([[np.nan]])
    with pytest.raises(ValueError, match='finite'):
        brightness_rows(frames)

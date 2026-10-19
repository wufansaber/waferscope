"""Physical camera chain, fixed-pattern stability, ADC and lossless export."""
from dataclasses import replace

import numpy as np
import pytest

from camera_readout import digitize, export_raw_camera, sensor_maps
from detector_model import DetectorParameters, calculate_detector, integrate_pixels


def uniform_channel(photons=0, size=21, wavelength=266):
    axis = np.arange(size)-(size-1)/2
    intensity = photons*(6.62607015e-34*299792458/(wavelength*1e-9))/.001
    return dict(image_axis_um=axis, irradiance_w_um2=np.full((size, size), intensity),
                power_w=float(intensity*size*size))


def test_exact_adc_transfer_offset_rounding_and_both_clips():
    p = DetectorParameters(gain_e_dn=2, offset_dn=10, bit_depth=8)
    electrons = np.array([[-30., -20., 0., 1.2, 3.2, 490, 600]])
    r = digitize(electrons, np.maximum(electrons, 0), p)
    np.testing.assert_array_equal(r['digital_dn'], [[0, 0, 10, 11, 12, 255, 255]])
    assert r['digital_dn'].dtype == np.uint16
    assert r['adc_clipped_low_pixels'] == 1
    assert r['adc_clipped_high_pixels'] == 1
    assert r['adc_saturated_pixels'] == 2


def test_adc_can_saturate_before_full_well():
    p = DetectorParameters(enabled=True, pixels_x=4, pixels_y=4, pixel_pitch_um=1,
                           quantum_efficiency=1, full_well_e=30000, gain_e_dn=.1,
                           bit_depth=8, offset_dn=10, read_noise_e=0, dark_current_e_s=0)
    r = calculate_detector(uniform_channel(1000), 266, p)
    assert r['saturated_pixels'] == 0
    assert r['adc_saturated_pixels'] == 16
    assert not r['digital_snr_valid'].any()


def test_fill_factor_is_active_pixel_area_not_peak_scaling():
    p = DetectorParameters(pixel_pitch_um=1.5, pixels_x=4, pixels_y=4, fill_factor=.36)
    image = uniform_channel(1)['irradiance_w_um2']
    _, _, power = integrate_pixels(np.arange(-10, 11), image, p)
    np.testing.assert_allclose(power, image[0, 0]*1.5**2*.36, rtol=1e-13)


def test_fixed_maps_stay_fixed_between_frames_and_change_between_sensors():
    p = DetectorParameters(pixel_pitch_um=1, pixels_x=8, pixels_y=8, prnu_fraction=.02,
                           dsnu_e=2, dark_current_std_e_s=.1, dark_current_e_s=1, fixed_pattern_seed=37)
    a = calculate_detector(uniform_channel(100), 266, p, sensor_index=0)
    b = calculate_detector(uniform_channel(100), 266, replace(p, seed=17, exposure_ms=2), sensor_index=0)
    c = calculate_detector(uniform_channel(100), 266, p, sensor_index=1)
    np.testing.assert_array_equal(a['qe_map'], b['qe_map'])
    np.testing.assert_array_equal(a['dark_current_map_e_s'], b['dark_current_map_e_s'])
    np.testing.assert_array_equal(a['dsnu_offset_e'], b['dsnu_offset_e'])
    np.testing.assert_allclose(b['mean_signal_e'], a['mean_signal_e']*2)
    assert not np.array_equal(a['digital_dn'], b['digital_dn'])
    assert not np.array_equal(a['qe_map'], c['qe_map'])
    assert np.all((c['qe_map'] >= 0) & (c['qe_map'] <= 1))
    assert np.all(c['dark_current_map_e_s'] >= 0)


def test_prnu_statistics_and_dsnu_clipping_are_explicit():
    p = DetectorParameters(quantum_efficiency=.5, prnu_fraction=.02,
                           dark_current_e_s=10, dark_current_std_e_s=1, dsnu_e=2)
    qe, dark, offset = sensor_maps(p, (512, 512))
    assert abs(qe.mean()/.5-1) < .0002
    assert abs(qe.std()/qe.mean()-.02) < .0002
    assert abs(dark.mean()-10) < .01
    assert abs(dark.std()-1) < .01
    assert abs(offset.mean()) < .01
    assert abs(offset.std()-2) < .02


def test_digital_photon_transfer_mean_and_variance():
    p = DetectorParameters(pixel_pitch_um=.05, pixels_x=300, pixels_y=300,
                           quantum_efficiency=.6, read_noise_e=3, gain_e_dn=.5,
                           background_e_s=10000, offset_dn=100)
    r = calculate_detector(uniform_channel(0), 266, p)
    measured = r['digital_dn'].astype(float)
    assert abs(measured.mean()-120) < .16
    assert abs(measured.var()-(19/.5**2+1/12)) < 1.2
    assert r['digital_snr_valid'].all()


def test_zero_na_still_has_black_level_and_dark_noise():
    p = DetectorParameters(pixel_pitch_um=1, pixels_x=4, pixels_y=4, read_noise_e=0,
                           dark_current_e_s=0, background_e_s=0, offset_dn=100)
    r = calculate_detector(uniform_channel(0), 266, p)
    np.testing.assert_array_equal(r['digital_dn'], np.full((4, 4), 100))
    assert r['signal_e'] == 0
    assert r['digital_ideal_roi_snr'] == 0


def test_lossless_tiff_roundtrip(tmp_path):
    from PIL import Image
    data = np.array([[0, 1], [4095, 65535]], dtype=np.uint16)
    export_raw_camera(tmp_path, [dict(digital_dn=data)])
    with Image.open(tmp_path/'camera_ch1.tiff') as image:
        np.testing.assert_array_equal(np.asarray(image), data)


@pytest.mark.parametrize('options', [dict(gain_e_dn=0), dict(offset_dn=-1), dict(offset_dn=256, bit_depth=8),
    dict(bit_depth=17), dict(bit_depth=12.5), dict(fill_factor=0), dict(fill_factor=1.1),
    dict(prnu_fraction=-.1), dict(prnu_fraction=.51), dict(dsnu_e=-1),
    dict(dark_current_std_e_s=-1), dict(fixed_pattern_seed=-1)])
def test_invalid_camera_inputs(options):
    with pytest.raises(ValueError):
        replace(DetectorParameters(), **options).validate()

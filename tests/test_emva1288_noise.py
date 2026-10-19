"""EMVA temporal pair statistics, conversion-gain and total-SNR conventions."""
from dataclasses import replace
import numpy as np

from detector_model import DetectorParameters, calculate_detector
from emva1288_noise import noise_metrics, pair_statistics
from tests.test_camera_readout import uniform_channel


def test_equations_14_and_69_gain_is_dn_per_electron():
    p = DetectorParameters(gain_e_dn=.5, read_noise_e=3, dsnu_e=2, prnu_fraction=.01)
    signal = np.full((3, 3), 100.)
    background = np.full((3, 3), 5.)
    result = noise_metrics(signal, background, p)
    expected_temporal = 100+5+9+.5**2/12
    np.testing.assert_allclose(result['temporal_variance_dn2'], expected_temporal/.5**2)
    np.testing.assert_allclose(result['emva_total_snr'], 100/np.sqrt(expected_temporal+4+1))
    assert result['emva_system_gain_dn_e'] == 2
    assert np.all(result['emva_total_snr'] < result['digital_ideal_snr'])


def test_pair_difference_rejects_static_pattern_and_global_frame_offsets():
    fixed = np.arange(100).reshape(10, 10)
    result = pair_statistics(fixed+100, fixed+150)
    assert result['temporal_variance_dn2'] == 0
    assert result['mean_dn'] == 174.5
    np.testing.assert_allclose(result['spatial_covariance_dn2'], fixed.var())


def test_dsnu_is_not_resampled_temporal_noise_or_dark_current():
    p = DetectorParameters(pixel_pitch_um=.05, pixels_x=256, pixels_y=256,
                           dark_current_e_s=0, read_noise_e=0, dsnu_e=4, gain_e_dn=.5)
    channel = uniform_channel(0)
    first = calculate_detector(channel, 266, p)
    second = calculate_detector(channel, 266, replace(p, seed=25, exposure_ms=100))
    np.testing.assert_array_equal(first['digital_dn'], second['digital_dn'])
    assert first['mean_total_e'].sum() == 0
    stats = pair_statistics(first['digital_dn'], second['digital_dn'])
    assert stats['temporal_variance_dn2'] == 0
    assert abs(np.sqrt(stats['spatial_covariance_dn2'])*p.gain_e_dn-4) < .1


def test_exact_gaussian_photons_to_pixel_electrons():
    from camera_validation import pixel_integral_check
    report = pixel_integral_check()
    assert report['passed']
    assert report['pixel_power_relative_l2'] < .001


def test_gain_and_dsnu_do_not_change_physical_full_well():
    from camera_readout import digitize
    p = DetectorParameters(full_well_e=1000)
    r = digitize(np.array([[800.]]), np.array([[800.]]), p, fixed_offset_e=100)
    assert r['digital_snr_valid'].all()
    assert r['nominal_dn'][0, 0] == 1900

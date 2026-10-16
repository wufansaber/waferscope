"""Signed defect contrast, haze shot noise, paired EMVA noise and clipping."""
from dataclasses import replace
import numpy as np
import pytest

from defect_detection import calculate_detection
from detector_model import DetectorParameters


def channel(electrons, p):
    axis = np.arange(-5, 6, dtype=float)
    photon_energy = 6.62607015e-34*299792458/(266e-9)
    intensity = np.full((11, 11), electrons*photon_energy/(p.exposure_ms*.001*p.quantum_efficiency))
    return dict(image_axis_um=axis, irradiance_w_um2=intensity, power_w=float(intensity.sum()))


def parameters(**kwargs):
    return replace(DetectorParameters(enabled=True, pixel_pitch_um=1, pixels_x=2, pixels_y=2,
        dark_current_e_s=0, read_noise_e=2, offset_dn=100, dsnu_e=0), **kwargs)


def test_analytic_particle_haze_camera_budget():
    p = parameters()
    r = calculate_detection(channel(110, p), channel(100, p), 266, p)
    expected = 110+100+2*4+2*.5**2/12
    np.testing.assert_allclose(r['defect_contrast_e'], 10)
    np.testing.assert_allclose(r['haze_signal_e'], 100)
    np.testing.assert_allclose(r['detection_variance_e2'], expected)
    assert r['detection_roi_snr'] == pytest.approx(40/np.sqrt(4*expected))
    assert r['detection_template_snr'] == pytest.approx(2*10/np.sqrt(expected))


def test_no_defect_is_zero_snr_despite_bright_haze():
    p = parameters()
    r = calculate_detection(channel(100, p), channel(100, p), 266, p)
    assert r['detection_roi_snr'] == 0
    assert r['detection_template_snr'] == 0
    assert r['ideal_roi_snr'] > 0
    assert np.any(r['difference_dn'] != 0)


def test_destructive_interference_retains_negative_contrast():
    p = parameters()
    r = calculate_detection(channel(90, p), channel(100, p), 266, p)
    np.testing.assert_allclose(r['defect_contrast_e'], -10)
    assert np.all(r['detection_snr'] < 0)
    assert r['detection_roi_snr'] > 0


def test_more_haze_reduces_detection_at_fixed_contrast():
    p = parameters()
    a = calculate_detection(channel(20, p), channel(10, p), 266, p)
    b = calculate_detection(channel(110, p), channel(100, p), 266, p)
    assert a['detection_template_snr'] > b['detection_template_snr']


def test_reference_stack_variance_and_haze_mismatch():
    p = parameters(reference_frames=4, haze_mismatch_fraction=.1)
    r = calculate_detection(channel(110, p), channel(100, p), 266, p)
    expected = 110+100/4+(4+.5**2/12)*1.25+100
    np.testing.assert_allclose(r['detection_variance_e2'], expected)
    assert r['reference_frames'] == 4
    assert r['haze_mismatch_roi_variance_e2'] == pytest.approx(400)


def test_fixed_dsnu_cancels_in_reference_subtraction():
    p = parameters(read_noise_e=0, dsnu_e=10)
    # Zero light and zero read noise: the two raw frames share the same offset exactly.
    r = calculate_detection(channel(0, p), channel(0, p), 266, p)
    np.testing.assert_array_equal(r['difference_dn'], 0)
    assert r['defect_contrast_roi_e'] == 0
    assert not r['detection_roi_valid']  # Unresolved ADC noise approximation.


def test_saturated_sample_or_reference_invalidates_whole_roi():
    p = parameters(full_well_e=50)
    r = calculate_detection(channel(110, p), channel(100, p), 266, p)
    assert not r['detection_roi_valid']
    assert r['detection_roi_snr'] is None
    assert r['detection_template_snr'] is None
    assert r['detection_invalid_pixels'] == 4


def test_disabled_optical_channel_keeps_camera_noise():
    p = parameters()
    r = calculate_detection(channel(0, p), channel(0, p), 266, p)
    assert r['detection_template_snr'] == 0
    assert np.any(r['difference_dn'] != 0)


def test_independent_reference_difference_matches_predicted_variance():
    p = parameters(pixel_pitch_um=.02, pixels_x=256, pixels_y=256,
                   prnu_fraction=0, dsnu_e=3, reference_frames=4)
    # Divide by area because channel() generates W/um^2 for 1 um^2 pixels.
    a, b = channel(110/.02**2, p), channel(100/.02**2, p)
    r = calculate_detection(a, b, 266, p)
    measured_e = r['difference_dn']*p.gain_e_dn
    assert abs(measured_e.mean()-10) < .15
    assert measured_e.var() == pytest.approx(r['detection_variance_e2'].mean(), rel=.02)


@pytest.mark.parametrize('options', [dict(reference_frames=0), dict(reference_frames=129),
    dict(reference_frames=True), dict(reference_frames=1.5), dict(haze_mismatch_fraction=-1),
    dict(haze_mismatch_fraction=float('nan')), dict(haze_mismatch_fraction=11)])
def test_detection_inputs_validate(options):
    with pytest.raises(ValueError):
        parameters(**options).validate()

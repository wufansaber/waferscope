import numpy as np
import pytest

from tools.gray_frame_snr import measure_paired_frames


def test_gray_estimator_uses_empirical_raw_dn_not_a_noise_budget():
    sample = np.array([110, 112, 108, 114], dtype=np.uint16)[:, None, None]
    reference = np.full_like(sample, 100)
    result, series = measure_paired_frames(sample, reference)
    differences = np.array([10, 12, 8, 14])
    assert result['mean_contrast_roi_dn'] == differences.mean()
    assert result['noise_roi_std_dn'] == pytest.approx(differences.std(ddof=1))
    assert result['snr'] == pytest.approx(abs(differences.mean())/differences.std(ddof=1))
    assert result['nsr']*result['snr'] == pytest.approx(1)
    np.testing.assert_array_equal(series[:, 2], differences)


def test_negative_uint16_difference_is_preserved_and_roi_is_fixed():
    sample = np.array([90, 92, 88, 94], dtype=np.uint16)[:, None, None]
    reference = np.full_like(sample, 100)
    result, _ = measure_paired_frames(sample, reference)
    assert result['mean_contrast_roi_dn'] < 0
    assert result['snr'] > 0
    a = np.tile(sample, (1, 2, 2))
    b = np.tile(reference, (1, 2, 2))
    mask = np.array([[True, False], [False, False]])
    subset, _ = measure_paired_frames(a, b, mask)
    assert subset['snr'] == result['snr']


def test_single_frame_display_scaled_or_mismatching_stacks_are_rejected():
    a = np.full((4, 2, 2), 100, dtype=np.uint16)
    for sample, reference in ((a[:1], a[:1]), (a, a[:, :, :1]), (a.astype(float), a)):
        with pytest.raises(ValueError):
            measure_paired_frames(sample, reference)
    with pytest.raises(ValueError):
        measure_paired_frames(a, a, np.zeros((2, 2), bool))


def test_adc_clipping_invalidates_the_roi_metric():
    a = np.array([100, 110, 120, 65535], dtype=np.uint16)[:, None, None]
    result, _ = measure_paired_frames(a, np.full_like(a, 100))
    assert not result['valid']
    assert result['snr'] is None
    assert result['adc_rail_values'] == 1


def test_gemi_area_parameters_have_no_tdi_gain_and_explicit_uncalibrated_gain():
    from tools.gemi_area_profile import profile
    values = profile()
    p = values['detector']
    assert p['exposure_ms'] == .01
    assert p['quantum_efficiency'] == .639
    assert p['pixel_pitch_um'] == 5
    assert p['read_noise_e'] == 12.1
    assert p['bit_depth'] == 10
    assert p['gain_e_dn'] == 16
    assert 'NO TDI' in values['provenance']['mode']
    assert any('NOT given' in a for a in values['provenance']['assumptions'])


def test_camera_override_reintegrates_pixels_and_writes_raw_adc_not_display(tmp_path):
    import json
    from detector_model import DetectorParameters
    from tools.gemi_area_profile import profile
    from tools.measure_camera_gray_snr import generate_stacks
    axis = np.arange(-12, 13)*5.
    arrays = {}
    for i in (1, 2, 3):
        arrays['ch{}_image_axis_um'.format(i)] = axis
        arrays['ch{}_image_with_w_um2'.format(i)] = np.full((25, 25), 3e-11)
        arrays['ch{}_image_without_w_um2'.format(i)] = np.full((25, 25), 1e-11)
    np.savez_compressed(tmp_path/'image_fields.npz', **arrays)
    report = dict(settings=dict(detector={}), scattering_parameters=dict(wavelength_nm=266),
                  channels=[dict(without=dict(image_power_w=1e-6), **{'with': dict(image_power_w=3e-6)})]*3)
    (tmp_path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
    p = DetectorParameters(**profile()['detector'])
    effective, status, path = generate_stacks(tmp_path, tmp_path, 4, 5, p)
    assert effective['settings']['detector']['pixel_pitch_um'] == 5
    assert effective['settings']['detector']['bit_depth'] == 10
    assert all(c['full_well_clipped_values'] == 0 for c in status)
    with np.load(path, allow_pickle=False) as saved:
        assert saved['ch1_sample_dn'].shape == (4, 20, 20)
        assert saved['ch1_sample_dn'].dtype == np.uint16
        assert saved['ch1_sample_dn'].max() < 1024

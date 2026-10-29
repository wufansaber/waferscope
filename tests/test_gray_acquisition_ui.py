"""Actual ADC acquisition, resource bounds and UI preset round trips."""
import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from dataclasses import asdict, replace

import numpy as np
import pytest
from PyQt5.QtWidgets import QApplication

from camera_profiles import gemi_area_parameters
from defect_detection import calculate_detection
from gray_acquisition import acquire_gray, validate_acquisition
from tools.gray_frame_snr import measure_paired_frames
from tests.test_defect_detection import channel


@pytest.mark.parametrize('reference_frames', (1, 4))
def test_gray_snr_can_be_recomputed_from_exported_integer_frames(reference_frames):
    p = replace(gemi_area_parameters(), pixel_pitch_um=1, pixels_x=2, pixels_y=2, reference_frames=reference_frames)
    detector = calculate_detection(channel(1200, p), channel(300, p), 266, p, np.random.default_rng(4))
    result = acquire_gray(detector, p, 1000, 42)
    assert result['gray_valid']
    raw = result['gray_sample_raw_dn']
    refs = result['gray_reference_raw_dn']
    assert raw.dtype == np.uint16 and refs.shape == (1000, reference_frames, 2, 2)
    assert raw.max() < 1024
    sample_sum = raw.sum(axis=(1, 2), dtype=float)
    ref_sum = refs.mean(axis=1).sum(axis=(1, 2))
    delta = sample_sum-ref_sum
    assert result['gray_snr'] == pytest.approx(abs(delta.mean())/delta.std(ddof=1))
    np.testing.assert_allclose(result['gray_roi_series_dn'][:, 2], delta)
    if reference_frames == 1:
        measured, _ = measure_paired_frames(raw, refs[:, 0], bit_depth=10)
        assert measured['snr'] == pytest.approx(result['gray_snr'])
    again = acquire_gray(detector, p, 1000, 42)
    np.testing.assert_array_equal(again['gray_sample_raw_dn'], raw)


def test_gray_clipping_and_mismatch_are_never_reported_as_valid_snr():
    p = replace(gemi_area_parameters(), pixels_x=2, pixels_y=2, full_well_e=100)
    d = calculate_detection(channel(1200, p), channel(300, p), 266, p)
    result = acquire_gray(d, p, 100, 42)
    assert not result['gray_valid'] and result['gray_snr'] is None
    assert result['gray_full_well_values'] > 0
    mismatch = acquire_gray(d, replace(p, haze_mismatch_fraction=.1), 100, 42)
    assert mismatch['gray_frames'] == 0 and mismatch['gray_snr'] is None
    with pytest.raises(ValueError, match='too large'):
        validate_acquisition(replace(p, pixels_x=512, pixels_y=512), 5000, 42)
    for pairs in (1, 99, 20001):
        with pytest.raises(ValueError):
            validate_acquisition(p, pairs, 42)


def test_gemi_ui_preset_preserves_small_transmission_and_restores_settings():
    app = QApplication.instance() or QApplication([])
    from fullwave_ui.controls import SimulationControls
    from fullwave_prototype.imaging_options import imaging_settings
    from fullwave_prototype.imaging_run import build_parser, scattering_settings
    controls = SimulationControls()
    controls.preset.setCurrentText('266 nm low-haze camera (10 nm gap)')
    args = build_parser().parse_args(controls.arguments())
    p, _ = scattering_settings(args)
    settings = imaging_settings(args, parameters=p)
    assert settings.transmission == .0001
    assert settings.camera_profile == 'gemi_area' and settings.gray_snr_pairs == 5000
    assert settings.detector == gemi_area_parameters()
    assert args.rms_nm == .2 and p.diameter_nm == 100 and p.wavelength_nm == 266
    controls.load_imaging_settings(asdict(settings), p.incidence_deg)
    assert controls.fields['camera-profile'].currentText() == 'gemi_area'
    assert controls.fields['transmission'].value() == .0001
    controls.fields['detector-gain-e-dn'].setValue(15)
    assert controls.fields['camera-profile'].currentText() == 'custom'
    controls.camera_controls.apply.click()
    assert controls.fields['detector-gain-e-dn'].value() == 16
    controls.close()
    app.processEvents()

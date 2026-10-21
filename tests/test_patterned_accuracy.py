import numpy as np
import pytest

from tools.render_patterned_ch3 import raw_triplet
from tools.patterned_numerical_audit import relative


def test_ch3_signed_difference_preserves_negative_dn_without_normalization():
    sample = np.array([[61, 72]], dtype=np.uint16)
    reference = np.array([[65, 63]], dtype=np.uint16)
    a, b, delta = raw_triplet(dict(detector_gray_first_sample_dn=sample,
                                  detector_gray_first_reference_dn=reference))
    np.testing.assert_array_equal(a, sample)
    np.testing.assert_array_equal(b, reference)
    np.testing.assert_array_equal(delta, [[-4, 9]])
    assert delta.dtype == np.int32


def test_ch3_rejects_different_camera_grids():
    with pytest.raises(ValueError):
        raw_triplet(dict(detector_gray_first_sample_dn=np.ones((3, 3)),
                         detector_gray_first_reference_dn=np.ones((4, 4))))


def test_convergence_error_does_not_fit_intensity_or_phase():
    reference = np.array([1+2j, -3j])
    assert relative(2*reference, reference) == pytest.approx(1)
    assert relative(1j*reference, reference) == pytest.approx(np.sqrt(2))


def test_collected_power_uses_solid_angle_not_projected_area():
    from tools.summarize_patterned_accuracy import power
    amplitude = np.array([[1., 0.], [1., 0.]])
    directions = np.array([[0., 0., 1.], [.6, 0., .8]])
    assert power(amplitude, directions, np.array([.2, .2])) == pytest.approx(.45)


def test_invalid_fields_cannot_be_used_as_convergence_evidence():
    with pytest.raises(ValueError):
        relative(np.array([np.nan]), np.array([1]))


def test_convergence_audit_rejects_changed_material():
    from dataclasses import asdict
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.patterned_wafer import LinePatch
    from tools.summarize_patterned_accuracy import check_same_physics
    original = dict(parameters=asdict(Parameters()), pattern=asdict(LinePatch()), settings={})
    changed = dict(original, parameters=dict(original['parameters'], substrate_k=4.4))
    with pytest.raises(ValueError, match='physical inputs'):
        check_same_physics(original, changed)


def test_ch3_direct_imaging_audit_checks_production_pupil_without_normalizing():
    from fullwave_prototype.parameters import Parameters
    from fullwave_prototype.imaging import ImagingSettings, calculate_channel
    from tools.render_patterned_ch3 import check_image_propagation
    p = Parameters()
    settings = ImagingSettings(numerical_apertures=(.6,)*3, magnification=40,
                               pupil_samples=32, padding=2, transmission=1e-6)
    def provider(k, state):
        amplitude = np.column_stack((.2+.1*k[:, 0]+.03j*k[:, 1], .3j+.2*k[:, 0]))
        return amplitude*(1.1 if state == 'with' else 1)
    arrays, _ = calculate_channel({}, p, settings, 2, amplitude_provider=provider)
    result = check_image_propagation(arrays, p, settings)
    assert result['passed']
    assert max(result['flux_relative_l2'].values()) < 1e-12

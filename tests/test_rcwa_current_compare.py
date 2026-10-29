"""Absolute amplitude and phase must survive continuous-angle comparison."""
import numpy as np
import pytest

from rcwa_prototype.current_compare import compare_arrays


def data():
    return dict(directions=np.array([[[0., 0., 1.], [.6, 0., .8]]]),
                solid_angle_weights=np.array([[.3, .4]]), far_electric_um=np.ones((1, 2, 3), complex))


def test_absolute_amplitude_change_is_not_fitted_away():
    a, b = data(), data()
    b['far_electric_um'] *= 2
    row, = compare_arrays(a, b)
    assert row['relative_collected_power_change'] == pytest.approx(3.)
    assert row['relative_weighted_complex_field_change'] == pytest.approx(1.)


def test_phase_change_is_retained_even_when_power_matches():
    a, b = data(), data()
    b['far_electric_um'] *= 1j
    row, = compare_arrays(a, b)
    assert row['relative_collected_power_change'] == pytest.approx(0.)
    assert row['relative_weighted_complex_field_change'] == pytest.approx(np.sqrt(2))


def test_different_angular_samples_cannot_be_compared():
    a, b = data(), data()
    b['directions'][0, 0] = [.6, 0., .8]
    with pytest.raises(AssertionError):
        compare_arrays(a, b)


def test_nonfinite_weight_is_rejected():
    a, b = data(), data()
    a['solid_angle_weights'][0, 0] = np.nan
    b['solid_angle_weights'][0, 0] = np.nan
    with pytest.raises(ValueError, match='Finite fields'):
        compare_arrays(a, b)


def test_identical_nan_directions_are_not_accepted_as_matching():
    a, b = data(), data()
    a['directions'][0, 0, 0] = np.nan
    b['directions'][0, 0, 0] = np.nan
    with pytest.raises(ValueError, match='Finite unit directions'):
        compare_arrays(a, b)

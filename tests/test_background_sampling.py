"""Vectorized incident/reflected/transmitted fields preserve the old ray phases."""
from dataclasses import replace

import numpy as np
import pytest

from fullwave_prototype.background_sampling import electric_background
from fullwave_prototype.local_spm_validation import comparison, probe_directions
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.planar import background_fields


@pytest.mark.parametrize('polarization', ['s', 'p'])
@pytest.mark.parametrize('finite', [False, True])
def test_batched_electric_background_matches_original_ray_sum(polarization, finite):
    p = Parameters(polarization=polarization, wavelength_nm=266)
    if finite:
        p = replace(p, beam_profile='gaussian', beam_length_um=2, beam_width_um=2, beam_order=31)
    points = np.random.default_rng(7).uniform(-.2, .2, (5, 7, 3))
    expected = background_fields(points, p)[0]
    for block in (1, 8, 512):
        np.testing.assert_allclose(electric_background(points, p, block), expected, rtol=1e-12, atol=1e-13)


def test_background_shape_validation_and_empty_points():
    p = Parameters()
    assert electric_background(np.empty((0, 3)), p).shape == (0, 3)
    for points, block in ((np.zeros((4, 2)), 1), (np.array([np.nan, 0, 1]), 1),
                          (np.zeros(3), 0), (np.zeros(3), 2.5), (np.array(0), 1)):
        with pytest.raises(ValueError):
            electric_background(points, p, block)


def test_limit_comparison_does_not_hide_a_phase_or_gain_error():
    spm = np.array([[1, 2j], [-3j, .5]])
    result = comparison(1.1j*spm, spm)
    assert result['complex_relative_l2'] == pytest.approx(abs(1.1j-1))
    assert result['sampled_intensity_ratio'] == pytest.approx(1.21)
    assert result['scalar_corrected_relative_l2'] < 1e-14
    assert result['best_scalar_gain'] == pytest.approx([0, 1.1])
    np.testing.assert_allclose(np.linalg.norm(probe_directions(), axis=-1), 1, atol=1e-14)

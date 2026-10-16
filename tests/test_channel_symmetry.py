"""Do not manufacture statistical symmetry or mis-reflect origin-centered grids."""
import numpy as np
import pytest
from scipy.stats import t

from tools.channel_symmetry import paired_summary, reflect_centered_grid, reflect_sources


def test_origin_centered_reflection_is_an_involution():
    array = np.arange(32).reshape(4, 8)
    reflected = reflect_centered_grid(array)
    np.testing.assert_array_equal(reflected[:, 4], array[:, 4])
    np.testing.assert_array_equal(reflect_centered_grid(reflected), array)
    np.testing.assert_array_equal(reflected[:, 3], array[:, 5])


def test_vector_source_reflection_changes_x_only():
    values = np.arange(12).reshape(4, 3)+1j
    data = {prefix+name: values for prefix in ('source_points_', 'source_polarization_weights_')
            for name in ('with', 'without')}
    data.update(local_spm_positions_um=values[:, :2], local_spm_weights=values,
                whole_weighted_fields=np.ones((4, 8, 3), complex))
    result = reflect_sources(data)
    np.testing.assert_array_equal(result['source_points_with'], values*[-1, 1, 1])
    np.testing.assert_array_equal(reflect_sources(result)['whole_weighted_fields'], data['whole_weighted_fields'])


def test_half_space_radiation_swaps_mirrored_sources():
    from fullwave_prototype.far_field import radiation_amplitude
    from fullwave_prototype.parameters import Parameters
    from optical_propagation import pupil_geometry
    rng = np.random.default_rng(17)
    points = rng.normal(size=(12, 3))*.02
    weights = rng.normal(size=(12, 3))+1j*rng.normal(size=(12, 3))
    geometry = pupil_geometry(.5, 45, 16)
    q = geometry[1][geometry[3]]
    p = Parameters(wavelength_nm=266)
    left = radiation_amplitude(points, weights, q, p)
    right = radiation_amplitude(points*[-1, 1, 1], weights*[-1, 1, 1], q*[-1, 1, 1], p)
    np.testing.assert_allclose(abs(left)**2, abs(right)**2, rtol=2e-13, atol=1e-15)


def test_pupil_to_image_mirror_and_roi_use_correct_grid_conventions():
    from detector_model import DetectorParameters
    from fullwave_prototype.imaging import ImagingSettings
    from fullwave_prototype.parameters import Parameters
    from optical_propagation import pupil_geometry
    from tools.channel_symmetry import channel_image
    p = Parameters(wavelength_nm=266)
    settings = ImagingSettings(pupil_samples=32, padding=4)
    camera = DetectorParameters(pixel_pitch_um=2, pixels_x=8, pixels_y=8)
    geometry = pupil_geometry(.5, 45, 32)
    q = geometry[1][geometry[3]]
    rng = np.random.default_rng(17)
    values = rng.normal(size=(len(q), 2))+1j*rng.normal(size=(len(q), 2))
    amplitude = np.zeros(geometry[3].shape+(2,), complex)
    amplitude[geometry[3]] = values
    mirrored = amplitude[:, ::-1]*[-1, 1]
    other = pupil_geometry(.5, -45, 32)
    a, left = channel_image(values, geometry, p, settings, camera)
    b, right = channel_image(mirrored[other[3]], other, p, settings, camera)
    assert a == pytest.approx(b)
    np.testing.assert_allclose(left['flux'], reflect_centered_grid(right['flux']), atol=1e-15, rtol=1e-12)
    np.testing.assert_allclose(left['pixels'], right['pixels'][:, ::-1], atol=1e-15, rtol=1e-12)


def test_paired_ci_uses_difference_variance_not_independent_channels():
    common = np.arange(1, 41.)
    noise = np.tile([-1., 1.], 20)*.01
    values = np.column_stack((common+noise, common-noise))
    result = paired_summary(values)
    expected = t.ppf(.975, 39)*(2*noise).std(ddof=1)/np.sqrt(40)
    assert result['paired_difference_ci_w'] == pytest.approx([-expected, expected])
    assert result['zero_in_ci']
    assert result['ci_within_equivalence_margin']


def test_inconclusive_is_not_equivalent_and_bias_is_not_hidden():
    rng = np.random.default_rng(15)
    broad = paired_summary(rng.uniform(.1, 2, (30, 2)))
    assert not broad['ci_within_equivalence_margin']
    biased = paired_summary(np.tile([1., 2.], (30, 1)))
    assert not biased['zero_in_ci']
    assert not biased['ci_within_equivalence_margin']


def test_brightness_comparison_preserves_absolute_ratio():
    from tools.channel_symmetry_local import brightness, comparison
    from detector_model import DetectorParameters
    axis = np.arange(-30., 31.)
    flux = np.ones((len(axis), len(axis)))*3
    image = dict(axis_um=axis, flux=flux, pixels=np.full((2, 2), 9.))
    result = brightness(image, DetectorParameters())
    assert result['central_peak_w_um2'] == 3
    assert result['central_20um_box_power_w'] == 3*20*20
    assert result['camera_roi_power_w'] == 36
    assert comparison(4, 3) == pytest.approx(dict(ratio=4/3, relative_to_second=1/3,
                                               symmetric_relative_difference=2/7))


@pytest.mark.parametrize('values', ([[1, 1]], [[-1, 1], [1, 1]], [[0, 0], [0, 0]],
                                  [[1, float('nan')], [1, 1]]))
def test_invalid_ensembles_rejected(values):
    with pytest.raises(ValueError):
        paired_summary(values)

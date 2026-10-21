import numpy as np
import pytest

from tools.plot_patterned_differences import difference_arrays


def fixture():
    data = {}
    for i in (1, 2, 3):
        p = 'ch{}_'.format(i)
        data[p+'image_axis_um'] = np.array([-1., 1.])
        data[p+'image_with_w_um2'] = np.array([[1., 3.], [0., 2.]])*i
        data[p+'image_without_w_um2'] = np.array([[2., 1.], [0., 2.]])*i
        data[p+'image_difference_w_um2'] = np.array([[-1., 2.], [0., 0.]])*i
        data[p+'detector_gray_first_sample_dn'] = np.array([[0, 100], [150, 65535]], np.uint16)
        data[p+'detector_gray_first_reference_dn'] = np.array([[65535, 150], [100, 0]], np.uint16)
    return data


def test_raw_difference_preserves_negative_values_without_uint_wrap():
    result = difference_arrays(fixture())
    for i in (1, 2, 3):
        delta = result['ch{}_raw_difference_dn'.format(i)]
        assert delta.dtype == np.int32
        np.testing.assert_array_equal(delta, [[-65535, -50], [50, 65535]])


def test_optical_difference_preserves_absolute_channel_ratios():
    result = difference_arrays(fixture())
    np.testing.assert_array_equal(result['ch3_optical_difference_w_um2'], 3*result['ch1_optical_difference_w_um2'])


def test_inconsistent_saved_optical_difference_is_rejected():
    data = fixture()
    data['ch1_image_difference_w_um2'] *= -1
    with pytest.raises(AssertionError):
        difference_arrays(data)

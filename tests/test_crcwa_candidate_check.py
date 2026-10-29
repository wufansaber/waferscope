import numpy as np
import pytest

from tools.crcwa_candidate_check import check_interface


def interface_data():
    theta = np.arcsin(np.arange(1, 10)*700e-9/(1.5*5e-6))
    transmitted_cosine = np.sqrt(1-(1.5/5.)**2*np.sin(theta)**2)
    reflection = ((1.5*np.cos(theta)-5*transmitted_cosine)
                  / (1.5*np.cos(theta)+5*transmitted_cosine))**2
    return np.column_stack((np.rad2deg(theta), reflection, reflection))


def test_upstream_exclusive_angle_range():
    theta, expected, r_error, angle_error = check_interface(interface_data())
    assert len(theta) == 9
    assert expected.shape == (9,)
    assert r_error < 1e-14
    assert angle_error < 1e-14


def test_incomplete_native_run_is_rejected():
    with pytest.raises(ValueError, match='nine'):
        check_interface(interface_data()[:-1])


@pytest.mark.parametrize('column', [0, 1])
def test_incorrect_angles_or_reflectance_rejected(column):
    data = interface_data()
    data[0, column] += .01
    with pytest.raises(RuntimeError, match='Fresnel or angle'):
        check_interface(data)


def test_nonfinite_output_rejected():
    data = interface_data()
    data[0, 1] = np.nan
    with pytest.raises(ValueError):
        check_interface(data)

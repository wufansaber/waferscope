"""Keep the continuous imaging oracle independent and its scope explicit."""
import numpy as np

from fullwave_prototype.direct_imaging import polar_cap
from fullwave_prototype.imaging_convergence_audit import (
    axial_dipole_cap, continuous_cap, dipole_amplitude, dipole_case, uniform_disk_case,
)


def test_exact_axial_collection_including_zero_and_small_na():
    assert axial_dipole_cap(0) == 0
    for na in (1e-6, .01, .5, .96):
        uv, weights, k = polar_cap(na, 0, 96, 64)
        amplitude = dipole_amplitude(k, (0, 0, 1))
        integral = np.sum(np.sum(abs(amplitude)**2, axis=1)*weights/k[:, 2])
        np.testing.assert_allclose(integral, axial_dipole_cap(na), rtol=1e-10)


def test_collected_flux_invariant_under_polarization_basis():
    _, _, k = polar_cap(.96, 45, 12, 48)
    moment = np.array((.3, .2j, 1))
    transverse = moment-k*(k@moment)[:, None]
    np.testing.assert_allclose(np.sum(abs(dipole_amplitude(k))**2, axis=1),
                               np.sum(abs(transverse)**2, axis=1), rtol=1e-13)


def test_discrete_power_conservation_is_not_continuous_accuracy():
    report, fields = dipole_case(.5, 0, 32)
    assert report['parseval_relative_error'] < 1e-12
    assert 1e-4 < report['collection_relative_error'] < .03
    assert report['reference_profile_relative_l2'] < 1e-10
    assert report['profile_relative_l2'] < .03
    assert np.isfinite(fields['fft']).all()


def test_high_image_na_against_bessel_oracle():
    report = uniform_disk_case(.96, 1, .4, 128)
    assert report['reference_relative_l2'] < 1e-10
    assert report['flux_relative_l2'] < .01


def test_clipped_horizon_reference_is_converged_and_mirror_symmetric():
    powers = []
    for tilt in (45, -45):
        uv, weights, k = continuous_cap(.96, tilt, 64, 128)
        assert np.all(k[:, 2] > 0)
        assert np.all(np.sum(uv**2, axis=1) < .96**2)
        powers.append(np.sum(weights/np.sqrt(1-np.sum(uv**2, axis=1))))
    np.testing.assert_allclose(powers[0], powers[1], rtol=1e-13)
    report, _ = dipole_case(.96, 45, 32)
    assert report['reference_collection_relative_error'] < 1e-10
    assert report['reference_profile_relative_l2'] < 1e-10

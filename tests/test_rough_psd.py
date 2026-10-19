"""The finite-patch reference must represent the same height field as FEM."""
from dataclasses import replace
import unittest

import numpy as np

from fullwave_prototype.haze_crosscheck import _surface_transform
from fullwave_prototype.rough_psd import expected_finite_psd, height_basis_transform
from fullwave_prototype.rough_statistics import (abc_amplitude, abc_correlation,
                                                 abc_node_covariance,
                                                 nyquist_variance_fraction)
from fullwave_prototype.rough_surface import (Roughness, _taper_grid, height_mesh,
                                              taper_effective_area)
from fullwave_prototype.parameters import Parameters


class RoughPsdTests(unittest.TestCase):
    def setUp(self):
        self.roughness = Roughness(model='gaussian', grid_points=9, taper_fraction=.125)

    def test_hat_basis_recovers_exact_height_realization_transform(self):
        theta = np.deg2rad([0., 45., 80.])
        phi = np.deg2rad([0., 90., 180.])
        polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
        frequency = np.column_stack((np.sin(polar).ravel()*np.cos(azimuth).ravel(),
            np.sin(polar).ravel()*np.sin(azimuth).ravel()-np.sin(np.deg2rad(75))))/.266
        vertices, triangles = height_mesh(self.roughness)
        predicted = height_basis_transform(self.roughness, frequency) @ vertices[..., 2].ravel()
        actual, _ = _surface_transform(triangles, theta, phi, .266, 75, order=6)
        np.testing.assert_allclose(predicted.reshape(actual.shape), actual, rtol=1e-11,
                                   atol=1e-18)

    def test_expected_psd_is_covariance_average_not_per_sample_rescaling(self):
        r = self.roughness
        theta = np.deg2rad([0., 45.])
        phi = np.deg2rad([0., 90.])
        polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
        frequency = np.column_stack((np.sin(polar).ravel()*np.cos(azimuth).ravel(),
                                     np.sin(polar).ravel()*np.sin(azimuth).ravel()))/.266
        basis = height_basis_transform(r, frequency)
        axis = np.linspace(-r.width_um/2, r.width_um/2, r.grid_points)
        kernel = np.exp(-((axis[:, None]-axis[None, :])/r.correlation_length_um)**2)
        weighted = basis*_taper_grid(r).ravel()[None, :]
        covariance = np.kron(kernel, kernel)*(r.rms_nm/1000)**2
        direct = np.einsum('fi,ij,fj->f', weighted.conj(), covariance, weighted).real
        expected = expected_finite_psd(r, theta, phi, .266, 0.)
        np.testing.assert_allclose(expected.ravel(), direct/taper_effective_area(r),
                                   rtol=1e-12)
        doubled = expected_finite_psd(replace(r, rms_nm=4.), theta, phi, .266, 0.)
        np.testing.assert_allclose(doubled, 4*expected, rtol=1e-12)
        self.assertTrue(np.all(direct >= 0))

    def test_abc_variance_normalization_and_exponential_special_case(self):
        from scipy.integrate import quad
        sigma_nm, length = 2., .2
        for exponent in (3., 4.5):
            amplitude = abc_amplitude(sigma_nm, length, exponent)
            variance, _ = quad(lambda f: 2*np.pi*f*amplitude/
                (1+(length*f)**2)**(exponent/2), 0, np.inf,
                epsabs=1e-13, epsrel=1e-10, limit=200)
            self.assertAlmostEqual(variance/(sigma_nm/1000)**2, 1., places=8)
            self.assertAlmostEqual(float(abc_correlation(0., length, exponent)), 1.)
        distance = np.linspace(0, .3, 7)
        np.testing.assert_allclose(abc_correlation(distance, length, 3.),
                                   np.exp(-2*np.pi*distance/length), rtol=1e-12)
        self.assertLess(nyquist_variance_fraction('abc', .025, .05, 3.), .3)

    def test_abc_covariance_matches_numerical_inverse_hankel_transform(self):
        from scipy.integrate import quad
        from scipy.special import j0
        length, exponent, distance = .2, 4.5, .015
        amplitude = abc_amplitude(2., length, exponent)
        variance = (.002)**2
        integral, _ = quad(lambda f: 2*np.pi*f*amplitude*j0(2*np.pi*f*distance)/
            (1+(length*f)**2)**(exponent/2), 0, np.inf,
            epsabs=1e-13, epsrel=1e-9, limit=300)
        self.assertAlmostEqual(integral/variance,
                               float(abc_correlation(distance, length, exponent)),
                               delta=2e-6)

    def test_abc_finite_psd_matches_direct_nodal_covariance(self):
        r = replace(self.roughness, model='abc', correlation_length_um=.2,
                    psd_exponent=3.5)
        theta = np.deg2rad([0., 40.])
        phi = np.deg2rad([0., 90.])
        polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
        frequency = np.column_stack((np.sin(polar).ravel()*np.cos(azimuth).ravel(),
                                     np.sin(polar).ravel()*np.sin(azimuth).ravel()))/.266
        weighted = height_basis_transform(r, frequency)*_taper_grid(r).ravel()[None, :]
        covariance = abc_node_covariance(r.grid_points, r.width_um,
                                          r.correlation_length_um, r.psd_exponent)
        direct = (r.rms_nm/1000)**2*np.einsum(
            'fi,ij,fj->f', weighted.conj(), covariance, weighted).real
        actual = expected_finite_psd(r, theta, phi, .266, 0.)
        np.testing.assert_allclose(actual.ravel(), direct/taper_effective_area(r),
                                   rtol=1e-11)
        self.assertTrue(np.all(actual >= 0))
        with self.assertRaisesRegex(ValueError, 'shared UI'):
            replace(r, psd_exponent=10.1).validate(Parameters(), particle_enabled=False)


if __name__ == '__main__':
    unittest.main()

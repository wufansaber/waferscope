"""Cone quadrature checked against analytic solid angles and smooth angular fields."""
import unittest

import numpy as np
from scipy.integrate import quad

from fullwave_prototype.aperture_integrals import cap_quadrature, integrate_maps


class ApertureIntegralTests(unittest.TestCase):
    def test_exact_solid_angle_without_horizon_clipping(self):
        for na in (.000001, .4, .5, .96):
            for tilt in (0, 5) if na > .9 else (0, 45, -45):
                rays, weights = cap_quadrature(na, tilt)
                self.assertAlmostEqual(weights.sum(), 2*np.pi*na*na/(1+np.sqrt(1-na*na)), places=12)
                np.testing.assert_allclose(np.linalg.norm(rays, axis=1), 1, atol=1e-14)
                self.assertGreaterEqual(rays[:, 2].min(), 0)

    def test_clipped_cap_agrees_with_independent_one_dimensional_integral(self):
        for tilt in (80, -80):
            na = .96
            s, c = np.sin(np.deg2rad(tilt)), np.cos(np.deg2rad(tilt))
            lower = np.sqrt(1-na*na)
            threshold = c*lower/(s*na)
            angle = np.arccos(threshold)
            def span(beta):
                a = max(s*np.cos(beta), 0)
                return 1-max(lower, a/np.hypot(a, c))
            expected = quad(span, 0, 2*np.pi, points=[angle, 2*np.pi-angle], epsabs=1e-12)[0]
            rays, weights = cap_quadrature(na, tilt, 64)
            self.assertAlmostEqual(weights.sum(), expected, places=10)
            self.assertGreaterEqual(rays[:, 2].min(), -1e-14)

    def test_constant_dsc_is_grid_independent_and_disabled_is_zero(self):
        expected = 3*2*np.pi*(1-np.sqrt(.75))
        for step in (10, 5, 2.5):
            theta = np.deg2rad(np.arange(0, 90+step/2, step))
            phi = np.deg2rad(np.arange(0, 360, step))
            totals = integrate_maps(np.full((2, len(theta), len(phi)), 3.), theta, phi,
                                    (.5, .5, 0), (0, 45, -45))
            np.testing.assert_allclose(totals[:, :2], expected, atol=1e-13)
            np.testing.assert_array_equal(totals[:, 2], 0)

    def test_axisymmetric_linear_field_analytic_cap_moment(self):
        rays, weights = cap_quadrature(.5, 45)
        expected = np.pi*.5**2*np.cos(np.pi/4)
        self.assertAlmostEqual(weights @ rays[:, 2], expected, places=12)

    def test_invalid_configuration(self):
        for na, tilt, order in ((1, 0, 20), (-.1, 0, 20), (.5, 90, 20), (.5, 0, 1)):
            with self.assertRaises(ValueError):
                cap_quadrature(na, tilt, order)

    def test_refined_quadrature_preserves_analytic_solid_angle_and_moment(self):
        for order in (144, 192, 256):
            rays, weights = cap_quadrature(.5, 45, order)
            self.assertAlmostEqual(weights.sum(), 2*np.pi*(1-np.sqrt(.75)), places=12)
            self.assertAlmostEqual(weights @ rays[:, 2], np.pi*.25*np.cos(np.pi/4), places=12)


if __name__ == '__main__':
    unittest.main()

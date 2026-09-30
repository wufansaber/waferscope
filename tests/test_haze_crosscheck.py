"""Geometry and integration checks for the read-only haze cross-check."""
import unittest

import numpy as np

from fullwave_prototype.haze_crosscheck import _hemisphere, _surface_transform


class HazeCrosscheckTests(unittest.TestCase):
    def test_constant_height_square_has_exact_zero_frequency_transform(self):
        height, half = .002, .1
        vertices = np.array(((-half, -half, height), (half, -half, height),
                             (half, half, height), (-half, half, height)))
        triangles = vertices[np.array(((0, 1, 2), (0, 2, 3)))]
        transformed, surface = _surface_transform(triangles, np.array((0.,)),
                                                   np.array((0., np.pi)), .266, 0., 4)
        np.testing.assert_allclose(transformed, height*(2*half)**2, rtol=1e-13)
        self.assertAlmostEqual(surface['area_um2'], .04)
        self.assertAlmostEqual(surface['rms_nm'], 2.)

    def test_unit_dsc_integrates_to_upper_hemisphere_solid_angle(self):
        theta = np.linspace(0, np.pi/2, 1001)
        dsc = np.ones((len(theta), 72))
        self.assertAlmostEqual(_hemisphere(dsc, theta), 2*np.pi, delta=2e-6)


if __name__ == '__main__':
    unittest.main()

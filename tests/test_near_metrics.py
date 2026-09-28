import unittest

import numpy as np

from fullwave_prototype.near_metrics import surface_distance, clearances, peak, comparison
from fullwave_prototype.parameters import Parameters


class NearMetricsTests(unittest.TestCase):
    def setUp(self):
        self.flat = np.array([[[-1, -1, 0], [1, -1, 0], [1, 1, 0]],
                              [[-1, -1, 0], [1, 1, 0], [-1, 1, 0]]], float)

    def test_flat_patch_and_exterior_form_plane(self):
        points = np.array([[0, 0, .2], [.9, .8, -.3], [2, 3, .4], [0, 1, .5]])
        np.testing.assert_allclose(surface_distance(points, self.flat), abs(points[:, 2]))

    def test_tilted_patch_uses_normal_not_vertical_distance(self):
        triangles = self.flat.copy()
        triangles[..., 2] = .2*triangles[..., 0]
        points = np.array([[0, 0, .1], [0, 0, -.1]])
        np.testing.assert_allclose(surface_distance(points, triangles), .1/np.sqrt(1.04))

    def test_clearance_to_sphere(self):
        p = Parameters(diameter_nm=100, gap_nm=10)
        points = np.array([[0, 0, .065], [0, 0, .12]])
        _, distance = clearances(points, p, self.flat)
        np.testing.assert_allclose(distance, [.045, .01])

    def test_peak_mask_does_not_leak_interface_values(self):
        points = np.array([[0, 0, 0], [0, 0, 1]], float)
        base = np.array([[100, 0, 0], [2, 0, 0]], complex)
        refined = np.array([[200, 0, 0], [2, 0, 0]], complex)
        report = comparison(points, base, refined, np.array([False, True]))
        self.assertEqual(report['sampled_peak_relative_change'], 0)
        self.assertEqual(report['complex_e_relative_l2'], 0)
        self.assertEqual(report['base_peak']['point_nm'], [0, 0, 1000])
        self.assertEqual(peak(points, refined)['intensity'], 40000)
        with self.assertRaises(ValueError):
            peak(points, base, np.array([False, False]))


if __name__ == '__main__':
    unittest.main()

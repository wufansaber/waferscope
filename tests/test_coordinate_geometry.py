import unittest
from types import SimpleNamespace

import numpy as np

from scatmech_ui.coordinate_geometry import (
    INCIDENCE_PLANE,
    ILLUMINATION_ROTATION_DEG,
    global_to_scatmech_xy,
    incident_propagation_direction,
    incident_source_direction,
    scatmech_to_global_xy,
)
from scatmech_ui.native_far_field import _format_far_field_result
from scatmech_ui.objective_math import pupil_polarization, sample_far_field


class CoordinateGeometryTests(unittest.TestCase):
    def test_incident_ray_comes_from_negative_ky(self):
        root_half = np.sqrt(0.5)

        self.assertEqual(INCIDENCE_PLANE, "KY-KZ")
        self.assertEqual(ILLUMINATION_ROTATION_DEG, 0.0)
        np.testing.assert_allclose(
            incident_propagation_direction(45.0),
            (0.0, root_half, -root_half),
        )
        np.testing.assert_allclose(
            incident_source_direction(45.0),
            (0.0, -root_half, root_half),
        )

    def test_scatmech_xy_rotation_is_invertible(self):
        global_x, global_y = scatmech_to_global_xy(1.0, 0.0)
        self.assertEqual((global_x, global_y), (-0.0, 1.0))
        local_x, local_y = global_to_scatmech_xy(global_x, global_y)
        self.assertEqual((local_x, local_y), (1.0, 0.0))

    def test_objective_sampling_uses_inverse_azimuth_rotation(self):
        theta_axis = np.array([0.0, 0.5, 1.0])
        phi_axis = np.linspace(0.0, 2.0 * np.pi, 5)
        phi_grid, theta_grid = np.meshgrid(phi_axis, theta_axis)
        azimuth_values = np.array([0.0, 1.0, 2.0, 3.0, 0.0])
        intensity = np.tile(azimuth_values, (theta_axis.size, 1))
        source = {
            "theta": theta_grid,
            "phi": phi_grid,
            "scatmech_phi": phi_grid,
            "intensity": intensity,
            "q": intensity,
            "u": intensity,
        }
        radius = np.sin(0.5)
        kx = radius * np.array([0.0, -1.0, 0.0, 1.0])
        ky = radius * np.array([1.0, 0.0, -1.0, 0.0])
        kz = np.full_like(kx, np.cos(0.5))

        sampled = sample_far_field(source, kx, ky, kz)

        np.testing.assert_allclose(sampled["intensity"], [0.0, 1.0, 2.0, 3.0])

    def test_objective_polarization_uses_global_scattering_basis(self):
        result = pupil_polarization(
            intensity=np.ones(1),
            far_field_per_sr=np.ones(1),
            phi=np.array([0.5 * np.pi]),
            kx=np.zeros(1),
            ky=np.array([np.sin(0.5)]),
            kz=np.array([np.cos(0.5)]),
            tilt_deg=0.0,
            aperture=np.ones(1, dtype=bool),
            q=np.ones(1),
            u=np.zeros(1),
        )

        np.testing.assert_allclose(np.abs(result["pol_x"]), [1.0])
        np.testing.assert_allclose(result["pol_y"], [0.0], atol=1e-12)

    def test_rotated_far_field_grid_remains_a_unit_hemisphere(self):
        theta_axis = np.linspace(0.0, 0.5 * np.pi, 5)
        phi_axis = np.linspace(0.0, 2.0 * np.pi, 9)
        phi_grid, theta_grid = np.meshgrid(phi_axis, theta_axis)
        values = np.ones_like(theta_grid)

        result = _format_far_field_result(
            SimpleNamespace(incidence_deg=75.0),
            "test",
            theta_grid,
            phi_grid,
            values,
            values,
            np.zeros_like(values),
            np.zeros_like(values),
        )

        radius_squared = result["x"] ** 2 + result["y"] ** 2 + result["z"] ** 2
        np.testing.assert_allclose(radius_squared, 1.0, atol=1e-12)
        self.assertAlmostEqual(np.min(result["x"]), np.min(result["y"]))
        self.assertAlmostEqual(np.max(result["x"]), np.max(result["y"]))


if __name__ == "__main__":
    unittest.main()

"""Haze ensemble statistics without running expensive Maxwell solves."""
from dataclasses import replace
import unittest
from pathlib import Path
import tempfile

import numpy as np

from fullwave_prototype.haze_ensemble import summarize
from fullwave_prototype.haze_ensemble_plotting import plot_ensemble
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness, height_mesh, raw_height_grid


class HazeEnsembleTests(unittest.TestCase):
    def test_mean_is_of_intensities_and_area_is_explicit(self):
        theta = np.deg2rad([0, 45, 90])
        samples = [np.ones((3, 4)), np.full((3, 4), 9.)]
        reference = np.full((3, 4), 10.)
        roughness = Roughness(model='gaussian', grid_points=9)
        mean, stderr, metrics = summarize(samples, reference, theta, .04, .08, roughness)
        np.testing.assert_allclose(mean, 10.)
        np.testing.assert_allclose(stderr, 8.)
        self.assertAlmostEqual(metrics['fem_to_nist_hemisphere_ratio'], 1.)
        self.assertAlmostEqual(metrics['angular_shape_cosine'], 1.)
        self.assertAlmostEqual(metrics['patch_width_to_correlation_length'], 4.)

    def test_single_sample_has_no_fabricated_standard_error(self):
        theta = np.deg2rad([0, 45, 90])
        mean, stderr, metrics = summarize([np.ones((3, 4))], None, theta,
                                           .04, .04, Roughness(model='gaussian', grid_points=9))
        np.testing.assert_allclose(mean, 1.)
        self.assertTrue(np.isnan(stderr).all())
        self.assertIsNone(metrics['fem_hemisphere_standard_error_um2'])

    def test_finite_patch_reference_and_four_panel_plot(self):
        theta = np.deg2rad([0, 45, 90])
        phi = np.deg2rad([0, 90, 180, 270])
        sample = np.full((3, 4), 2.)
        finite = np.full((3, 4), 2.)
        infinite = np.full((3, 4), 4.)
        mean, stderr, metrics = summarize([sample, sample], infinite, theta,
            .04, .04, Roughness(model='gaussian', grid_points=9), finite)
        self.assertAlmostEqual(metrics['fem_to_finite_patch_hemisphere_ratio'], 1.)
        self.assertAlmostEqual(metrics['finite_patch_to_nist_hemisphere_ratio'], .5)
        self.assertAlmostEqual(metrics['finite_patch_angular_shape_cosine'], 1.)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'haze_ensemble.png'
            report = dict(metrics=metrics, parameters=dict(incidence_deg=75))
            plot_ensemble(path, theta, phi, mean, stderr, infinite, report, finite)
            for name in ('haze_ensemble.png', 'haze_ensemble_fem.png',
                         'haze_ensemble_finite.png', 'haze_ensemble_nist.png',
                         'haze_ensemble_stderr.png'):
                self.assertTrue((path.parent/name).is_file())

    def test_edge_taper_default_preserves_existing_surface(self):
        roughness = Roughness(model='gaussian', grid_points=9)
        raw = raw_height_grid(roughness)
        old = np.sin(np.linspace(0, np.pi, 9))**2
        vertices, _ = height_mesh(roughness)
        np.testing.assert_allclose(vertices[..., 2], raw*np.outer(old, old)*.002, atol=1e-18)
        narrow, _ = height_mesh(replace(roughness, taper_fraction=.125))
        self.assertGreater(np.sqrt(np.mean(narrow[..., 2]**2)),
                           np.sqrt(np.mean(vertices[..., 2]**2)))
        with self.assertRaisesRegex(ValueError, 'at least one'):
            replace(roughness, taper_fraction=.1).validate(Parameters())

    def test_haze_only_does_not_require_particle_clearance(self):
        roughness = Roughness(model='gaussian', grid_points=9, rms_nm=20,
                              taper_fraction=.125)
        parameters = Parameters(gap_nm=0)
        with self.assertRaisesRegex(ValueError, 'Sphere bottom'):
            roughness.validate(parameters)
        roughness.validate(parameters, particle_enabled=False)


if __name__ == '__main__':
    unittest.main()

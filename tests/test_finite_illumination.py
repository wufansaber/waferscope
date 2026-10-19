"""Coherent beam geometry, transversality and absolute-power invariants."""
from dataclasses import asdict, replace
import unittest

import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.illumination import illumination_rays
from fullwave_prototype.finite_illumination import incident_fields, power_area, patch_coverage
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.haze_ensemble import summarize
from fullwave_prototype.haze_crosscheck import statistical_reference
from fullwave_prototype.rough_surface import Roughness


class FiniteIlluminationTests(unittest.TestCase):
    def beam(self, **values):
        return replace(Parameters(wavelength_nm=266, beam_profile='gaussian',
                       beam_length_um=2, beam_width_um=2, beam_order=21), **values)

    def test_rays_are_transverse_downward_and_coherently_summed(self):
        p = self.beam()
        rays = illumination_rays(p)
        for direction, _, _, e, _, _ in rays:
            self.assertLess(direction[2], 0)
            self.assertAlmostEqual(np.linalg.norm(direction), 1.)
            self.assertAlmostEqual(abs(direction@e), 0., places=12)
        e, _ = incident_fields(np.zeros((1, 3)), p)
        np.testing.assert_allclose(e[0], sum(ray[3] for ray in rays), atol=1e-13)
        self.assertGreater(np.linalg.norm(e[0])**2, 5*sum(np.linalg.norm(ray[3])**2 for ray in rays))

    def test_width_changes_field_not_just_power(self):
        p = self.beam(incidence_deg=0)
        points = np.array([[0, 0, 0], [1, 0, 0]])
        narrow = incident_fields(points, p)[0]
        wide = incident_fields(points, replace(p, beam_width_um=4))[0]
        self.assertGreater(np.linalg.norm(wide[1]), 2*np.linalg.norm(narrow[1]))

    def test_power_parseval_matches_spatial_poynting_integral(self):
        p = self.beam(incidence_deg=0)
        axis = np.linspace(-3, 3, 101)
        x, y = np.meshgrid(axis, axis)
        e, h = incident_fields(np.stack((x, y, np.zeros_like(x)), axis=-1), p)
        flux = -np.real(np.cross(e, h.conj()))[..., 2]
        w = np.ones(len(axis))
        w[[0, -1]] = .5
        integrated = float(np.sum(flux*w[:, None]*w)*(axis[1]-axis[0])**2)
        self.assertAlmostEqual(integrated/power_area(p), 1, delta=.001)
        self.assertLess(patch_coverage(p, .2), .04)

    def test_cli_dimensions_and_power_are_coupled_to_scattering(self):
        args = build_parser().parse_args(['--beam-profile', 'gaussian', '--wavelength-nm', '266',
            '--illumination-power-w', '1', '--illumination-length-um', '2', '--illumination-width-um', '3'])
        p, _ = scattering_settings(args)
        self.assertEqual((p.beam_length_um, p.beam_width_um), (2, 3))
        settings = imaging_settings(args, parameters=p)
        self.assertAlmostEqual(settings.irradiance_w_um2*power_area(p), 1)
        args.illumination_length_um = 4
        with self.assertRaisesRegex(ValueError, 'new FEM solve'):
            imaging_settings(args, parameters=p)

    def test_finite_haze_is_not_area_extrapolated_or_a_plane_wave_reference(self):
        theta, phi = np.deg2rad([0, 45, 90]), np.deg2rad([0, 90, 180, 270])
        mean, _, metrics = summarize(np.ones((2, 3, 4)), None, theta, .01, 4,
                                    Roughness(model='abc'), finite_beam=True)
        np.testing.assert_array_equal(mean, 1)
        self.assertFalse(metrics['area_extrapolation_applied'])
        with self.assertRaisesRegex(ValueError, 'finite coherent-beam'):
            statistical_reference(theta, phi, asdict(self.beam()), asdict(Roughness(model='abc')), .04)

    def test_reimage_defaults_use_saved_beam_dimensions_and_incidence(self):
        p = self.beam(incidence_deg=60, beam_length_um=4, beam_width_um=3)
        settings = imaging_settings(build_parser().parse_args([]), parameters=p)
        self.assertEqual((settings.illumination_length_um, settings.illumination_width_um,
                          settings.illumination_incidence_deg), (4, 3, 60))
        self.assertAlmostEqual(settings.irradiance_w_um2*power_area(p), 1)

    def test_default_grazing_beam_is_stable_over_the_displayed_footprint(self):
        p = self.beam(beam_order=31)
        axis = np.linspace(-3, 3, 41)
        x, y = np.meshgrid(axis, axis)
        points = np.stack((x, y, np.zeros_like(x)), axis=-1)
        field, _ = incident_fields(points, p)
        higher, _ = incident_fields(points, replace(p, beam_order=41))
        self.assertLess(np.linalg.norm(field-higher)/np.linalg.norm(higher), 1e-4)
        flux = -np.real(np.cross(field, incident_fields(points, p)[1].conj()))[..., 2]
        weights = np.ones(len(axis))
        weights[[0, -1]] = .5
        captured = np.sum(flux*weights[:, None]*weights)*(axis[1]-axis[0])**2/power_area(p)
        self.assertGreater(captured, .98)
        self.assertLess(captured, 1.001)

    def test_convergence_plan_retains_the_same_finite_incident_field(self):
        from fullwave_prototype.haze_convergence import default_plan, validate_plan, case_settings
        p = self.beam(beam_order=31)
        rough = Roughness(model='abc', correlation_length_um=.2, psd_exponent=4.5)
        plan = default_plan(p, rough, samples=2)
        plan['scans'] = {'angular_step_deg': [10, 5]}
        validate_plan(plan)
        changed, _, angular, _ = case_settings(plan, 'angular_step_deg', 10)
        self.assertEqual(changed, p)
        self.assertEqual(angular, 10)


if __name__ == '__main__':
    unittest.main()

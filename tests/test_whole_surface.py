"""Whole-beam roughness coverage, absolute PSD and vector scattering checks."""
from dataclasses import asdict, replace
import importlib.util
import unittest

import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness
from fullwave_prototype.finite_illumination import incident_fields
from fullwave_prototype.planar import background_fields
from fullwave_prototype.surface_scattering import radiation_kernel, surface_fields, SurfaceRadiation
from fullwave_prototype.whole_surface import plan_surface, random_surface, window


class WholeSurfaceTests(unittest.TestCase):
    def beam(self):
        return Parameters(wavelength_nm=266, beam_profile='gaussian', beam_length_um=2,
                          beam_width_um=2, beam_order=31)

    def test_factored_grid_matches_coherent_ray_fields(self):
        p = self.beam()
        axis = np.linspace(-2, 2, 9)
        x, y = np.meshgrid(axis, axis)
        points = np.stack((x, y, np.zeros_like(x)), axis=-1)
        total, e, h = surface_fields(axis, axis, p)
        expected_e, expected_h = incident_fields(points, p)
        np.testing.assert_allclose(e, expected_e, atol=1e-13)
        np.testing.assert_allclose(h, expected_h, atol=1e-13)
        np.testing.assert_allclose(total, background_fields(points, p)[0], atol=1e-13)

    def test_power_coverage_expands_geometry_without_area_multiplier(self):
        p, plan = plan_surface(self.beam(), Roughness(model='gaussian', correlation_length_um=.05))
        self.assertGreaterEqual(plan['planning_power_coverage'], .99)
        self.assertGreater(plan['width_um'], 2)
        self.assertGreater(plan['grid_points'], 41)
        self.assertGreaterEqual(plan['nyquist_variance_fraction'], .9)
        self.assertLess(plan['field_quadrature_relative_l2_change'], .001)
        self.assertGreaterEqual(p.beam_order, 41)

    def test_no_silent_coarsening_when_grid_limit_is_exceeded(self):
        with self.assertRaisesRegex(ValueError, 'Do not coarsen'):
            plan_surface(self.beam(), Roughness(model='gaussian', correlation_length_um=.001))
        with self.assertRaisesRegex(ValueError, 'Gaussian footprint'):
            plan_surface(Parameters(), Roughness(model='gaussian'))
        with self.assertRaisesRegex(ValueError, 'RMS / wavelength'):
            plan_surface(self.beam(), Roughness(model='gaussian', rms_nm=30))

    def test_absolute_psd_and_no_per_sample_rms_rescaling(self):
        rough = Roughness(model='gaussian', correlation_length_um=.1, rms_nm=2)
        plan = dict(grid_points=128, spacing_um=.02)
        samples = [random_surface(rough, plan, seed) for seed in range(8)]
        rms = np.array([np.sqrt(np.mean(h*h))*1000 for h in samples])
        self.assertAlmostEqual(rms.mean(), 2, delta=.15)
        self.assertGreater(rms.std(), .01)
        np.testing.assert_allclose(random_surface(replace(rough, rms_nm=4), plan, 0), samples[0]*2)

    def test_fft_phase_and_absolute_transform_match_direct_surface_integral(self):
        p = self.beam()
        axis = (np.arange(64)-32)*.04
        x, y = np.meshgrid(axis, axis)
        field = np.stack((np.exp(-x*x-y*y)*np.exp(3j*y), x*np.exp(-3*(x*x+y*y)),
                          y*np.exp(-3*(x*x+y*y))), axis=-1)*.001
        radiation = SurfaceRadiation(axis, axis, field, p, padding=4)
        q = np.array([[.3, .4, np.sqrt(.75)], [-.4, .1, np.sqrt(.83)]])
        transform = np.array([np.sum(field*np.exp(-1j*p.k0*(x*k[0]+y*k[1]))[..., None],
                                           axis=(0, 1))*.04**2 for k in q])
        expected = np.einsum('nij,nj->ni', radiation_kernel(q, p), transform)
        self.assertLess(np.linalg.norm(radiation.amplitude(q)-expected)/np.linalg.norm(expected), .01)

    @unittest.skipUnless(importlib.util.find_spec('pySCATMECH'), 'NIST environment required')
    def test_unit_psd_kernel_matches_nist_for_both_polarizations(self):
        from fullwave_prototype.haze_crosscheck import _nist_dsc
        theta, phi = np.deg2rad([0, 20, 50, 80]), np.deg2rad([0, 60, 120, 180, 270])
        t, f = np.meshgrid(theta, phi, indexing='ij')
        q = np.stack((np.sin(t)*np.cos(f), np.sin(t)*np.sin(f), np.cos(t)), axis=-1)
        for polarization in ('s', 'p'):
            p = Parameters(polarization=polarization)
            e = background_fields(np.zeros(3), p)[0]
            values = np.sum(abs(radiation_kernel(q, p)@e)**2, axis=-1)
            nist = _nist_dsc(theta, phi, asdict(p), {None: 'Unit_PSD_Function'}, 1)
            np.testing.assert_allclose(values, nist, rtol=1e-11, atol=1e-15)

    def test_taper_does_not_remove_the_surface_core(self):
        axis = np.array([-3., -2., 0, 2., 3.])
        w = window(axis, axis, 6, .1)
        self.assertEqual(w[2, 2], 1)
        np.testing.assert_array_equal(w[1:4, 1:4], 1)
        np.testing.assert_array_equal(w[0], 0)

    def test_abc_sampling_keeps_absolute_band_limited_psd_amplitude(self):
        rough = Roughness(model='abc', correlation_length_um=.2, psd_exponent=4.5, rms_nm=2)
        plan = dict(grid_points=128, spacing_um=.02)
        heights = [random_surface(rough, plan, seed) for seed in range(8)]
        self.assertAlmostEqual(np.mean([np.sqrt(np.mean(h*h))*1000 for h in heights]), 2, delta=.2)
        np.testing.assert_allclose(random_surface(replace(rough, rms_nm=4), plan, 0), heights[0]*2)

    def test_direct_solid_angle_integral_and_disabled_channel(self):
        from fullwave_prototype.surface_observables import integrals
        class ConstantRadiation:
            def amplitude(self, directions):
                return np.ones(np.asarray(directions).shape[:-1]+(2,), complex)
        values = integrals(ConstantRadiation(), (0, .5, .5), (0, 45, -45), 32)
        self.assertAlmostEqual(values[0], 4*np.pi)
        self.assertEqual(values[1], 0)
        self.assertAlmostEqual(values[2], 4*np.pi*(1-np.sqrt(1-.5**2)))
        self.assertAlmostEqual(values[3], values[2])


if __name__ == '__main__':
    unittest.main()

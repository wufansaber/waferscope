"""Geometry and integration checks for the read-only haze cross-check."""
import unittest
import importlib.util

import numpy as np

from fullwave_prototype.haze_crosscheck import (_hemisphere, _nist_dsc,
                                                _surface_transform, abc_reference, _channel_integrals)
from fullwave_prototype.rough_statistics import abc_amplitude


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

    def test_channel_crosscheck_uses_cone_integral_not_grid_mask(self):
        theta = np.deg2rad(np.arange(0, 91, 5))
        phi = np.deg2rad(np.arange(0, 360, 5))
        dsc = np.ones((len(theta), len(phi)))
        result = _channel_integrals(dsc, theta, phi, (.5, .5, 0), (0, 45, 0))
        expected = 2*np.pi*(1-np.sqrt(.75))
        np.testing.assert_allclose(result, [expected, expected, 0], atol=1e-13)

    @unittest.skipUnless(importlib.util.find_spec('pySCATMECH'), 'NIST2 environment required')
    def test_abc_reference_uses_nist_psd_at_global_bragg_frequency(self):
        theta = np.deg2rad([0., 30.])
        phi = np.deg2rad([0., 90.])
        parameters = dict(wavelength_nm=266, substrate_n=5.57, substrate_k=.387,
                          incidence_deg=75, polarization='p')
        roughness = dict(model='abc', rms_nm=2., correlation_length_um=.2,
                         psd_exponent=3.5)
        unit = _nist_dsc(theta, phi, parameters, {None: 'Unit_PSD_Function'}, .04)
        abc = abc_reference(theta, phi, parameters, roughness, .04)
        polar, azimuth = np.meshgrid(theta, phi, indexing='ij')
        frequency = np.hypot(np.sin(polar)*np.cos(azimuth),
                             np.sin(polar)*np.sin(azimuth)-np.sin(np.deg2rad(75)))/.266
        expected = abc_amplitude(2., .2, 3.5)/(
            1+(.2*frequency)**2)**(3.5/2)
        np.testing.assert_allclose(abc/unit, expected, rtol=1e-11)


if __name__ == '__main__':
    unittest.main()

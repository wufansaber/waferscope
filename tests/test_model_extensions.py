"""Independent invariants for finite illumination and bounded radiation evaluation."""
from dataclasses import replace
import unittest
import tempfile
from pathlib import Path
import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.illumination import illumination_rays
from fullwave_prototype.planar import background_fields, background_vectors
from fullwave_prototype.far_field import radiation_amplitude_scalar
from fullwave_prototype.radiation import radiation_amplitude
from fullwave_prototype.imaging import ImagingSettings, calculate_images
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.imaging_run import build_parser
from fullwave_prototype.wavefront import pupil_phase
from optical_propagation import pupil_geometry
from detector_model import DetectorParameters


class ExtensionsTest(unittest.TestCase):
    def test_radiation_matches_scalar(self):
        rng = np.random.default_rng(18)
        points = rng.uniform(-.1, .1, (109, 3))
        weights = rng.normal(size=points.shape)+1j*rng.normal(size=points.shape)
        k = rng.normal(size=(75, 3))
        k[:, 2] = abs(k[:, 2])
        k /= np.linalg.norm(k, axis=1)[:, None]
        k[:2] = (0, 0, 1), (1, 0, 0)
        for p in (Parameters(), Parameters(substrate_n=1, substrate_k=0)):
            if p.eps_substrate == 1:
                k[2:20, 2] *= -1
            expected = radiation_amplitude_scalar(points, weights, k, p)
            for block in (1, 17, 100):
                np.testing.assert_allclose(radiation_amplitude(points, weights, k, p, block), expected,
                                           rtol=2e-13, atol=2e-12)

    def test_beam_maxwell_and_interface(self):
        for pol in ('p', 's'):
            p = Parameters(beam_waist_um=2, polarization=pol)
            rays = illumination_rays(p)
            e, h = np.zeros(3, complex), np.zeros(3, complex)
            for kd, ku, kt, ei, er, et in rays:
                for k, field in ((kd, ei), (ku, er), (kt, et)):
                    self.assertLess(abs(k@field), 1e-14)
                np.testing.assert_allclose((ei+er)[:2], et[:2], atol=1e-14)
                np.testing.assert_allclose((np.cross(kd, ei)+np.cross(ku, er))[:2],
                                           np.cross(kt, et)[:2], atol=1e-14)
                e += ei
                h += np.cross(kd, ei)
            self.assertAlmostEqual(float(np.real(np.cross(e, h.conj()))@background_vectors(p)[0]), 1)

    def test_air_radiation_below_origin(self):
        p = Parameters(substrate_n=1, substrate_k=0)
        points = np.array([[0., 0., -.1]])
        weights = np.array([[1., 0., 0.]])
        k = np.array([[0., 0., -1.]])
        actual = radiation_amplitude(points, weights, k, p)[0]
        expected = p.k0**2/(4*np.pi)*np.exp(-1j*p.k0*.1)*np.array([-1., 0.])
        np.testing.assert_allclose(actual, expected, atol=1e-14)

    def test_beam_limit_and_quadrature(self):
        points = np.random.default_rng(4).uniform(-.1, .1, (100, 3))
        p = Parameters()
        for actual, expected in zip(background_fields(points, replace(p, beam_waist_um=1e6)), background_fields(points, p)):
            np.testing.assert_allclose(actual, expected, atol=2e-12)
        a = background_fields(points, replace(p, beam_waist_um=2, beam_order=5))[0]
        b = background_fields(points, replace(p, beam_waist_um=2, beam_order=7))[0]
        self.assertLess(np.linalg.norm(a-b)/np.linalg.norm(b), 1e-6)
        with self.assertRaises(ValueError):
            illumination_rays(replace(p, beam_waist_um=.1))

    def test_wavefront_analytic(self):
        local, _, _, _, du = pupil_geometry(.5, 0, 32)
        phase = pupil_phase(local, du, .355, (10, 0, 0))
        expected = np.exp(2j*np.pi*10*np.sqrt(6)*(local[..., 0]**2-local[..., 1]**2)/(.5**2*355))
        np.testing.assert_allclose(phase, expected, atol=1e-14)
        np.testing.assert_allclose(abs(phase), 1, atol=1e-14)

    def test_compact_detector_and_aberration(self):
        data = {}
        for name in ('with', 'without'):
            data['source_points_'+name] = np.array([[0., 0., .07]])
            data['source_polarization_weights_'+name] = np.array([[1e-5, 0, 0]])*(name == 'with')
        p = Parameters()
        s = ImagingSettings(pupil_samples=32, padding=2,
            detector=DetectorParameters(enabled=True, pixels_x=4, pixels_y=4))
        full, report = calculate_images(data, p, s)
        compact, report2 = calculate_images(data, p, s, False)
        for key, value in compact.items():
            np.testing.assert_array_equal(value, full[key])
        self.assertNotIn('ch1_image_e_with', compact)
        changed, report3 = calculate_images(data, p, replace(s, aberrations_nm=(30, 20, 10)), False)
        self.assertGreater(np.linalg.norm(changed['ch1_image_with_w_um2']-compact['ch1_image_with_w_um2']), 0)
        for a, b in zip(report['channels'], report3['channels']):
            self.assertAlmostEqual(a['with']['image_power_w']/b['with']['image_power_w'], 1)
        self.assertGreater(report2['channels'][0]['detector']['signal_e'], 0)

    def test_cli_settings(self):
        args = build_parser().parse_args(['--detector', '--coma-nm', '15'])
        settings = imaging_settings(args)
        self.assertTrue(settings.detector.enabled)
        self.assertEqual(settings.aberrations_nm, (0, 15, 0))

    def test_cache_invalidation(self):
        from fullwave_prototype.pupil_cache import PupilCache
        data = {}
        for name in ('with', 'without'):
            data['source_points_'+name] = np.array([[0., 0., .07]])
            data['source_polarization_weights_'+name] = np.array([[1e-5, 0, 0]], complex)
        with tempfile.TemporaryDirectory() as directory:
            p = Parameters()
            c = PupilCache(directory, data, p)
            values = {name: np.ones((32, 32, 2), complex) for name in ('with', 'without')}
            c.save(.5, 45, 32, values)
            self.assertIsNotNone(c.load(.5, 45, 32))
            self.assertIsNone(c.load(.4, 45, 32))
            self.assertIsNone(c.load(.5, 0, 32))
            self.assertIsNone(PupilCache(directory, data, replace(p, wavelength_nm=266)).load(.5, 45, 32))
            data['source_polarization_weights_with'][0, 0] *= 2
            self.assertIsNone(PupilCache(directory, data, p).load(.5, 45, 32))
            c.path(.5, 45, 32).write_bytes(b'broken cache')
            self.assertIsNone(c.load(.5, 45, 32))


if __name__ == '__main__':
    unittest.main()

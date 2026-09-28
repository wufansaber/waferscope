import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from scatmech_ui.detector import DetectorParameters, calculate_detector, integrate_pixels
from scatmech_ui.imaging import ImagingParameters, calculate_imaging, export_imaging
from scatmech_ui.parameters import SimulationParameters, ObjectiveCollectionParameters
from scatmech_ui.surface_model import SCATMECH_IMPORT_ERROR


class DetectorTests(unittest.TestCase):
    def setUp(self):
        self.axis = np.arange(-10, 11, dtype=float)
        self.p = DetectorParameters(enabled=True, pixel_pitch_um=1.3, pixels_x=6, pixels_y=4,
                                    dark_current_e_s=0, read_noise_e=0)

    def channel(self, intensity):
        image = np.full((len(self.axis), len(self.axis)), intensity)
        return dict(image_axis_um=self.axis, irradiance_w_um2=image, power_w=float(image.sum()))

    def test_fractional_pixel_overlap_integrates_uniform_field(self):
        x, y, power = integrate_pixels(self.axis, np.full((21, 21), 3.0), self.p)
        self.assertEqual(power.shape, (4, 6))
        np.testing.assert_allclose(power, 3*1.3**2, atol=1e-13)
        self.assertAlmostEqual(x.mean(), 0)
        self.assertAlmostEqual(y.mean(), 0)

    def test_matching_domain_conserves_nonuniform_power(self):
        image = np.random.default_rng(13).random((21, 21))
        p = replace(self.p, pixel_pitch_um=3, pixels_x=7, pixels_y=7)
        _, _, power = integrate_pixels(self.axis, image, p)
        self.assertAlmostEqual(power.sum(), image.sum(), places=11)

    def test_field_boundary_and_invalid_image(self):
        for p in (replace(self.p, pixels_x=100), replace(self.p, pixels_y=100)):
            with self.assertRaisesRegex(ValueError, 'exceeds'):
                integrate_pixels(self.axis, np.ones((21, 21)), p)
        for image in (np.full((21, 21), -1), np.full((21, 21), np.nan)):
            with self.assertRaises(ValueError):
                integrate_pixels(self.axis, image, self.p)
        with self.assertRaises(ValueError):
            integrate_pixels(self.axis[::-1], np.ones((21, 21)), self.p)

    def test_photon_conversion_exposure_qe_and_wavelength(self):
        c = self.channel(1e-15)
        r = calculate_detector(c, 266, self.p)
        expected = 1e-15*1.3**2*.001/(6.62607015e-34*299792458/(266e-9))
        np.testing.assert_allclose(r['mean_photons'], expected)
        np.testing.assert_allclose(r['mean_signal_e'], expected*.6)
        doubled = calculate_detector(c, 532, replace(self.p, exposure_ms=2, quantum_efficiency=.3))
        np.testing.assert_allclose(doubled['mean_signal_e'], r['mean_signal_e']*2)
        self.assertAlmostEqual(r['ideal_roi_snr'], np.sqrt(r['signal_e']))

    def test_zero_signal_keeps_noise_and_background(self):
        p = replace(self.p, dark_current_e_s=2000, background_e_s=3000, read_noise_e=2)
        r = calculate_detector(self.channel(0), 355, p)
        self.assertEqual(r['signal_e'], 0)
        self.assertEqual(r['ideal_roi_snr'], 0)
        self.assertEqual(r['roi_power_fraction'], 0)
        np.testing.assert_allclose(r['mean_total_e'], 5)
        self.assertTrue(np.any(r['measured_e'] != 0))

    def test_seed_reproducibility_and_channel_independence(self):
        p = replace(self.p, read_noise_e=2)
        c = self.channel(1e-15)
        a = calculate_detector(c, 355, p)
        b = calculate_detector(c, 355, p)
        np.testing.assert_array_equal(a['measured_e'], b['measured_e'])
        streams = np.random.SeedSequence(0).spawn(2)
        a = calculate_detector(c, 355, p, np.random.default_rng(streams[0]))
        b = calculate_detector(c, 355, p, np.random.default_rng(streams[1]))
        self.assertFalse(np.array_equal(a['measured_e'], b['measured_e']))

    def test_poisson_and_read_noise_variance(self):
        p = replace(self.p, pixel_pitch_um=.05, pixels_x=300, pixels_y=300,
                    background_e_s=10000, read_noise_e=3)
        measured = calculate_detector(self.channel(0), 355, p)['measured_e']
        self.assertAlmostEqual(measured.mean(), 10, delta=.08)
        self.assertAlmostEqual(measured.var(), 19, delta=.3)

    def test_saturation_before_read_noise(self):
        p = replace(self.p, full_well_e=100, read_noise_e=3)
        r = calculate_detector(self.channel(1e-10), 355, p)
        self.assertEqual(r['saturated_pixels'], 24)
        self.assertEqual(r['expected_overfull_pixels'], 24)
        self.assertTrue(np.any(r['measured_e'] > 100))
        r = calculate_detector(self.channel(1e-10), 355, replace(p, read_noise_e=0))
        np.testing.assert_array_equal(r['measured_e'], np.full((4, 6), 100))

    def test_validation(self):
        for options in (dict(pixel_pitch_um=0), dict(exposure_ms=-1), dict(quantum_efficiency=2),
                        dict(full_well_e=0), dict(read_noise_e=-1), dict(seed=-1),
                        dict(pixels_x=2.5), dict(pixels_y=513), dict(background_e_s=np.nan)):
            with self.assertRaises(ValueError):
                replace(self.p, **options).validate()


@unittest.skipIf(SCATMECH_IMPORT_ERROR is not None, 'SCATMECH unavailable')
class NativeDetectorTests(unittest.TestCase):
    def test_three_channel_pipeline_and_export(self):
        from matplotlib.figure import Figure
        imaging = ImagingParameters(pupil_samples=32, detector=DetectorParameters(enabled=True))
        objective = ObjectiveCollectionParameters(channel_numerical_apertures=(.4, .4, 0))
        result = calculate_imaging(SimulationParameters(wavelength_nm=266), objective, imaging)
        channels = result['channels']
        self.assertGreater(channels[0]['detector']['signal_e'], 0)
        self.assertEqual(channels[2]['detector']['signal_e'], 0)
        for channel in channels:
            self.assertLessEqual(channel['detector']['roi_power_fraction'], 1+1e-12)
        with tempfile.TemporaryDirectory() as directory:
            export_imaging(result, directory, Figure(), Figure())
            payload = json.loads((Path(directory)/'imaging_parameters.json').read_text())
            self.assertTrue(payload['imaging_parameters']['detector']['enabled'])
            self.assertEqual(payload['channels'][2]['detector']['signal_e'], 0)
            with np.load(str(Path(directory)/'coherent_imaging.npz')) as data:
                np.testing.assert_array_equal(data['ch1_detector_measured_e'], channels[0]['detector']['measured_e'])
            self.assertTrue((Path(directory)/'detector_pixels.png').is_file())


if __name__ == '__main__':
    unittest.main()

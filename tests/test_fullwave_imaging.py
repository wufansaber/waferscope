"""Backend-independent checks for the FEM-to-image adapter."""
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import io

import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.imaging import ImagingSettings, calculate_channel


class FullwaveImagingTests(unittest.TestCase):
    def test_cached_sources_reject_new_solver_controls(self):
        from fullwave_prototype.imaging_run import main
        for option in (['--accuracy', 'refined'], ['--pml-power', '2'], ['--pml-mesh-scale', '.8'],
                       ['--direct-solver', 'sparsecholesky'], ['--diameter-nm', '120'],
                       ['--polarization', 's'], ['--roughness-seed', '9'],
                       ['--pml-strength', '1'], ['--integration-bonus', '4'],
                       ['--pml-order', '2'], ['--pml-inner-fraction', '.5'], ['--local-mesh-scale', '.8'],
                       ['--particle-mesh-um', '.0015'], ['--quadrature-order', '7'],
                       ['--particle-mesh-um', '.0015'], ['--quadrature-order', '7'],
                       ['--threads', '4'], ['--near-diagnostics']):
            with patch('sys.argv', ['imaging_run', '--source', 'not-read']+option), patch('sys.stderr', new_callable=io.StringIO):
                with self.assertRaises(SystemExit) as error:
                    main()
                self.assertEqual(error.exception.code, 2)

    def test_physical_cli_settings_reach_model(self):
        from fullwave_prototype.imaging_run import build_parser, scattering_settings
        args = build_parser().parse_args(['--accuracy', 'refined', '--direct-solver', 'sparsecholesky',
            '--wavelength-nm', '266', '--diameter-nm', '120', '--particle-n', '1.7', '--particle-k', '.01',
            '--substrate-n', '4.5', '--substrate-k', '.2', '--incidence-deg', '65', '--polarization', 's',
            '--gap-nm', '20', '--rms-nm', '1', '--roughness-seed', '8'])
        p, r = scattering_settings(args)
        self.assertEqual((p.wavelength_nm, p.diameter_nm, p.particle_n, p.particle_k), (266, 120, 1.7, .01))
        self.assertEqual((p.substrate_n, p.substrate_k, p.incidence_deg, p.polarization), (4.5, .2, 65, 's'))
        self.assertEqual((p.order, p.direct_solver, p.gap_nm, r.rms_nm, r.seed), (3, 'sparsecholesky', 20, 1, 8))

    def test_thread_settings_are_validated(self):
        from fullwave_prototype.imaging_run import build_parser, scattering_settings
        parser = build_parser()
        for count in (1, 2, 4):
            p, _ = scattering_settings(parser.parse_args(['--threads', str(count)]))
            self.assertEqual(p.threads, count)
        for count in (0, -1):
            with self.assertRaises(ValueError):
                scattering_settings(parser.parse_args(['--threads', str(count)]))

    def test_illumination_footprint_sets_irradiance(self):
        from fullwave_prototype.imaging_run import build_parser
        from fullwave_prototype.imaging_options import imaging_settings
        parser = build_parser()
        default = imaging_settings(parser.parse_args([]))
        self.assertEqual((default.illumination_power_w, default.illumination_length_um,
                          default.illumination_width_um), (1, 2, 2))
        self.assertAlmostEqual(default.irradiance_w_um2,
                               1/(4*np.cos(np.deg2rad(75))))
        oblique = imaging_settings(parser.parse_args(['--incidence-deg', '60']))
        self.assertAlmostEqual(oblique.irradiance_w_um2, .5)
        legacy = imaging_settings(parser.parse_args(['--irradiance-w-um2', '1e-6']))
        self.assertEqual(legacy.irradiance_w_um2, 1e-6)
        self.assertIsNone(legacy.illumination_power_w)
        with self.assertRaisesRegex(ValueError, 'together'):
            imaging_settings(parser.parse_args(['--illumination-power-w', '1']))
        with self.assertRaisesRegex(ValueError, 'conflicts'):
            imaging_settings(parser.parse_args(['--illumination-power-w', '1',
                '--illumination-length-um', '2', '--illumination-width-um', '2',
                '--irradiance-w-um2', '1e-6']))

    def test_reimage_default_uses_saved_incidence(self):
        from fullwave_prototype.imaging_run import main
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)
            (source/'report.json').write_text(json.dumps(dict(parameters=dict(
                incidence_deg=60))), encoding='utf-8')
            with patch('sys.argv', ['imaging_run', '--source', str(source)]), patch(
                    'fullwave_prototype.imaging_run.image_saved_sources') as reimage:
                main()
            self.assertAlmostEqual(reimage.call_args.args[2].irradiance_w_um2, .5)

    def setUp(self):
        self.p = Parameters(substrate_n=1, substrate_k=0)
        self.settings = ImagingSettings(pupil_samples=32, padding=4,
                                       numerical_apertures=(.4, .4, 0), tilt_degrees=(0, 45, -45))
        self.data = dict(source_points_with=np.zeros((1, 3)), source_points_without=np.zeros((1, 3)),
                         source_polarization_weights_with=np.array(((0., 1., 0.),), complex),
                         source_polarization_weights_without=np.array(((0., .2j, 0.),), complex))

    def test_disabled_channel_and_finite_validation(self):
        arrays, metrics = calculate_channel(self.data, self.p, self.settings, 2)
        self.assertEqual(metrics['with']['image_power_w'], 0)
        np.testing.assert_array_equal(arrays['image_with_w_um2'], 0)
        for options in (dict(numerical_apertures=(.4, .4)), dict(numerical_apertures=(1, 0, 0)),
                        dict(tilt_degrees=(0, 90, 0)), dict(focal_length_mm=float('inf')),
                        dict(transmission=-1), dict(pupil_samples=20), dict(magnification=.5)):
            with self.assertRaises(ValueError):
                replace(self.settings, **options).validate()

    def test_pupil_image_power_and_coherent_decomposition(self):
        arrays, metrics = calculate_channel(self.data, self.p, self.settings, 1)
        self.assertLess(metrics['with']['flux_relative_error'], 1e-12)
        self.assertLess(metrics['without']['flux_relative_error'], 1e-12)
        self.assertLess(metrics['decomposition_relative_max_error'], 1e-12)
        np.testing.assert_allclose(arrays['image_with_w_um2'], arrays['image_without_w_um2']+
            arrays['image_change_w_um2']+arrays['image_interference_w_um2'], atol=1e-20)
        self.assertFalse(np.allclose(arrays['image_difference_w_um2'], arrays['image_change_w_um2'], atol=1e-20))

    def test_absolute_irradiance_transmission_and_analyzers(self):
        base, metrics = calculate_channel(self.data, self.p, self.settings, 0)
        scaled, _ = calculate_channel(self.data, self.p,
            replace(self.settings, irradiance_w_um2=3e-6, transmission=.5), 0)
        np.testing.assert_allclose(scaled['image_with_w_um2'], 1.5*base['image_with_w_um2'], atol=1e-20)
        _, x = calculate_channel(self.data, self.p, replace(self.settings, analyzer='x'), 0)
        _, y = calculate_channel(self.data, self.p, replace(self.settings, analyzer='y'), 0)
        self.assertAlmostEqual((x['with']['image_power_w']+y['with']['image_power_w'])/
                               metrics['with']['image_power_w'], 1, places=12)

    def test_magnification_and_image_defocus_preserve_power(self):
        base, metrics = calculate_channel(self.data, self.p, self.settings, 0)
        larger, lm = calculate_channel(self.data, self.p, replace(self.settings, magnification=40), 0)
        np.testing.assert_allclose(larger['image_axis_um'], 2*base['image_axis_um'])
        self.assertAlmostEqual(lm['with']['image_power_w']/metrics['with']['image_power_w'], 1, places=12)
        defocused, dm = calculate_channel(self.data, self.p, replace(self.settings, image_defocus_um=400), 0)
        self.assertAlmostEqual(dm['with']['image_power_w']/metrics['with']['image_power_w'], 1, places=12)
        self.assertFalse(np.allclose(defocused['image_with_w_um2'], base['image_with_w_um2'], rtol=.001, atol=1e-20))

    def test_axial_dipole_cap_integral(self):
        self.data['source_polarization_weights_with'] = np.array(((0., 0., 1.),), complex)
        _, metrics = calculate_channel(self.data, self.p, replace(self.settings, pupil_samples=64), 0)
        c = np.sqrt(1-.4**2)
        exact = self.p.k0**4/(16*np.pi**2)*2*np.pi*(2/3-c+c**3/3)
        self.assertLess(abs(metrics['with']['collected_dsc_um2']/exact-1), .01)

    def test_source_translation_inverts_in_image(self):
        base, _ = calculate_channel(self.data, self.p, self.settings, 0)
        axis = base['image_axis_um']
        shift = 5*(axis[1]-axis[0])/self.settings.magnification
        self.data['source_points_with'] = np.array(((shift, 0, 0),))
        moved, _ = calculate_channel(self.data, self.p, self.settings, 0)
        expected = np.roll(base['image_with_w_um2'], -5, axis=1)
        np.testing.assert_allclose(moved['image_with_w_um2'], expected, rtol=1e-10, atol=1e-20)

    def test_small_nonzero_na_remains_active(self):
        arrays, metrics = calculate_channel(self.data, self.p,
            replace(self.settings, numerical_apertures=(.0001, 0, 0)), 0)
        self.assertTrue(arrays['pupil_mask'].any())
        self.assertGreater(metrics['with']['image_power_w'], 0)
        self.assertLess(metrics['with']['flux_relative_error'], 1e-11)


if __name__ == '__main__':
    unittest.main()

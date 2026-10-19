import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from scatmech_ui.imaging import (
    ImagingParameters, calculate_imaging, collimated_pupil,
    export_imaging, mode_weights, propagate_pupil, pupil_geometry,
)
from scatmech_ui.parameters import SimulationParameters, ObjectiveCollectionParameters
from scatmech_ui.surface_model import SCATMECH_IMPORT_ERROR, build_scattering_model


class VectorPropagationTests(unittest.TestCase):
    def setUp(self):
        self.parameters = ImagingParameters(pupil_samples=32)
        self.local, self.global_k, self.rotation, self.mask, self.du = pupil_geometry(.4, 45, 32)
        self.pupil = np.zeros((32, 32, 2, 2), dtype=complex)
        self.pupil[..., 0, 0] = self.mask

    def test_vector_rotation_preserves_solid_angle_weighted_power(self):
        jones = np.zeros((32, 32, 2, 2), dtype=complex)
        jones[..., 0, 0] = self.mask
        jones[..., 1, 1] = self.mask * (0.4+0.7j)
        pupil = collimated_pupil(jones, self.global_k, self.local, self.rotation)
        np.testing.assert_allclose(np.sum(abs(pupil)**2, axis=(-2, -1)),
                                   np.sum(abs(jones)**2, axis=(-2, -1))/self.local[..., 2], atol=1e-12)

    def test_parseval_axial_flux_at_focus_and_defocus(self):
        expected = np.sum(abs(self.pupil)**2) * self.du**2
        for defocus in (0, 300):
            params = replace(self.parameters, image_defocus_um=defocus)
            axis, e, h = propagate_pupil(self.pupil, self.local, self.du, .266, params)
            flux = np.real(e[..., 0, :]*h[..., 1, :].conj()-e[..., 1, :]*h[..., 0, :].conj())
            self.assertAlmostEqual(float(flux.sum()*(axis[1]-axis[0])**2), expected, places=12)

    def test_magnification_changes_size_and_density_not_power(self):
        axis1, e1, _ = propagate_pupil(self.pupil, self.local, self.du, .266, self.parameters)
        axis2, e2, _ = propagate_pupil(self.pupil, self.local, self.du, .266,
                                     replace(self.parameters, magnification=40))
        np.testing.assert_allclose(axis2, 2*axis1)
        # Vector refocusing slightly changes components; total integrated flux is tested separately.
        self.assertAlmostEqual(float(np.max(abs(e1)**2)/np.max(abs(e2)**2)), 4, delta=.003)

    def test_particle_translation_produces_inverted_magnified_image(self):
        axis, _, _ = propagate_pupil(self.pupil, self.local, self.du, .266, self.parameters)
        shift_pixels = 5
        object_shift = shift_pixels*(axis[1]-axis[0])/self.parameters.magnification
        phase = np.exp(-2j*np.pi/.266*self.local[..., 0]*object_shift)
        _, e, _ = propagate_pupil(self.pupil*phase[..., None, None], self.local,
                                 self.du, .266, self.parameters)
        peak = np.unravel_index(np.sum(abs(e)**2, axis=(-2, -1)).argmax(), e.shape[:2])
        self.assertEqual(peak, (len(axis)//2, len(axis)//2-shift_pixels))

    def test_invalid_and_nonfinite_inputs(self):
        for options in (dict(magnification=0), dict(transmission=2),
                        dict(irradiance_mw_um2=float('nan')), dict(pupil_samples=33)):
            with self.assertRaises(ValueError):
                ImagingParameters(**options).validate()


@unittest.skipIf(SCATMECH_IMPORT_ERROR is not None, "SCATMECH unavailable")
class NativeImagingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.simulation = SimulationParameters(wavelength_nm=266)
        cls.objective = ObjectiveCollectionParameters(channel_numerical_apertures=(.4, .4, 0))
        cls.imaging = ImagingParameters(pupil_samples=32)
        cls.result = calculate_imaging(cls.simulation, cls.objective, cls.imaging)

    def test_native_pupil_jones_matches_mueller(self):
        model = build_scattering_model(self.simulation).model
        from pySCATMECH.mueller import Polarization
        channel = self.result['channels'][0]
        local, k, _, mask, du = pupil_geometry(.4, 45, 32)
        row, col = np.argwhere(mask)[len(np.argwhere(mask))//3]
        theta = np.arccos(k[row, col, 2])
        phi = np.arctan2(-k[row, col, 0], k[row, col, 1])
        for index, pol in enumerate(('s', 'p')):
            dsc = model.DSC(np.deg2rad(self.simulation.incidence_deg), theta, phi,
                            inc=Polarization(pol))
            self.assertAlmostEqual(float(np.sum(abs(channel['jones_dsc_um'][row, col, :, index])**2)), dsc, places=13)

    def test_disabled_channel_and_power_conservation(self):
        self.assertEqual(self.result['tube_focal_length_mm'], 400)
        self.assertEqual(self.result['channels'][2]['power_w'], 0)
        for c in self.result['channels'][:2]:
            self.assertGreater(c['power_w'], 0)
            self.assertLess(c['flux_relative_error'], 1e-10)

    def test_unpolarized_is_incoherent_average_and_irradiance_scales(self):
        objective = replace(self.objective, channel_numerical_apertures=(.4, 0, 0))
        s = calculate_imaging(replace(self.simulation, polarization='s'), objective, self.imaging)
        u = calculate_imaging(replace(self.simulation, polarization='unpolarized'), objective,
                              replace(self.imaging, irradiance_mw_um2=2*self.imaging.irradiance_mw_um2))
        expected = s['channels'][0]['irradiance_w_um2'] + self.result['channels'][0]['irradiance_w_um2']
        np.testing.assert_allclose(u['channels'][0]['irradiance_w_um2'], expected, atol=1e-25)

    def test_size_change_and_analyzer_change_image(self):
        objective = replace(self.objective, channel_numerical_apertures=(.4, 0, 0))
        larger = calculate_imaging(replace(self.simulation, particle_diameter_nm=150), objective, self.imaging)
        self.assertFalse(np.allclose(larger['channels'][0]['irradiance_w_um2'],
                                    self.result['channels'][0]['irradiance_w_um2'], rtol=.01, atol=1e-25))
        x = calculate_imaging(self.simulation, objective, replace(self.imaging, analyzer='x'))
        y = calculate_imaging(self.simulation, objective, replace(self.imaging, analyzer='y'))
        self.assertAlmostEqual(x['channels'][0]['power_w']+y['channels'][0]['power_w'],
                               self.result['channels'][0]['power_w'], delta=1e-20)

    def test_roughness_does_not_silently_produce_coherent_image(self):
        with self.assertRaisesRegex(ValueError, 'roughness'):
            calculate_imaging(replace(self.simulation, roughness_enabled=True), self.objective, self.imaging)

    def test_export_contains_complex_fields_and_parameters(self):
        from matplotlib.figure import Figure
        with tempfile.TemporaryDirectory() as directory:
            export_imaging(self.result, directory, Figure())
            with np.load(str(Path(directory)/'coherent_imaging.npz')) as arrays:
                self.assertTrue(np.iscomplexobj(arrays['ch1_image_e_modes']))
            self.assertTrue((Path(directory)/'imaging_parameters.json').is_file())


if __name__ == '__main__':
    unittest.main()

"""Phase translation is geometric and must never renormalize power."""
import unittest

import numpy as np

from scatmech_ui.field_reference import sphere_to_wafer_phase
from scatmech_ui.surface_model import SCATMECH_IMPORT_ERROR


class PhaseReferenceTests(unittest.TestCase):
    def test_phase_is_incident_and_outgoing_translation_with_unit_modulus(self):
        directions = np.array([[0., 0., 1.], [.6, 0., .8], [0., .8, .6]])
        phase = sphere_to_wafer_phase(directions, .266, .05, 75.)
        inc = np.exp(-2j*np.pi/.266*.05*np.cos(np.deg2rad(75.)))
        out = np.exp(-2j*np.pi/.266*.05*directions[:, 2])
        np.testing.assert_allclose(phase, inc*out, atol=1e-15)
        np.testing.assert_allclose(abs(phase), 1., atol=1e-15)
        np.testing.assert_array_equal(sphere_to_wafer_phase(directions, .266, 0., 0.), np.ones(3))

    def test_invalid_physical_inputs_are_rejected(self):
        for args in ((0., .05, 0.), (.266, -.01, 0.), (.266, .05, 90.)):
            with self.assertRaises(ValueError):
                sphere_to_wafer_phase([[0., 0., 1.]], *args)

    @unittest.skipIf(SCATMECH_IMPORT_ERROR is not None, 'SCATMECH required')
    def test_independent_free_space_dipole_at_known_nonzero_height(self):
        from pySCATMECH.local import Local_BRDF_Model
        from fullwave_prototype.parameters import Parameters
        from fullwave_prototype.far_field import radiation_amplitude
        p = Parameters(wavelength_nm=266., diameter_nm=1., gap_nm=200.,
                       incidence_deg=35., substrate_n=1., substrate_k=0.)
        model = Local_BRDF_Model('Bobbert_Vlieger_BRDF_Model', {
            'lambda': .266, 'radius': .0005, 'sphere': complex(1.59),
            'substrate': complex(1.), 'delta': .2, 'lmax': 8, 'order': -1,
            'Norm_Inc_Approx': 0, 'improve': 3})
        directions = np.array([[0., 0., 1.], [.3, .4, np.sqrt(.75)],
                               [-.6, 0., .8], [0., -.8, .6]])
        raw = np.array([model.JonesDSC(np.deg2rad(35.), np.arccos(k[2]),
                        np.arctan2(-k[0], k[1]), coords='psps')[:, 1] for k in directions])
        height = .2005
        translated = raw*sphere_to_wafer_phase(directions, .266, height, 35.)[:, None]
        ei = np.array([0., np.cos(np.deg2rad(35.)), np.sin(np.deg2rad(35.))])
        volume = 4*np.pi*.0005**3/3
        moment = 3*volume*(1.59**2-1)/(1.59**2+2)*ei*np.exp(-1j*p.k0*height*np.cos(np.deg2rad(35.)))
        expected = radiation_amplitude(np.array([[0., 0., height]]), moment[None, :], directions, p)
        phase_only = np.exp(-1j*np.angle(np.vdot(expected, translated)))
        relative_l2 = np.linalg.norm(translated*phase_only-expected)/np.linalg.norm(expected)
        self.assertLess(relative_l2, 1e-4)
        self.assertGreater(np.linalg.norm(raw*phase_only-expected)/np.linalg.norm(expected), .1)

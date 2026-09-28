import unittest
from unittest.mock import patch
import numpy as np

from fullwave_prototype.physical_audit import sphere_rule
from fullwave_prototype.audit_optics import airy_check


class IndependentValidationTests(unittest.TestCase):
    def test_invalid_pupil_sampling_rejected_before_fem(self):
        from fullwave_prototype.physical_audit import flat_pupil
        from fullwave_prototype.parameters import Parameters
        with patch('fullwave_prototype.physical_audit.solve') as solve:
            with self.assertRaises(ValueError):
                flat_pupil(Parameters(), None, None, pupil_samples=31)
            solve.assert_not_called()

    def test_spherical_measure_and_moments(self):
        d, w = sphere_rule(12)
        self.assertAlmostEqual(w.sum(), 4*np.pi, places=12)
        np.testing.assert_allclose((d*w[:, None]).sum(axis=0), 0, atol=1e-13)
        np.testing.assert_allclose((d*d*w[:, None]).sum(axis=0), 4*np.pi/3, atol=1e-13)
        d, w = sphere_rule(12, hemisphere=True)
        self.assertAlmostEqual(w.sum(), 2*np.pi, places=12)
        self.assertGreater(d[:, 2].min(), 0)

    def test_airy_absolute_profile_in_paraxial_limit(self):
        report, _ = airy_check()
        self.assertLess(report['absolute_slice_relative_l2'], .01)
        self.assertLess(report['peak_relative_error'], .01)


if __name__ == '__main__':
    unittest.main()

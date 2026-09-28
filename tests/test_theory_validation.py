"""Independent formula checks and regressions found by the theory audit."""
from dataclasses import replace
import importlib.util
import unittest

import numpy as np

from optical_propagation import pupil_geometry, propagate_pupil
from fullwave_prototype.imaging import ImagingSettings, axial_flux
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.theory_optics import pec_image_check, vector_focus_check
from fullwave_prototype.theory_reference import bessel_mie_cross_sections, mie_fields


class TheoryTests(unittest.TestCase):
    def test_huygens_surface_against_analytic_dipole(self):
        from fullwave_prototype.physical_audit import sphere_rule
        from fullwave_prototype.huygens import surface_amplitude
        n, weights = sphere_rule(20)
        radius, k0 = .08, 2*np.pi/.355
        moment = np.array((.3, .2j, 1))
        center = np.array((.01, -.02, .03))
        ndotp = (n@moment)[:, None]
        phase = np.exp(1j*k0*radius)/(4*np.pi)
        e = phase*(k0**2/radius*(moment-n*ndotp)+
                   (1/radius**3-1j*k0/radius**2)*(3*n*ndotp-moment))
        h = phase*(k0**2/radius+1j*k0/radius**2)*np.cross(n, moment)
        directions = n[::8]
        actual = surface_amplitude(center+radius*n, n, weights*radius**2, e, h, directions, k0)
        expected = k0**2/(4*np.pi)*(moment-directions*(directions@moment)[:, None])
        expected *= np.exp(-1j*k0*(directions@center))[:, None]
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)

    def test_high_na_empty_corners_and_power_conservation(self):
        settings = ImagingSettings(magnification=1, pupil_samples=32, padding=4)
        local, _, _, mask, du = pupil_geometry(.95, 0, 32)
        pupil = np.zeros(mask.shape+(2, 2), complex)
        pupil[..., 0, 0] = mask
        with np.errstate(all='raise'):
            axis, e, h = propagate_pupil(pupil, local, du, .355, settings)
        self.assertTrue(np.isfinite(e).all())
        power = axial_flux(e[..., 0], h[..., 0]).sum()*(axis[1]-axis[0])**2
        self.assertAlmostEqual(power, mask.sum()*du**2, places=11)
        pupil[0, 0, 0, 0] = 1
        with self.assertRaisesRegex(ValueError, 'non-propagating'):
            propagate_pupil(pupil, local, du, .355, settings)

    def test_nonparaxial_hankel_reference(self):
        reports, _ = vector_focus_check()
        for report in reports:
            self.assertLess(report['electric_relative_l2'], .004)
            self.assertLess(report['magnetic_relative_l2'], .004)
            self.assertLess(report['flux_relative_l2'], .004)
            self.assertLess(report['oracle_quadrature_relative'], 1e-10)

    def test_pec_image_limit(self):
        reports = pec_image_check()
        for coarse, fine in (reports[:2], reports[2:]):
            self.assertLess(fine['relative_complex_l2'], 4e-5)
            self.assertLess(fine['relative_complex_l2'], coarse['relative_complex_l2']/90)

    def test_mie_rayleigh_limit_and_absorption(self):
        p = Parameters(diameter_nm=1, substrate_n=1, substrate_k=0)
        alpha = 4*np.pi*(p.diameter_nm/2000)**3*(p.eps_particle-1)/(p.eps_particle+2)
        rayleigh = p.k0**4/(6*np.pi)*abs(alpha)**2
        exact = bessel_mie_cross_sections(p)
        self.assertLess(abs(exact['csca_um2']/rayleigh-1), 2e-5)
        self.assertAlmostEqual(exact['cext_um2']/exact['csca_um2'], 1, places=10)
        lossy = bessel_mie_cross_sections(replace(p, particle_k=.1))
        self.assertGreater(lossy['cext_um2'], lossy['csca_um2'])

    def test_integration_bonus_validation(self):
        for bad in (-1, 1.5, 13):
            with self.assertRaises(ValueError):
                Parameters(integration_bonus=bad).validate()


@unittest.skipUnless(importlib.util.find_spec('miepython'), 'Independent Mie oracle required')
class MieOracleTests(unittest.TestCase):
    def test_frame_and_phase_against_no_contrast_plane_wave(self):
        from fullwave_prototype.planar import background_fields
        points = np.array(((.01, .015, .06), (.08, -.03, .12), (-.06, .01, .03)))
        for pol in ('s', 'p'):
            p = Parameters(particle_n=1, substrate_n=1, substrate_k=0, gap_nm=10, polarization=pol)
            e, h = mie_fields(points, p)
            eb, hb = background_fields(points, p)
            np.testing.assert_allclose(e, eb, atol=2e-6, rtol=2e-6)
            np.testing.assert_allclose(h, hb, atol=2e-6, rtol=2e-6)

    def test_mie_coefficients_independent_cross_check(self):
        import miepython
        for size, loss in ((30, 0), (100, .1), (200, 0)):
            p = Parameters(diameter_nm=size, particle_k=loss)
            exact = bessel_mie_cross_sections(p)
            ext, sca, _, _ = miepython.efficiencies(complex(p.particle_n, -loss), size, p.wavelength_nm)
            area = np.pi*(size/2000)**2
            self.assertAlmostEqual(ext*area/exact['cext_um2'], 1, places=10)
            self.assertAlmostEqual(sca*area/exact['csca_um2'], 1, places=10)


if __name__ == '__main__':
    unittest.main()

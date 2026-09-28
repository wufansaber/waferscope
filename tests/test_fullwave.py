"""Analytic checks also runnable in NIST2; optional native FEM smoke test."""
import importlib.util
import unittest
import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.planar import fresnel, background_fields, background_vectors
from fullwave_prototype.far_field import sphere_quadrature, radiation_amplitude, angular_grid


class FullwaveAnalyticTests(unittest.TestCase):
    def test_parameter_validation(self):
        for kwargs in (dict(gap_nm=-1), dict(wavelength_nm=0), dict(order=4),
                       dict(polarization='unpolarized'), dict(substrate_k=-1),
                       dict(particle_mesh_um=float('nan')), dict(domain_radius_um=.05)):
            with self.assertRaises(ValueError):
                Parameters(**kwargs).validate()

    def test_fresnel_normal_incidence_and_no_contrast(self):
        rs, rp, q = fresnel(1., 2.25)
        self.assertAlmostEqual(rs, -.2)
        self.assertAlmostEqual(rp, .2)
        rs, rp, q = fresnel(.3, 1.)
        self.assertAlmostEqual(rs, 0)
        self.assertAlmostEqual(rp, 0)
        rs, rp, _ = fresnel(0., 1.)
        self.assertEqual(rs, 0)
        self.assertEqual(rp, 0)
        self.assertGreaterEqual(fresnel(.3, (5.57+.387j)**2)[2].imag, 0)

    def test_incidence_and_boundary_tangential_continuity(self):
        for pol in ('s', 'p'):
            p = Parameters(polarization=pol)
            down, up, kt, ei, er, et = background_vectors(p)
            self.assertGreater(down[1], 0)
            self.assertLess(down[2], 0)
            self.assertAlmostEqual(np.dot(down, ei), 0)
            np.testing.assert_allclose((ei+er)[:2], et[:2], atol=1e-14)
            np.testing.assert_allclose((np.cross(down, ei)+np.cross(up, er))[:2],
                                       np.cross(kt, et)[:2], atol=1e-14)

    def test_substrate_wave_decays_inward(self):
        e, _ = background_fields(np.array(((0, 0, -.01), (0, 0, -.2))), Parameters())
        self.assertLess(np.linalg.norm(e[1]), np.linalg.norm(e[0]))

    def test_sphere_quadrature_volume_center_and_second_moment(self):
        p = Parameters()
        points, w = sphere_quadrature(p, 5)
        radius = p.diameter_nm/2000
        volume = 4*np.pi*radius**3/3
        self.assertAlmostEqual(w.sum(), volume, places=15)
        np.testing.assert_allclose((points*w[:, None]).sum(axis=0)/volume,
                                   (0, 0, radius+p.gap_nm/1000), atol=1e-14)
        self.assertAlmostEqual(np.sum(points[:, 0]**2*w)/volume, radius**2/5, places=14)

    def test_free_space_point_dipole_radiation(self):
        p = Parameters(substrate_n=1, substrate_k=0)
        theta, phi, directions = angular_grid(10)
        k = directions.reshape(-1, 3)
        moment = np.array((1., 0, 0))
        amplitude = radiation_amplitude(np.zeros((1, 3)), moment[None, :], k, p)
        expected = p.k0**4/(16*np.pi**2)*(1-k[:, 0]**2)
        np.testing.assert_allclose(np.sum(abs(amplitude)**2, axis=1), expected, atol=1e-10)

    def test_translation_phase_in_free_space(self):
        p = Parameters(substrate_n=1, substrate_k=0)
        k = np.array(((.6, 0, .8), (0, .6, .8)))
        moment = np.array(((0, 1+1j, 0),))
        shift = np.array((.02, .03, .06))
        a = radiation_amplitude(np.zeros((1, 3)), moment, k, p)
        b = radiation_amplitude(shift[None, :], moment, k, p)
        np.testing.assert_allclose(b, a*np.exp(-1j*p.k0*(k@shift))[:, None], atol=1e-12)


@unittest.skipIf(importlib.util.find_spec('ngsolve') is None, 'Native NGSolve environment required')
class NativeFullwaveTests(unittest.TestCase):
    def test_zero_particle_contrast_has_zero_scattered_field(self):
        from fullwave_prototype.solver import solve, evaluate_fields
        p = Parameters(particle_n=1, particle_k=0, order=1)
        mesh, field, metrics = solve(p)
        points = np.array(((0, 0, .052), (.1, 0, .1), (0, .1, -.02)))
        e, h, esc = evaluate_fields(mesh, field, points, p)
        ebg, hbg = background_fields(points, p)
        np.testing.assert_allclose(esc, 0, atol=1e-13)
        np.testing.assert_allclose(e, ebg, atol=1e-13)
        np.testing.assert_allclose(h, hbg, atol=1e-13)
        self.assertLess(metrics['algebraic_residual'], 1e-10)


if __name__ == '__main__':
    unittest.main()

"""Roughness geometry, reciprocal radiation and coherent paired-solve checks."""
import importlib.util
import unittest
import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness, height_mesh
from fullwave_prototype.far_field import radiation_amplitude, sphere_quadrature
from fullwave_prototype.rough_run import decompose


class RoughnessAnalyticTests(unittest.TestCase):
    def test_actual_triangle_mean_rms_and_repeatability(self):
        r = Roughness()
        vertices, triangles = height_mesh(r)
        heights = triangles[..., 2]
        second = np.mean((np.sum(heights**2, axis=1)+np.sum(heights, axis=1)**2)/12)
        self.assertAlmostEqual(heights.mean(), 0, places=15)
        self.assertAlmostEqual(np.sqrt(second)*1000, r.rms_nm, places=12)
        np.testing.assert_array_equal(vertices, height_mesh(r)[0])
        np.testing.assert_array_equal(vertices[[0, -1], :, 2], 0)
        np.testing.assert_array_equal(vertices[:, [0, -1], 2], 0)
        np.testing.assert_array_equal(height_mesh(Roughness(rms_nm=0))[0][..., 2], 0)
        self.assertFalse(np.allclose(vertices, height_mesh(Roughness(seed=8))[0]))

    def test_geometry_guards(self):
        for r, p in ((Roughness(), Parameters()),
                     (Roughness(width_um=1), Parameters(gap_nm=10)),
                     (Roughness(rms_nm=-1), Parameters(gap_nm=10)),
                     (Roughness(grid_points=3), Parameters(gap_nm=10)),
                     (Roughness(seed=-1), Parameters(gap_nm=10))):
            with self.assertRaises(ValueError):
                r.validate(p)
        Roughness().validate(Parameters(gap_nm=10))

    def test_lower_source_becomes_free_space_without_substrate_contrast(self):
        p = Parameters(substrate_n=1, substrate_k=0)
        directions = np.array(((.6, 0, .8), (0, .6, .8), (0, 0, 1), (1, 0, 0)))
        source = np.array(((1+1j, 2, 3j),))
        origin = radiation_amplitude(np.zeros((1, 3)), source, directions, p)
        shift = np.array((.01, -.03, -.05))
        moved = radiation_amplitude(shift[None, :], source, directions, p)
        np.testing.assert_allclose(moved, origin*np.exp(-1j*p.k0*(directions@shift))[:, None], atol=1e-12)

    def test_reciprocal_kernel_interface_conditions(self):
        p = Parameters()
        directions = np.array(((.6, 0, .8), (0, .6, .8), (0, 0, 1)))
        for component in range(3):
            source = np.eye(3)[component:component+1]
            upper = radiation_amplitude(np.array(((0, 0, 1e-12),)), source, directions, p)
            lower = radiation_amplitude(np.array(((0, 0, -1e-12),)), source, directions, p)
            # Tangential reciprocal E is continuous; normal eps*E is continuous.
            factor = p.eps_substrate if component == 2 else 1
            np.testing.assert_allclose(upper, factor*lower, rtol=1e-8, atol=1e-8)

    def test_signed_interference_not_incoherent_power_addition(self):
        a0 = np.array(((1+1j, 2-1j),))
        a1 = .5*a0
        d = decompose(a1, a0)
        np.testing.assert_allclose(d['total'], d['background']+d['change_field']+d['interference'])
        self.assertLess(d['interference'][0], 0)
        self.assertGreater(d['change_field'][0], 0)
        self.assertLess((d['total']-d['background'])[0], 0)

    def test_transmitted_p_kernel_at_critical_angle(self):
        from fullwave_prototype.planar import background_vectors
        p = Parameters(substrate_n=.6, substrate_k=0, incidence_deg=np.rad2deg(np.arcsin(.6)))
        *_, et = background_vectors(p)
        self.assertTrue(np.isfinite(et).all())
        self.assertAlmostEqual(et[2].real, 2*.6/p.eps_substrate.real, places=6)
        k = np.array(((.6, 0, .8),))
        a = radiation_amplitude(np.array(((0, 0, -1e-12),)), np.array(((0, 0, 1),)), k, p)
        expected = p.k0**2/(4*np.pi)*2*.6/p.eps_substrate
        self.assertAlmostEqual(a[0, 1].real/expected.real, 1, places=6)


@unittest.skipIf(importlib.util.find_spec('ngsolve') is None, 'Native NGSolve environment required')
class NativeRoughnessTests(unittest.TestCase):
    def test_closed_solid_and_contrast_volume(self):
        import ngsolve as ng
        from fullwave_prototype.rough_surface import surface_solid, build_rough_mesh
        r = Roughness()
        self.assertAlmostEqual(surface_solid(r, .03).mass, r.width_um**2*.03, places=12)
        mesh = build_rough_mesh(Parameters(gap_nm=10, order=1), r)
        positive = ng.Integrate(1, mesh, definedon=mesh.Materials('bump'))
        negative = ng.Integrate(1, mesh, definedon=mesh.Materials('pit'))
        self.assertGreater(positive, 0)
        self.assertAlmostEqual(positive, negative, places=12)

    def test_identical_material_pair_has_zero_field_change(self):
        from fullwave_prototype.rough_surface import build_rough_mesh
        from fullwave_prototype.solver import solve
        p = Parameters(gap_nm=10, particle_n=1, order=1)
        mesh = build_rough_mesh(p, Roughness())
        _, a, _ = solve(p, mesh=mesh, particle_enabled=True)
        _, b, _ = solve(p, mesh=mesh, particle_enabled=False)
        np.testing.assert_allclose(a.vec.FV().NumPy(), b.vec.FV().NumPy(), atol=1e-12)

    def test_flat_limit_and_element_source_quadrature(self):
        from fullwave_prototype.rough_surface import build_rough_mesh
        from fullwave_prototype.solver import solve, evaluate_fields
        from fullwave_prototype.volume_sources import polarization_sources
        p = Parameters(gap_nm=10, order=1)
        mesh = build_rough_mesh(p, Roughness(rms_nm=0))
        self.assertNotIn('bump', mesh.GetMaterials())
        _, field, _ = solve(p, mesh=mesh)
        points, sources, volumes = polarization_sources(mesh, field, p)
        self.assertAlmostEqual(volumes['particle']/(4*np.pi*(p.diameter_nm/2000)**3/3), 1, places=3)
        directions = np.array(((0, 0, 1), (0, .6, .8), (.6, 0, .8)))
        actual = radiation_amplitude(points, sources, directions, p)
        qp, qw = sphere_quadrature(p, 14)
        electric, _, _ = evaluate_fields(mesh, field, qp, p)
        expected = radiation_amplitude(qp, (p.eps_particle-1)*electric*qw[:, None], directions, p)
        self.assertLess(np.linalg.norm(actual-expected)/np.linalg.norm(expected), .005)


if __name__ == '__main__':
    unittest.main()

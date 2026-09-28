import importlib.util
from dataclasses import replace
import unittest

import numpy as np

from fullwave_prototype.parameters import Parameters
from fullwave_prototype.absorbing_boundary import stretch_factors
from fullwave_prototype.numerical_settings import numerical_parameters


class AbsorbingBoundaryTests(unittest.TestCase):
    def test_profiles_and_explicit_overrides(self):
        self.assertEqual(numerical_parameters('refined').order, 3)
        self.assertTrue(numerical_parameters('refined').static_condensation)
        balanced = numerical_parameters('balanced')
        self.assertEqual((balanced.order, balanced.pml_strength, balanced.pml_mesh_scale), (3, 1., 1.25))
        self.assertTrue(balanced.static_condensation)
        self.assertEqual(balanced.direct_solver, 'sparsecholesky')
        self.assertEqual(numerical_parameters('balanced', pml_strength=1.5).pml_strength, 1.5)
        self.assertEqual(numerical_parameters(), Parameters())
        self.assertEqual(numerical_parameters('preview'), Parameters())
        p = numerical_parameters(mesh_scale=.8, gap_nm=10, order=2, pml_power=2)
        self.assertEqual(p.order, 2)
        self.assertEqual(p.pml_power, 2)
        self.assertAlmostEqual(p.air_mesh_um, .064)
        self.assertEqual(p.gap_nm, 10)
        with self.assertRaises(ValueError):
            numerical_parameters(mesh_scale=-1)

    def test_interface_and_outer_stretch(self):
        for power in (1, 2, 3):
            p = Parameters(pml_power=power)
            r = p.domain_radius_um
            tangential, radial = stretch_factors(r, p)
            self.assertEqual(tangential, 1)
            if power > 1:
                self.assertEqual(radial, 1)
            outer = r+p.pml_thickness_um
            t, _ = stretch_factors(outer, p)
            self.assertAlmostEqual((outer*t).imag, p.pml_strength*p.pml_thickness_um)

    def test_full_jacobian_by_finite_differences(self):
        x = np.array((.2, .1, .22))
        for power in (1, 2, 3):
            p = Parameters(pml_power=power)
            r = np.linalg.norm(x)
            t, radial = stretch_factors(r, p)
            expected = t*np.eye(3)+(radial-t)*np.outer(x, x)/r**2
            def transform(point):
                return point*stretch_factors(np.linalg.norm(point), p)[0]
            delta = np.eye(3)*1e-7
            actual = np.column_stack([(transform(x+d)-transform(x-d))/2e-7 for d in delta])
            np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-8)

    def test_new_controls_validate_and_old_reports_load(self):
        Parameters(**{'order': 2, 'gap_nm': 10}).validate()
        for update in (dict(pml_power=0), dict(pml_power=1.5), dict(pml_mesh_scale=0),
                       dict(pml_mesh_scale=float('nan')), dict(static_condensation='yes'),
                       dict(direct_solver='invalid')):
            with self.assertRaises((ValueError, TypeError)):
                Parameters(**update).validate()


@unittest.skipIf(importlib.util.find_spec('ngsolve') is None, 'Native NGSolve required')
class NativeCondensationTests(unittest.TestCase):
    def test_symmetric_solver_matches_umfpack_on_lossy_wafer(self):
        from fullwave_prototype.solver import solve
        p = Parameters(order=3, gap_nm=10, static_condensation=True,
                       domain_radius_um=.14, pml_thickness_um=.10,
                       air_mesh_um=.1, substrate_mesh_um=.08, particle_mesh_um=.035)
        mesh, reference, _ = solve(p)
        _, symmetric, info = solve(replace(p, direct_solver='sparsecholesky'), mesh=mesh)
        np.testing.assert_allclose(symmetric.vec.FV().NumPy(), reference.vec.FV().NumPy(), rtol=1e-8, atol=1e-10)
        self.assertLess(info['algebraic_residual'], 1e-9)

    def test_condensed_and_full_fields_agree(self):
        from fullwave_prototype.solver import solve, evaluate_fields
        # Coarse, low-index case tests algebra, not physical convergence.
        p = Parameters(order=3, substrate_n=1, substrate_k=0, gap_nm=10,
                       domain_radius_um=.14, pml_thickness_um=.10,
                       air_mesh_um=.1, substrate_mesh_um=.1, particle_mesh_um=.035)
        mesh, full, _ = solve(p)
        _, condensed, info = solve(replace(p, static_condensation=True), mesh=mesh)
        np.testing.assert_allclose(condensed.vec.FV().NumPy(), full.vec.FV().NumPy(), rtol=1e-8, atol=1e-10)
        self.assertLess(info['algebraic_residual'], 1e-9)
        self.assertLess(info['factored_dofs'], info['dofs'])
        points = np.array(((.06, .02, .04), (0, -.03, .07)))
        for a, b in zip(evaluate_fields(mesh, full, points, p), evaluate_fields(mesh, condensed, points, p)):
            np.testing.assert_allclose(a, b, rtol=1e-8, atol=1e-10)


if __name__ == '__main__':
    unittest.main()

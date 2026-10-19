"""Joint-process refinement preserves common nodes and parent statistics."""
from dataclasses import replace
import unittest

import numpy as np

from fullwave_prototype.coupled_surface import (height_grid, joint_factor, joint_nodes,
                                               realization_signature)
from fullwave_prototype.haze_convergence import compare_cases, default_plan, case_settings
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness


class CoupledSurfaceTests(unittest.TestCase):
    def test_nested_common_nodes_are_identical_not_interpolated(self):
        for model in ('gaussian', 'abc'):
            r = Roughness(model=model, grid_points=9, sampling_grid_levels=(9, 17),
                          correlation_length_um=.2, psd_exponent=4.5)
            fine = replace(r, grid_points=17)
            np.testing.assert_array_equal(height_grid(r), height_grid(fine)[::2, ::2])
            self.assertEqual(realization_signature(r), realization_signature(fine))
            self.assertNotEqual(realization_signature(r), realization_signature(replace(r, seed=8)))

    def test_joint_covariance_matches_analytic_parent(self):
        from scipy.spatial.distance import cdist
        from fullwave_prototype.rough_statistics import abc_correlation
        nodes, _ = joint_nodes((5, 9))
        distance = cdist(nodes*.2, nodes*.2)
        for model in ('gaussian', 'abc'):
            factor = joint_factor(model, (5, 9), .2, .05, 4.5)
            expected = (np.exp(-(distance/.05)**2) if model == 'gaussian' else
                        abc_correlation(distance, .05, 4.5))
            np.testing.assert_allclose(factor@factor.T, expected, atol=2e-12)

    def test_non_nested_grids_share_coordinates(self):
        r = Roughness(model='abc', grid_points=9, sampling_grid_levels=(9, 13, 17))
        a, b = height_grid(r), height_grid(replace(r, grid_points=13))
        np.testing.assert_array_equal(a[::2, ::2], b[::3, ::3])

    def test_grid_cases_can_use_paired_delta_without_same_mesh(self):
        plan = default_plan(Parameters(gap_nm=10), Roughness(model='abc', grid_points=9))
        r = case_settings(plan, 'height_grid_points', 13)[1]
        self.assertEqual(r.sampling_grid_levels, [9, 13, 17])
        def case(value, height_hash):
            return dict(surface_sha256=[height_hash]*3, realization_sha256=['parent']*3,
                statistics=dict(sample_observables_um2=[[value]*4]*3, numerical_apertures=[.5]*3))
        result = compare_cases(case(1, 'coarse'), case(1.01, 'fine'), 'height_grid_points', .05)
        self.assertFalse(result['paired_same_surface'])
        self.assertTrue(result['paired_same_realization'])
        self.assertEqual(result['status'], 'within_tested_tolerance')

    def test_joint_node_budget_and_invalid_grid(self):
        with self.assertRaisesRegex(ValueError, '2048'):
            joint_nodes((37, 39, 41))
        with self.assertRaises(ValueError):
            joint_nodes((3, 9))


if __name__ == '__main__':
    unittest.main()

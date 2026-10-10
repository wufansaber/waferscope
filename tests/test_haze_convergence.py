"""Confidence intervals, exact-height pairing and scan-plan semantics."""
from dataclasses import asdict
import unittest
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np

from fullwave_prototype.haze_convergence import (case_settings, compare_cases, default_plan,
                                               run_plan, SampleLimitError, validate_plan)
from fullwave_prototype.haze_convergence import _sample
from fullwave_prototype.haze_sample_cache import completed_sample
from fullwave_prototype.haze_statistics import confidence_interval, observable_samples, statistical_report
from fullwave_prototype.haze_statistics import estimated_sample_budget
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.rough_surface import Roughness


class HazeConvergenceTests(unittest.TestCase):
    def test_student_t_not_normal_small_sample_interval(self):
        interval = confidence_interval([1, 2, 3])
        self.assertAlmostEqual(interval['mean_um2'], 2.)
        self.assertAlmostEqual(interval['standard_error_um2'], 1/np.sqrt(3))
        self.assertAlmostEqual(interval['ci95_um2'][1], 2+4.30265273/np.sqrt(3), places=7)
        self.assertIsNone(confidence_interval([1])['ci95_um2'])

    def test_apertures_solid_angle_and_disabled_channel(self):
        theta = np.deg2rad(np.arange(0, 90.25, .5))
        phi = np.deg2rad(np.arange(360))
        values = observable_samples(np.ones((1, len(theta), len(phi))), theta, phi,
                                   apertures=(.5, .5, 0), tilts=(0, 45, 0))[0]
        self.assertAlmostEqual(values[0], 2*np.pi, places=4)
        self.assertAlmostEqual(values[1], 2*np.pi*(1-np.sqrt(.75)), delta=.018)
        self.assertAlmostEqual(values[1], values[2], delta=.018)
        self.assertEqual(values[3], 0.)

    def test_statistics_do_not_report_one_surface_as_converged(self):
        theta, phi = np.deg2rad([0, 45, 90]), np.deg2rad([0, 90, 180, 270])
        one = statistical_report(np.ones((1, 3, 4)), theta, phi)
        self.assertFalse(one['statistical_target_met'])
        identical = statistical_report(np.ones((4, 3, 4)), theta, phi)
        self.assertTrue(identical['statistical_target_met'])
        self.assertEqual(estimated_sample_budget(identical)['estimated_total_surfaces'], 3)
        self.assertIsNone(estimated_sample_budget(one)['estimated_total_surfaces'])

    def test_polynomial_order_is_an_explicit_paired_scan(self):
        plan = default_plan(Parameters(gap_nm=10, order=1, pml_order=1),
                            Roughness(model='abc', grid_points=9), 3)
        validate_plan(plan)
        self.assertEqual(plan['scans']['fem_order'], [1, 2])
        p, r, _, _ = case_settings(plan, 'fem_order', 2)
        self.assertEqual(p.order, 2)
        self.assertEqual(p.pml_order, 1)
        self.assertEqual(r, Roughness(**plan['roughness']))

    def test_coupled_grid_cache_request_is_json_stable(self):
        r = Roughness(model='abc', grid_points=9, sampling_grid_levels=(9, 13, 17))
        with tempfile.TemporaryDirectory() as directory:
            def guarded(command, report_path, **kwargs):
                folder = Path(report_path).parent
                (folder/'sample_report.json').write_text('{}', encoding='utf-8')
                (folder/'haze_sources.npz').touch()
                Path(report_path).write_text(json.dumps(dict(status='completed', returncode=0,
                    surviving_owned_pids=[])), encoding='utf-8')
                return dict(status='completed')
            with patch('fullwave_prototype.haze_convergence.run_guarded', side_effect=guarded) as run:
                first = _sample(Parameters(), r, Path(directory), 5, 10, 180)
                second = _sample(Parameters(), r, Path(directory), 5, 10, 180)
            self.assertEqual(first, second)
            self.assertEqual(run.call_count, 1)

    def test_local_mesh_and_pml_scans_preserve_physical_surface(self):
        plan = default_plan(Parameters(gap_nm=10, order=2, pml_order=1),
                            Roughness(model='abc', grid_points=9, taper_fraction=.125), 3)
        validate_plan(plan)
        p, r, _, _ = case_settings(plan, 'rough_mesh_scale', .8)
        self.assertEqual(asdict(p), plan['parameters'])
        self.assertAlmostEqual(r.mesh_um, plan['roughness']['mesh_um']*.8)
        self.assertEqual(r.seed, plan['roughness']['seed'])
        p, r, _, _ = case_settings(plan, 'pml_mesh_scale', .85)
        self.assertEqual(p.air_mesh_um, plan['parameters']['air_mesh_um'])
        self.assertEqual(asdict(r), plan['roughness'])
        p, _, _, _ = case_settings(plan, 'pml_order', 2)
        self.assertEqual(p.pml_order, p.order)
        plan['scans'] = dict(pml_order=[1, 1.5])
        with self.assertRaisesRegex(ValueError, 'Invalid'):
            validate_plan(plan)

    def test_second_order_boundary_plan_does_not_retain_reduced_pml(self):
        from tools.prepare_second_order_study import prepare
        plan = default_plan(Parameters(gap_nm=10, order=1, pml_order=1),
                            Roughness(model='abc', grid_points=9, taper_fraction=.125), 3)
        boundary = prepare(plan, 'boundary')
        self.assertEqual(boundary['parameters']['order'], 2)
        self.assertEqual(boundary['parameters']['pml_order'], 2)
        self.assertEqual(plan['parameters']['order'], 1)

    def test_cache_rejects_incomplete_failed_or_mismatched_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            request = dict(seed=7)
            (folder/'input.json').write_text(json.dumps(request), encoding='utf-8')
            (folder/'sample_report.json').write_text('{}', encoding='utf-8')
            (folder/'haze_sources.npz').touch()
            self.assertFalse(completed_sample(folder, request))
            resource = dict(status='completed', returncode=0, surviving_owned_pids=[])
            for change in (dict(status='stopped_timeout'), dict(returncode=1),
                           dict(surviving_owned_pids=[123])):
                (folder/'resources.json').write_text(json.dumps(dict(resource, **change)), encoding='utf-8')
                self.assertFalse(completed_sample(folder, request))
            (folder/'resources.json').write_text(json.dumps(resource), encoding='utf-8')
            self.assertTrue(completed_sample(folder, request))
            self.assertFalse(completed_sample(folder, dict(seed=8)))
            (folder/'resources.json').write_text('{bad', encoding='utf-8')
            self.assertFalse(completed_sample(folder, request))

    def test_guard_start_failure_invalidates_stale_success(self):
        with tempfile.TemporaryDirectory() as directory:
            def guarded(command, report_path, **kwargs):
                folder = Path(report_path).parent
                (folder/'sample_report.json').write_text('{}', encoding='utf-8')
                (folder/'haze_sources.npz').touch()
                Path(report_path).write_text(json.dumps(dict(status='completed', returncode=0,
                    surviving_owned_pids=[])), encoding='utf-8')
                return dict(status='completed')
            with patch('fullwave_prototype.haze_convergence.run_guarded', side_effect=guarded):
                folder = _sample(Parameters(), Roughness(model='abc'), Path(directory), 5, 10, 180)
            (folder/'input.json').write_text('{}', encoding='utf-8')
            with patch('fullwave_prototype.haze_convergence.run_guarded',
                       side_effect=RuntimeError('Insufficient available memory')):
                with self.assertRaisesRegex(RuntimeError, 'Insufficient'):
                    _sample(Parameters(), Roughness(model='abc'), Path(directory), 5, 10, 180)
            request = json.loads((folder/'input.json').read_text(encoding='utf-8'))
            self.assertFalse(completed_sample(folder, request))

    def test_paired_mesh_change_cancels_random_surface_variance(self):
        def case(samples, hashes=('a', 'b', 'c')):
            return dict(surface_sha256=list(hashes), statistics=dict(
                sample_observables_um2=np.repeat(np.array(samples)[:, None], 4, axis=1).tolist(),
                numerical_apertures=[.5, .5, .5]))
        a, b = case([1, 20, 100]), case([1.01, 20.01, 100.01])
        result = compare_cases(a, b, 'fem_mesh_scale', .05)
        self.assertTrue(result['paired_same_surface'])
        self.assertEqual(result['status'], 'within_tested_tolerance')
        unpaired = compare_cases(a, b, 'height_grid_points', .05)
        self.assertFalse(unpaired['paired_same_surface'])
        self.assertEqual(unpaired['status'], 'unresolved')
        changed = compare_cases(a, case([1, 20, 100], ('d', 'e', 'f')), 'fem_mesh_scale', .05)
        self.assertFalse(changed['paired_same_surface'])

    def test_plan_preserves_surface_for_mesh_and_physical_taper_for_patch(self):
        r = Roughness(model='gaussian', grid_points=9, taper_fraction=.125)
        plan = default_plan(Parameters(gap_nm=10, order=1), r, 3)
        validate_plan(plan)
        _, mesh_rough, _, _ = case_settings(plan, 'fem_mesh_scale', .7)
        a, b = asdict(r), asdict(mesh_rough)
        del a['mesh_um'], b['mesh_um']
        self.assertEqual(a, b)
        _, wider, _, _ = case_settings(plan, 'patch_width_um', .24)
        self.assertAlmostEqual(wider.width_um*wider.taper_fraction, .025)
        self.assertLessEqual(wider.width_um/(wider.grid_points-1), .025)
        plan['scans']['height_grid_points'] = [9, 10.5]
        with self.assertRaisesRegex(ValueError, 'Invalid'):
            validate_plan(plan)

    def test_resource_limited_case_is_excluded_without_biasing_remaining_scan(self):
        plan = default_plan(Parameters(gap_nm=10, order=1),
                            Roughness(model='abc', grid_points=9, taper_fraction=.125), 3)
        plan['scans'] = dict(pml_strength=[2., 2.5])
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)/'sample'
            folder.mkdir()
            np.savez(folder/'haze_sources.npz', points=np.zeros((1, 3)), weights=np.ones((1, 3)))
            (folder/'sample_report.json').write_text(json.dumps(dict(surface_sha256='same',
                solver=dict(dofs=10), elapsed_seconds=1)), encoding='utf-8')
            (folder/'resources.json').write_text(json.dumps(dict(peak_accounted_gib=1)), encoding='utf-8')
            def sample(parameters, roughness, *args):
                if parameters.pml_strength == 2.5 and roughness.seed == 8:
                    raise SampleLimitError(folder/'resources.json', 'stopped_memory_limit')
                return folder
            cache = Path(directory)/'shared_cache'
            with patch('fullwave_prototype.haze_convergence._sample', side_effect=sample) as run, patch(
                'fullwave_prototype.haze_convergence.radiation_amplitude',
                side_effect=lambda points, weights, directions, p: np.ones((len(directions), 2))):
                report = run_plan(plan, Path(directory)/'result', sample_cache=cache)
            self.assertEqual(run.call_args.args[2], cache.resolve())
            self.assertEqual(report['sample_cache'], str(cache.resolve()))
            self.assertEqual(report['status'], 'complete_with_limits')
            self.assertEqual(len(report['axes']['pml_strength']['cases']), 1)
            self.assertEqual(report['axes']['pml_strength']['comparison']['status'], 'resource_limited')
            self.assertEqual(report['axes']['pml_strength']['limits'][0]['level'], 2.5)


if __name__ == '__main__':
    unittest.main()

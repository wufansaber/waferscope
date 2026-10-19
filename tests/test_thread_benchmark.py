"""A fast trial must never be recommended when the numerical result changes."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from fullwave_prototype.thread_benchmark import summarize


class ThreadBenchmarkTests(unittest.TestCase):
    def trial(self, root, repeat, threads, seconds, perturbation=0, status='completed'):
        folder = Path(root)/'trial_{}_{}t'.format(repeat, threads)
        folder.mkdir()
        (folder/'resources.json').write_text(json.dumps(dict(status=status, peak_accounted_gib=8)))
        if status != 'completed':
            return
        (folder/'result.json').write_text(json.dumps(dict(parameters=dict(threads=threads),
            metrics=dict(elements=10, dofs=30, factored_dofs=20, solve_seconds=seconds))))
        np.savez(folder/'fields.npz', e=np.ones(3)+perturbation,
                 h=np.ones(3), amplitude=np.ones(2))

    def test_select_fewer_threads_in_timing_noise_band(self):
        with tempfile.TemporaryDirectory() as root:
            for repeat in (1, 2):
                self.trial(root, repeat, 2, 100)
                self.trial(root, repeat, 4, 60)
                self.trial(root, repeat, 8, 59)
            report = summarize(root)
            self.assertTrue(report['correctness_passed'])
            self.assertEqual(report['recommended_threads'], 4)

    def test_reject_inaccurate_or_failed_fast_trials(self):
        for failed in (False, True):
            with tempfile.TemporaryDirectory() as root:
                for repeat in (1, 2):
                    self.trial(root, repeat, 2, 100)
                    self.trial(root, repeat, 4, 30, perturbation=.001,
                               status='stopped_memory_limit' if failed else 'completed')
                report = summarize(root)
                self.assertFalse(report['correctness_passed'])
                self.assertIsNone(report['recommended_threads'])

    def test_no_recommendation_for_one_run(self):
        with tempfile.TemporaryDirectory() as root:
            self.trial(root, 1, 2, 100)
            self.assertIsNone(summarize(root)['recommended_threads'])


if __name__ == '__main__':
    unittest.main()

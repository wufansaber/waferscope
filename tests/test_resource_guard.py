import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

from fullwave_prototype.guarded_run import _wait_for_owned, stop_reason, run_guarded


class GuardPolicyTests(unittest.TestCase):
    def test_limits(self):
        self.assertIsNone(stop_reason(1, 5, 10, 4, 2, 20))
        self.assertEqual(stop_reason(5, 5, 10, 4, 2, 20), 'memory_limit')
        self.assertEqual(stop_reason(1, 1, 10, 4, 2, 20), 'system_memory_reserve')
        self.assertEqual(stop_reason(1, 5, 30, 4, 2, 20), 'timeout')


@unittest.skipIf(importlib.util.find_spec('psutil') is None, 'Resource monitor not installed')
class GuardProcessTests(unittest.TestCase):
    def test_exited_windows_child_is_not_waited_on_again(self):
        import psutil
        dead = Mock(pid=1234)
        dead.is_running.return_value = False
        with patch('psutil.wait_procs', side_effect=psutil.AccessDenied(1234)) as wait:
            self.assertEqual(_wait_for_owned([dead]), [])
            wait.assert_not_called()

    def test_access_denied_wait_rechecks_and_keeps_unknown_process(self):
        import psutil
        gone = Mock(pid=1234)
        gone.is_running.side_effect = [True, False]
        unknown = Mock(pid=5678)
        unknown.is_running.side_effect = [True, psutil.AccessDenied(5678)]
        with patch('psutil.wait_procs', side_effect=psutil.AccessDenied(1234)), \
             patch('fullwave_prototype.guarded_run.time.sleep'):
            self.assertEqual(_wait_for_owned([gone, unknown]), [5678])

    def test_cancellation_cleanup(self):
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder)/'report.json'
            cancel = Path(folder)/'cancel.request'
            cancel.touch()
            result = run_guarded([sys.executable, '-c', 'import time; time.sleep(30)'],
                                 target, reserve=.01, timeout=30, cancel_file=cancel)
            self.assertEqual(result['status'], 'stopped_cancelled')
            self.assertEqual(result['surviving_owned_pids'], [])

    def test_completion_and_timeout_cleanup(self):
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder)/'report.json'
            result = run_guarded([sys.executable, '-c', 'print(42)'], target, reserve=.01, timeout=30)
            self.assertEqual(result['status'], 'completed')
            result = run_guarded([sys.executable, '-c', 'import time; time.sleep(30)'],
                                 target, reserve=.01, timeout=1)
            self.assertEqual(result['status'], 'stopped_timeout')
            self.assertEqual(result['surviving_owned_pids'], [])
            self.assertEqual(json.loads(target.read_text())['status'], 'stopped_timeout')


if __name__ == '__main__':
    unittest.main()

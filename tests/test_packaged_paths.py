import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fullwave_ui.runtime_paths import application_root, default_engine, results_root


class PackagedPathsTests(unittest.TestCase):
    def test_engine_prefers_adjacent_bundle(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.assertEqual(default_engine(root), root/'build/ngsolve_env/Scripts/python.exe')
            engine = root/'fem/FarField_FEM.exe'
            engine.parent.mkdir()
            engine.touch()
            self.assertEqual(default_engine(root), engine)

    def test_frozen_root_and_writable_results(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            with patch('sys.frozen', True, create=True), patch('sys.executable', str(root/'app.exe')):
                self.assertEqual(application_root(), root)
                with patch.dict(os.environ, {'LOCALAPPDATA': str(root/'user'), 'FARFIELD_RESULTS_DIR': ''}):
                    self.assertEqual(results_root(root), root/'user/FarField_Simulation/results')

    def test_results_override(self):
        with patch.dict(os.environ, {'FARFIELD_RESULTS_DIR': 'build/package_test'}):
            self.assertEqual(results_root(Path('.')), Path('build/package_test').resolve())

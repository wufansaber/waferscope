"""Qt bridge checks without importing NGSolve or starting a heavy solve."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
HAS_QT = importlib.util.find_spec('PyQt5') is not None
if HAS_QT:
    from PyQt5.QtCore import QProcess
    from PyQt5.QtGui import QPixmap
    from PyQt5.QtWidgets import QApplication
    from fullwave_ui.controls import SimulationControls
    from fullwave_ui.results import ResultsTabs
    from fullwave_ui.widget import FullwaveWidget
    from fullwave_prototype.imaging_run import build_parser, scattering_settings


@unittest.skipUnless(HAS_QT, 'Qt environment required')
class FullwaveUiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_controls_forward_materials_channels_and_polarization(self):
        controls = SimulationControls()
        controls.fields['particle-k'].setValue(.12)
        controls.fields['wavelength-nm'].setValue(266)
        controls.fields['polarization'].setCurrentText('s')
        controls.fields['na2'].setValue(0)
        controls.fields['tilt1'].setValue(30)
        args = build_parser().parse_args(controls.arguments())
        p, _ = scattering_settings(args)
        self.assertEqual((args.accuracy, args.pupil_samples, p.order), ('balanced', 128, 3))
        self.assertEqual((p.pml_strength, p.pml_mesh_scale), (1., 1.25))
        self.assertEqual((p.pml_order, p.pml_inner_fraction), (0, .5))
        self.assertEqual((p.particle_k, p.wavelength_nm, p.polarization), (.12, 266, 's'))
        self.assertEqual(args.na, [.5, 0, .5])
        self.assertEqual(args.tilts, [30, -45, 0])
        controls.fields['pml-order'].setCurrentText('2')
        controls.fields['pml-inner-fraction'].setValue(.25)
        controls.fields['local-mesh-scale'].setValue(.9)
        tuned, rough = scattering_settings(build_parser().parse_args(controls.arguments()))
        self.assertEqual(tuned.pml_order, 2)
        self.assertEqual(tuned.pml_inner_fraction, .25)
        self.assertAlmostEqual(tuned.particle_mesh_um, .0135)
        self.assertAlmostEqual(rough.mesh_um, .0225)
        controls.fields['gap-nm'].setValue(0)
        with self.assertRaises(ValueError):
            controls.arguments()

    def test_load_is_atomic_and_keeps_result_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            report = dict(channels=[{'with': dict(image_power_w=1e-12*i)}
                                    for i in (1, 2, 3)])
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            tabs = ResultsTabs()
            pixmap = QPixmap(80, 60)
            pixmap.fill()
            for _, name in tabs.FIGURES:
                self.assertTrue(pixmap.save(str(path/name)))
            self.assertEqual(tabs.load(path), [1, 2, 3])
            old = tabs.parameters.toPlainText()
            (path/tabs.FIGURES[2][1]).unlink()
            with self.assertRaises(ValueError):
                tabs.load(path)
            self.assertEqual(tabs.parameters.toPlainText(), old)
            self.assertIn(str(path.resolve()), old)
            report['channels'][0]['with']['image_power_w'] = float('nan')
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'invalid image powers'):
                tabs.load(path)

    def test_process_command_and_cooperative_cancel(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            python = root/'python.exe'
            python.touch()
            widget = FullwaveWidget()
            widget.interpreter.setText(str(python))
            with patch('fullwave_ui.widget.ROOT', root), patch.object(widget.process, 'start') as start:
                widget.start_job()
                interpreter, command = start.call_args[0]
                self.assertEqual(interpreter, str(python))
                self.assertIn('fullwave_prototype.guarded_run', command)
                self.assertIn('--cancel-file', command)
                self.assertIn('--particle-k', command)
                self.assertFalse(widget.controls.isEnabled())
                with patch.object(widget.process, 'state', return_value=QProcess.Running):
                    widget.cancel_job()
                self.assertTrue(widget.cancel_file.exists())
                self.assertFalse(widget.cancel.isEnabled())
                widget.process_error(QProcess.FailedToStart)
                self.assertTrue(widget.controls.isEnabled())

    def test_window_keeps_scatmech_and_fem_separate(self):
        from scatmech_ui.app import ScatmechStandaloneWindow
        window = ScatmechStandaloneWindow()
        self.assertIs(window.tabs.widget(0), window.simulation_widget)
        self.assertIs(window.tabs.widget(1), window.fullwave_widget)
        self.assertTrue(window.fullwave_widget.wait_for_shutdown())
        window.close()


if __name__ == '__main__':
    unittest.main()

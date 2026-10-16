"""Saved-source selection and explicit, independent-channel image views."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fullwave_ui.saved_source import reusable_source

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
HAS_QT = importlib.util.find_spec('PyQt5') is not None
if HAS_QT:
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtGui import QPixmap
    from fullwave_ui.widget import FullwaveWidget
    from fullwave_ui.results import ResultsTabs


def source_at(path):
    path.mkdir(parents=True, exist_ok=True)
    (path/'report.json').write_text('{"source_format_version": 1}', encoding='utf-8')
    (path/'paired_fields.npz').touch()
    return path.resolve()


class SourceTests(unittest.TestCase):
    def test_portable_local_source_precedes_recorded_original(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original = source_at(root/'original')
            local = source_at(root/'moved'/'scattering')
            (root/'moved'/'imaging_report.json').write_text(json.dumps(
                dict(source_directory=str(original))), encoding='utf-8')
            self.assertEqual(reusable_source(root/'moved'), local)
            self.assertEqual(reusable_source(original), original)

    def test_reimage_result_uses_recorded_source_and_missing_is_disabled(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = source_at(root/'saved')
            report = root/'imaging_report.json'
            report.write_text('{"source_directory": "saved"}', encoding='utf-8')
            self.assertEqual(reusable_source(root), source)
            (source/'paired_fields.npz').unlink()
            self.assertIsNone(reusable_source(root))
            report.write_text('[]', encoding='utf-8')
            self.assertIsNone(reusable_source(root))


@unittest.skipUnless(HAS_QT, 'Qt environment required')
class ImagingEntryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def result_at(self, path):
        report = dict(channels=[{'with': dict(image_power_w=i*1e-12)} for i in (1, 2, 3)],
                      image_channel_figures=[1, 2, 3])
        (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
        pixmap = QPixmap(80, 60)
        pixmap.fill()
        for name in [n for _, n in ResultsTabs.FIGURES]+['image_channel_{}.png'.format(i) for i in (1, 2, 3)]:
            self.assertTrue(pixmap.save(str(path/name)))
        return report

    def test_current_image_button_uses_loaded_source_and_respects_busy(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            self.result_at(path)
            source = source_at(path/'scattering')
            widget = FullwaveWidget()
            self.assertFalse(widget.image_current_button.isEnabled())
            self.assertTrue(widget.load_result(path))
            self.assertEqual(widget.loaded_source, source)
            self.assertTrue(widget.image_current_button.isEnabled())
            with patch.object(widget, 'start_job') as start:
                widget.image_current_button.click()
                start.assert_called_once_with(source=source)
            widget.set_running(True)
            self.assertFalse(widget.image_current_button.isEnabled())
            widget.set_running(False)
            self.assertTrue(widget.image_current_button.isEnabled())
            (source/'paired_fields.npz').unlink()
            self.assertTrue(widget.load_result(path))
            self.assertFalse(widget.image_current_button.isEnabled())
            widget.close()

    def test_individual_images_atomic_load_and_legacy_clear(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            report = self.result_at(path)
            tabs = ResultsTabs()
            self.assertEqual(tabs.load(path), [1., 2., 3.])
            for view in tabs.image_views.values():
                self.assertIsNotNone(view.item)
                self.assertTrue(tabs.image_tabs.isTabEnabled(tabs.image_tabs.indexOf(view)))
            old = tabs.parameters.toPlainText()
            (path/'image_channel_2.png').unlink()
            with self.assertRaisesRegex(ValueError, 'image channel figure: 2'):
                tabs.load(path)
            self.assertEqual(tabs.parameters.toPlainText(), old)
            report.pop('image_channel_figures')
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            tabs.load(path)
            for view in tabs.image_views.values():
                self.assertIsNone(view.item)
                self.assertFalse(tabs.image_tabs.isTabEnabled(tabs.image_tabs.indexOf(view)))
            tabs.close()

    def test_paired_camera_metrics_atomic_loading_and_legacy_clear(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            report = self.result_at(path)
            report.update(settings=dict(detector=dict(enabled=True, reference_frames=4)),
                          camera_detection_model_version=1)
            for channel in report['channels']:
                channel['detector'] = dict(detection_roi_valid=True, detection_roi_snr=2.,
                                           detection_template_snr=3.)
            pixmap = QPixmap(80, 60)
            pixmap.fill()
            for name in ('detector.png', 'detector_components.png', 'detector_snr.png'):
                pixmap.save(str(path/name))
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            tabs = ResultsTabs()
            tabs.load(path)
            roi_row = [key for _, key in tabs.camera_metrics.ROWS].index('detection_roi_snr')
            self.assertEqual(tabs.camera_metrics.table.item(roi_row, 0).text(), '2')
            self.assertTrue(tabs.camera_tabs.isTabEnabled(1))
            (path/'detector_snr.png').unlink()
            with self.assertRaisesRegex(ValueError, 'paired camera figure'):
                tabs.load(path)
            self.assertEqual(tabs.camera_metrics.table.item(roi_row, 0).text(), '2')
            report.pop('camera_detection_model_version')
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            tabs.load(path)
            self.assertFalse(tabs.camera_tabs.isTabEnabled(1))
            self.assertIsNone(tabs.camera_metrics.table.item(roi_row, 0))
            tabs.close()


if __name__ == '__main__':
    unittest.main()

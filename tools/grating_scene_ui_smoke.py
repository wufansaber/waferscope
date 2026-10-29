"""Offscreen UI workflow check using an actual completed native scene."""
import argparse
import json
import os
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def main():
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtGui import QFont, QFontDatabase
    from grating_ui.widget import GratingWidget
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--camera', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    font_path = Path('C:/Windows/Fonts/segoeui.ttf')
    if font_path.is_file():
        QFontDatabase.addApplicationFont(str(font_path))
        app.setFont(QFont('Segoe UI', 10))
    widget = GratingWidget()
    widget.resize(1500, 940)
    widget.show()
    if not widget.load_result(args.scene):
        raise RuntimeError('Completed native scene could not be loaded by Qt')
    app.processEvents()
    for index, name in ((0, 'geometry'), (2, 'near_xz'), (6, 'images'), (7, 'powers'), (9, 'diagnostics')):
        widget.results.setCurrentIndex(index)
        app.processEvents()
        widget.grab().save(str(args.output/(name+'.png')))
    if args.camera:
        with patch('grating_ui.widget.QFileDialog.getExistingDirectory', return_value=str(args.camera)):
            widget.choose_result()
        if widget.camera_result != args.camera.resolve():
            raise RuntimeError('Camera folder could not be loaded by the actual UI loader')
        app.processEvents()
        widget.grab().save(str(args.output/'camera.png'))
        widget.results.setCurrentIndex(11)
        app.processEvents()
        widget.grab().save(str(args.output/'camera_metrics.png'))
    widget.resize(1200, 760)
    widget.results.setCurrentIndex(6)
    app.processEvents()
    widget.grab().save(str(args.output/'compact.png'))
    report = dict(loaded=str(widget.loaded_result), tabs=widget.results.count(),
        na=[c.value() for c in widget.controls.nas], gap_nm=widget.controls.scene_fields['gap_um'][0].value(),
        update_optics_enabled=widget.reimage.isEnabled(), export_enabled=widget.export.isEnabled())
    (args.output/'ui_smoke.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    widget.close()
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

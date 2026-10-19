"""Capture current Qt controls and saved-result views for the operation manual."""
import json
import os
from pathlib import Path

os.environ['QT_QPA_PLATFORM'] = 'offscreen'


def main():
    from PyQt5.QtCore import QPoint
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication, QScrollArea, QSplitter
    from scatmech_ui.app import ScatmechStandaloneWindow
    root = Path(__file__).resolve().parents[1]
    out = root/'far_field_standalone/manual_assets'
    out.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    QFontDatabase.addApplicationFont('C:/Windows/Fonts/segoeui.ttf')
    app.setFont(QFont('Segoe UI', 10))
    window = ScatmechStandaloneWindow()
    window.resize(1560, 1000)
    window.tabs.setCurrentIndex(1)
    widget = window.fullwave_widget
    window.show()
    assert widget.load_result(root/'build/channel_symmetry_20261016/low_haze_02nm')
    widget.results.setCurrentWidget(widget.results.image_tabs)
    app.processEvents()
    widget.results.fit_current()
    assert window.grab().save(str(out/'fem_overview.png'))

    camera = root/'build/channel_symmetry_20261016/packaged_camera/results/20261016_174505_240531'
    assert widget.load_result(camera)
    widget.results.setCurrentWidget(widget.results.camera_tabs)
    widget.results.camera_tabs.setCurrentWidget(widget.results.camera_metrics)
    app.processEvents()
    assert widget.results.camera_metrics.grab().save(str(out/'camera_metrics.png'))
    scroll = widget.findChildren(QScrollArea)[0]
    scroll.setMinimumWidth(550)
    window.findChildren(QSplitter)[0].setSizes([560, 950])
    app.processEvents()
    fields = widget.controls.fields
    scroll.verticalScrollBar().setValue(max(0, fields['detector'].mapTo(scroll.widget(), QPoint(0,0)).y()-25))
    app.processEvents()
    # Crop a readable control excerpt; the manual table documents all fields.
    viewport = scroll.viewport()
    location = fields['detector'].mapTo(viewport, fields['detector'].rect().topLeft())
    top = max(0, location.y()-28)
    end = fields['detector-fill-factor'].mapTo(viewport, fields['detector-fill-factor'].rect().bottomLeft()).y()+8
    assert viewport.grab().copy(0, top, viewport.width(), min(end-top, viewport.height()-top)).save(str(out/'camera_controls.png'))
    widget.controls.preset.setCurrentText('266 nm cross preview')
    widget.controls.preset_button.click()
    scroll.verticalScrollBar().setValue(max(0, fields['defect-shape'].mapTo(scroll.widget(), QPoint(0,0)).y()-25))
    app.processEvents()
    top = max(0, fields['defect-shape'].mapTo(viewport, fields['defect-shape'].rect().topLeft()).y()-25)
    end = fields['gap-nm'].mapTo(viewport, fields['gap-nm'].rect().bottomLeft()).y()+8
    assert viewport.grab().copy(0, top, viewport.width(), min(end-top, viewport.height()-top)).save(str(out/'pattern_controls.png'))
    snapshots = {}
    for label in ('10 nm PSL preview', '266 nm low-haze camera', '266 nm cross preview'):
        widget.controls.preset.setCurrentText(label)
        widget.controls.preset_button.click()
        snapshots[label] = {name: (f.currentText() if hasattr(f, 'currentText') else f.isChecked()
            if hasattr(f, 'isChecked') else f.value()) for name, f in fields.items()}
    evidence = root/'build/manual_pdf_20261017'
    evidence.mkdir(parents=True, exist_ok=True)
    (evidence/'ui_inputs.json').write_text(json.dumps(snapshots, indent=2), encoding='utf-8')
    print('Captured four current UI views and three preset snapshots')
    window.close()


if __name__ == '__main__':
    main()

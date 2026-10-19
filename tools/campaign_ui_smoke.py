"""Inspect the completed ensemble and finite-envelope controls in the real Qt panel."""
import argparse
import os
from pathlib import Path
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def main():
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication, QScrollArea
    from fullwave_ui.widget import FullwaveWidget
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True)
    args = parser.parse_args()
    root = Path(args.root)
    out = root/'ui'
    out.mkdir(parents=True, exist_ok=True)
    app = QApplication([])
    QFontDatabase.addApplicationFont('C:/Windows/Fonts/segoeui.ttf')
    app.setFont(QFont('Segoe UI', 9))
    widget = FullwaveWidget()
    widget.resize(1480, 900)
    widget.show()
    assert widget.load_whole_haze_result(root/'ensemble128')
    assert widget.loaded_whole is not None and widget.loaded_hybrid is None
    widget.refresh_validation_actions()
    app.processEvents()
    widget.results.fit_current()
    assert widget.grab().save(str(out/'ensemble128_ui.png'))
    widget.controls.fields['beam-definition'].setCurrentText('measured_wafer')
    widget.controls.beam_measurement_path.setText(str((root/'beam_inputs/synthetic_intensity_demo.npz').resolve()))
    assert widget.controls.fields['beam-phase-assumption'].isEnabled()
    assert 'Invalid' not in widget.controls.irradiance_label.text()
    scroll = widget.controls.parentWidget()
    while scroll is not None and not isinstance(scroll, QScrollArea):
        scroll = scroll.parentWidget()
    if scroll is not None:
        scroll.ensureWidgetVisible(widget.controls.beam_measurement_path)
    app.processEvents()
    assert widget.grab().save(str(out/'measured_beam_controls.png'))
    widget.resize(1280, 800)
    app.processEvents()
    if scroll is not None:
        scroll.ensureWidgetVisible(widget.controls.beam_measurement_path)
    app.processEvents()
    widget.results.fit_current()
    assert widget.grab().save(str(out/'measured_beam_controls_compact.png'))
    widget.close()
    print('128-surface result and measured-beam controls passed')


if __name__ == '__main__':
    main()

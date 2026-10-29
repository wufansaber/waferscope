"""Load actual continuous radiation into the source Qt page and save a screenshot."""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtGui import QFont
    from grating_ui.widget import GratingWidget
    app = QApplication.instance() or QApplication([])
    app.setStyle('Fusion')
    app.setFont(QFont('Arial', 10))
    widget = GratingWidget()
    widget.resize(1600, 950)
    assert widget.open_result(args.result)
    widget.show()
    app.processEvents()
    args.output.mkdir(parents=True, exist_ok=True)
    assert widget.grab().save(str(args.output/'continuous_ui.png'))
    widget.results.continuous_view.surface.setChecked(True)
    app.processEvents()
    assert widget.grab().save(str(args.output/'continuous_ui_3d.png'))
    (args.output/'ui_smoke.json').write_text(json.dumps(dict(
        result=str(args.result.resolve()), opened=True, two_d_and_three_d_rendered=True,
        export_enabled=widget.export.isEnabled(), tab=widget.results.tabText(widget.results.currentIndex())),
        indent=2), encoding='utf-8')
    assert widget.wait_for_shutdown()
    widget.close()


if __name__ == '__main__':
    main()

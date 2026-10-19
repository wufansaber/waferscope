"""Offscreen source-window rendering and optional real guarded cancellation check."""
import argparse
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', default='build/fullwave_imaging/complete_preview')
    parser.add_argument('--output', default='build/followup_validation')
    parser.add_argument('--cancel-test', action='store_true')
    args = parser.parse_args()
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt5.QtCore import QTimer
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication, QScrollArea
    from scatmech_ui.app import ScatmechStandaloneWindow
    app = QApplication.instance() or QApplication([])
    # Offscreen Qt on Windows does not populate its normal system font database.
    font = Path(os.environ.get('WINDIR', 'C:/Windows'))/'Fonts/segoeui.ttf'
    if font.exists():
        QFontDatabase.addApplicationFont(str(font))
        app.setFont(QFont('Segoe UI', 9))
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    window = ScatmechStandaloneWindow()
    window.tabs.setCurrentIndex(1)
    widget = window.fullwave_widget
    widget.load_result(args.result)
    window.show()
    for width, height, tab in ((1200, 800, 0), (1680, 980, 2)):
        window.resize(width, height)
        widget.results.setCurrentIndex(tab)
        app.processEvents()
        if not window.grab().save(str(output/('fem_ui_'+str(width)+'.png'))):
            raise RuntimeError('Screenshot export failed')
        for scroll in widget.findChildren(QScrollArea):
            if scroll.widget() and scroll.widget().isAncestorOf(widget.controls):
                scroll.verticalScrollBar().setValue(scroll.verticalScrollBar().maximum())
                app.processEvents()
                if not window.grab().save(str(output/('fem_controls_'+str(width)+'.png'))):
                    raise RuntimeError('Controls screenshot export failed')
                scroll.verticalScrollBar().setValue(0)
    if args.cancel_test:
        widget.process.finished.connect(app.quit)
        widget.start_job()
        QTimer.singleShot(1000, widget.cancel_job)
        QTimer.singleShot(30000, app.quit)
        app.exec_()
        if not widget.wait_for_shutdown():
            raise RuntimeError('Guard did not stop its owned job')
        report = json.loads((widget.output/'resources.json').read_text(encoding='utf-8'))
        if report['status'] != 'stopped_cancelled' or report['surviving_owned_pids']:
            raise RuntimeError('Unexpected cancellation result: '+str(report))
        print('Real UI cancellation verified:', widget.output)
    window.close()
    print('Offscreen UI checks complete:', output.resolve())


if __name__ == '__main__':
    main()

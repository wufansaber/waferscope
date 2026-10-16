"""Offscreen verification of saved paired-camera tables at two window sizes."""
import argparse
import os
from pathlib import Path


def run(source, output):
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtGui import QFont, QFontDatabase
    from .results import ResultsTabs
    app = QApplication.instance() or QApplication([])
    font = Path(os.environ.get('WINDIR', 'C:/Windows'))/'Fonts/segoeui.ttf'
    if font.is_file():
        QFontDatabase.addApplicationFont(str(font))
        app.setFont(QFont('Segoe UI', 9))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    tabs = ResultsTabs()
    tabs.load(source)
    tabs.setCurrentWidget(tabs.camera_tabs)
    tabs.camera_tabs.setCurrentWidget(tabs.camera_metrics)
    for width in (1200, 640):
        tabs.resize(width, 800)
        tabs.show()
        app.processEvents()
        if not tabs.grab().save(str(output/'camera_metrics_{}.png'.format(width))):
            raise RuntimeError('Unable to save camera UI screenshot')
    tabs.camera_tabs.setCurrentWidget(tabs.camera_components)
    tabs.resize(1200, 900)
    app.processEvents()
    if not tabs.grab().save(str(output/'camera_components_ui.png')):
        raise RuntimeError('Unable to save camera decomposition screenshot')
    tabs.close()
    app.processEvents()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    run(args.source, args.output)

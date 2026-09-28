"""Aspect-preserving scientific result views and saved-parameter provenance."""
import json
import math
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (QGraphicsScene, QGraphicsView, QPlainTextEdit,
                            QTabWidget)


class FigureView(QGraphicsView):
    def __init__(self):
        super().__init__()
        self.setScene(QGraphicsScene(self))
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setToolTip('Result figure')
        self.item = None

    def load(self, path):
        pixmap = QPixmap(str(path))
        if pixmap.isNull():
            raise ValueError('Cannot load figure: '+str(path))
        self.scene().clear()
        self.item = self.scene().addPixmap(pixmap)
        self.item.setTransformationMode(Qt.SmoothTransformation)
        self.scene().setSceneRect(self.item.boundingRect())
        self.fit()

    def fit(self):
        if self.item is not None:
            self.fitInView(self.item, Qt.KeepAspectRatio)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.fit()

    def wheelEvent(self, event):
        self.scale(1.2 if event.angleDelta().y() > 0 else 1/1.2,
                   1.2 if event.angleDelta().y() > 0 else 1/1.2)
        event.accept()


class ResultsTabs(QTabWidget):
    FIGURES = (('Near field', 'near_field.png'), ('Far field', 'far_field.png'),
               ('Three image planes', 'image_planes.png'),
               ('Coherent decomposition', 'image_decomposition.png'))

    def __init__(self):
        super().__init__()
        self.views = []
        for title, _ in self.FIGURES:
            view = FigureView()
            self.views.append(view)
            self.addTab(view, title)
        self.parameters = QPlainTextEdit()
        self.parameters.setReadOnly(True)
        self.addTab(self.parameters, 'Result parameters')
        self.detector_view = FigureView()
        self.addTab(self.detector_view, 'Camera')
        self.setTabEnabled(self.indexOf(self.detector_view), False)
        self.currentChanged.connect(self.fit_current)

    def fit_current(self, *_):
        view = self.currentWidget()
        if isinstance(view, FigureView):
            view.fit()

    def load(self, directory):
        directory = Path(directory)
        report = json.loads((directory/'imaging_report.json').read_text(encoding='utf-8'))
        channels = report['channels']
        if len(channels) != 3:
            raise ValueError('Result must contain three channels')
        powers = [float(c['with']['image_power_w'])*1e12 for c in channels]
        if not all(math.isfinite(power) and power >= 0 for power in powers):
            raise ValueError('Result contains invalid image powers')
        # Validate every image before replacing the current result.
        detector = report.get('settings', {}).get('detector', {}).get('enabled', False)
        if detector and QPixmap(str(directory/'detector.png')).isNull():
            raise ValueError('Missing or invalid detector figure')
        for _, name in self.FIGURES:
            if QPixmap(str(directory/name)).isNull():
                raise ValueError('Missing or invalid result figure: '+name)
        for view, (_, name) in zip(self.views, self.FIGURES):
            view.load(directory/name)
        if detector:
            self.detector_view.load(directory/'detector.png')
        else:
            self.detector_view.scene().clear()
            self.detector_view.item = None
        self.setTabEnabled(self.indexOf(self.detector_view), detector)
        self.parameters.setPlainText('Loaded result: '+str(directory.resolve())+'\n\n'+
                                     json.dumps(report, indent=2, ensure_ascii=False))
        return powers

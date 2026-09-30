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
            if not self.count():
                self.near_tabs = QTabWidget()
                self.near_tabs.addTab(view, 'Overview')
                self.near_views = {}
                for plane in ('xz', 'xy', 'yz'):
                    self.near_views[plane] = FigureView()
                    index = self.near_tabs.addTab(self.near_views[plane], plane.upper())
                    self.near_tabs.setTabEnabled(index, False)
                self.near_tabs.currentChanged.connect(self.fit_current)
                self.addTab(self.near_tabs, title)
            elif title == 'Three image planes':
                self.image_tabs = QTabWidget()
                self.image_tabs.addTab(view, 'Overview')
                self.image_views = {}
                for channel in (1, 2, 3):
                    self.image_views[channel] = FigureView()
                    index = self.image_tabs.addTab(self.image_views[channel], 'Ch{}'.format(channel))
                    self.image_tabs.setTabEnabled(index, False)
                self.image_tabs.currentChanged.connect(self.fit_current)
                self.addTab(self.image_tabs, title)
            else:
                self.addTab(view, title)
        self.comparison_tabs = QTabWidget()
        self.comparison_views = {}
        for key, title in (('far', 'Far field'), ('far_log', 'Far log'),
                           ('image', 'Image plane'), ('image_log', 'Image log')):
            view = FigureView()
            self.comparison_views[key] = view
            self.comparison_tabs.addTab(view, title)
        self.comparison_totals = QPlainTextEdit()
        self.comparison_totals.setReadOnly(True)
        self.comparison_tabs.addTab(self.comparison_totals, 'Totals')
        self.addTab(self.comparison_tabs, 'Particle / haze')
        self.setTabEnabled(self.indexOf(self.comparison_tabs), False)
        self.parameters = QPlainTextEdit()
        self.parameters.setReadOnly(True)
        self.addTab(self.parameters, 'Result parameters')
        self.detector_view = FigureView()
        self.addTab(self.detector_view, 'Camera')
        self.setTabEnabled(self.indexOf(self.detector_view), False)
        self.currentChanged.connect(self.fit_current)

    def fit_current(self, *_):
        view = self.currentWidget()
        if isinstance(view, QTabWidget):
            view = view.currentWidget()
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
        planes = {info['plane'].lower() for info in report.get('near_field_planes', [])}
        if not planes.issubset(self.near_views):
            raise ValueError('Unknown near-field plane')
        for plane in planes:
            if QPixmap(str(directory/('near_field_'+plane+'.png'))).isNull():
                raise ValueError('Missing or invalid near-field figure: '+plane)
        image_channels = report.get('image_channel_figures', [])
        if image_channels not in ([], [1, 2, 3]):
            raise ValueError('Invalid image channel figure list')
        for channel in image_channels:
            if QPixmap(str(directory/('image_channel_{}.png'.format(channel)))).isNull():
                raise ValueError('Missing or invalid image channel figure: '+str(channel))
        comparison = None
        manifest_path = directory/'component_comparison_manifest.json'
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
            comparison_dir = (directory/manifest['directory']).resolve()
            comparison_dir.relative_to(directory.resolve())
            comparison_report = json.loads((comparison_dir/manifest['report']).read_text(encoding='utf-8'))
            if comparison_report.get('status') != 'complete':
                raise ValueError('Incomplete component comparison')
            figures = manifest['figures']
            for key in self.comparison_views:
                name = figures[key]
                if Path(name).name != name or QPixmap(str(comparison_dir/name)).isNull():
                    raise ValueError('Missing or invalid component comparison figure: '+key)
            totals = manifest['totals']
            if Path(totals).name != totals:
                raise ValueError('Invalid component comparison totals path')
            comparison = (comparison_dir, figures,
                          (comparison_dir/totals).read_text(encoding='utf-8'))
        for view, (_, name) in zip(self.views, self.FIGURES):
            view.load(directory/name)
        self.near_tabs.setCurrentIndex(0)
        for plane, view in self.near_views.items():
            if plane in planes:
                view.load(directory/('near_field_'+plane+'.png'))
            else:
                view.scene().clear()
                view.item = None
            self.near_tabs.setTabEnabled(self.near_tabs.indexOf(view), plane in planes)
        self.image_tabs.setCurrentIndex(0)
        for channel, view in self.image_views.items():
            if channel in image_channels:
                view.load(directory/('image_channel_{}.png'.format(channel)))
            else:
                view.scene().clear()
                view.item = None
            self.image_tabs.setTabEnabled(self.image_tabs.indexOf(view), channel in image_channels)
        if comparison is not None:
            comparison_dir, figures, totals = comparison
            for key, view in self.comparison_views.items():
                view.load(comparison_dir/figures[key])
            self.comparison_totals.setPlainText(totals)
        else:
            for view in self.comparison_views.values():
                view.scene().clear()
                view.item = None
            self.comparison_totals.clear()
        self.setTabEnabled(self.indexOf(self.comparison_tabs), comparison is not None)
        if detector:
            self.detector_view.load(directory/'detector.png')
        else:
            self.detector_view.scene().clear()
            self.detector_view.item = None
        self.setTabEnabled(self.indexOf(self.detector_view), detector)
        self.parameters.setPlainText('Loaded result: '+str(directory.resolve())+'\n\n'+
                                     json.dumps(report, indent=2, ensure_ascii=False))
        return powers

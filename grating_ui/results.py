"""Lazy Qt figure tabs, preserving unnormalized optical values."""
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QHeaderView, QScrollArea, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget

from grating_scene.plotting import geometry_figure, near_figure, far_figure, image_figure
from grating_scene.diagnostics import scene_diagnostics
from .camera_view import CameraView
from .continuous_view import ContinuousView


class SceneResults(QTabWidget):
    def __init__(self):
        super().__init__()
        self.data = None
        self.logarithmic = True
        self.figures, self.canvases = {}, {}
        for key, title in (('geometry', 'Geometry'), ('xy', 'Near XY'), ('xz', 'Near XZ'),
            ('yz', 'Near YZ'), ('far', 'Far field 2D'), ('far3d', 'Far field 3D'), ('images', 'Images')):
            panel = QWidget()
            layout = QVBoxLayout(panel)
            figure = Figure(figsize=(12, 7))
            canvas = FigureCanvasQTAgg(figure)
            canvas.setMinimumSize(1020, 820 if key == 'images' else 450)
            layout.addWidget(NavigationToolbar2QT(canvas, panel))
            layout.addWidget(canvas)
            self.figures[key], self.canvases[key] = figure, canvas
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(panel)
            self.addTab(scroll, title)
        self.table = QTableWidget(4, 8)
        self.table.setHorizontalHeaderLabels(('CH', 'NA', 'Tilt (deg)', 'Orders',
            'Reference (W)', 'Sample (W)', 'Sample-ref (W)', 'Difference field (W)'))
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.addTab(self.table, 'Collected powers')
        self.validation = QTableWidget(0, 6)
        self.validation.setHorizontalHeaderLabels(('Case', 'State', 'CH1 field (%)', 'CH2 field (%)',
                                                   'CH3 field (%)', 'Local step'))
        self.validation.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.addTab(self.validation, 'Validation')
        self.diagnostics = QTableWidget(0, 4)
        self.diagnostics.setHorizontalHeaderLabels(('Diagnostic', 'Reference', 'Sample', 'Difference field'))
        self.diagnostics.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.addTab(self.diagnostics, 'Diagnostics')
        self.camera_view = CameraView()
        self.camera_canvas = self.camera_view.canvas
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.camera_view)
        self.addTab(scroll, 'Camera DN')
        self.camera_metrics = QTableWidget(3, 6)
        self.camera_metrics.setHorizontalHeaderLabels(('CH', 'Raw DN ROI SNR', 'Gray pairs', 'Gray status',
            'Expected overfull pixels', 'Periodic copies X / Y'))
        self.camera_metrics.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.addTab(self.camera_metrics, 'Camera metrics')
        self.addTab(self.camera_view.table, 'Camera brightness')
        self.continuous_view = ContinuousView()
        self.addTab(self.continuous_view, 'Continuous angles')
        self.currentChanged.connect(self.render_current)
        self.dirty = set(self.figures)

    def set_scene(self, config, report, arrays):
        # Disabling the active near tab emits currentChanged synchronously.
        # Do not render the new far-only data until all tab states are updated.
        self.dirty.clear()
        self.data = (config, report, arrays)
        self.validation.setRowCount(0)
        self.camera_metrics.clearContents()
        self.camera_view.clear()
        self.continuous_view.clear()
        self.show_diagnostics(scene_diagnostics(arrays['far'], report['native'], config))
        for index in (1, 2, 3):
            self.setTabEnabled(index, bool(arrays['near']))
        self.dirty = set(self.figures)
        for row, channel in enumerate(report['collection']['channels']):
            values = [channel[key] for key in ('channel', 'na', 'tilt_deg', 'orders',
                'reference_power_w', 'sample_power_w', 'difference_power_w', 'delta_power_w')]
            for col, value in enumerate(values):
                self.table.setItem(row, col, QTableWidgetItem('{:.7g}'.format(value)))
        totals = report['collection']['totals']
        union = ['Union', '', '', int(arrays['far']['union_mask'].sum())]+[
            totals[key+'_union_power_w'] for key in ('reference', 'sample', 'difference', 'delta')]
        for col, value in enumerate(union):
            self.table.setItem(3, col, QTableWidgetItem(value if isinstance(value, str) else '{:.7g}'.format(value)))
        self.render_current()

    def show_camera(self, report, arrays):
        self.camera_view.show_camera(report, arrays)
        for i, row in enumerate(report['channels']):
            snr = row.get('gray_snr')
            periods = row['ideal_relay_periods_in_roi']
            values = [str(i+1), 'N/A' if snr is None else '{:.6g}'.format(snr), str(row.get('gray_frames', 0)),
                row.get('gray_status', 'Not acquired'), str(row['expected_overfull_pixels']),
                '{:.4g} / {:.4g}'.format(*periods)]
            for j, value in enumerate(values):
                self.camera_metrics.setItem(i, j, QTableWidgetItem(value))
        self.setCurrentIndex(10)

    def show_diagnostics(self, report):
        specs = [('Reflectance R', 'reflectance'), ('Flux into substrate T', 'flux_into_substrate'),
            ('1-R-T (not independent absorption)', 'finite_stack_absorption_by_balance'),
            ('Passive flux check', 'passive_flux_check_passed'),
            ('Modal E/native flux max abs error', 'modal_e_vs_native_flux_max_absolute_error'),
            ('Transverse E weighted relative L2', 'transverse_relative_l2'),
            ('X-mirror E weighted relative L2', 'mirror_relative_l2')]
        self.diagnostics.setRowCount(len(specs)+1)
        for i, (title, key) in enumerate(specs):
            self.diagnostics.setItem(i, 0, QTableWidgetItem(title))
            for j, label in enumerate(('reference', 'sample', 'delta'), 1):
                value = report['fields'][label].get(key)
                text = 'N/A' if value is None else str(value) if isinstance(value, bool) else '{:.7g}'.format(value)
                self.diagnostics.setItem(i, j, QTableWidgetItem(text))
        i = len(specs)
        self.diagnostics.setItem(i, 0, QTableWidgetItem('CH1/CH2 power asymmetry: 2|P1-P2|/(P1+P2)'))
        for j, label in enumerate(('reference', 'sample', 'delta'), 1):
            value = report['side_channel_power_asymmetry'][label]
            self.diagnostics.setItem(i, j, QTableWidgetItem('N/A' if value is None else '{:.7g}'.format(value)))

    def show_validation(self, report):
        cases = report['cases']
        self.validation.setRowCount(len(cases))
        for i, (name, row) in enumerate(cases.items()):
            comparison = row.get('comparison', {})
            changes = [c.get('weighted_complex_field_change') for c in comparison.get('channels', [])]
            values = [name, row['status']]+[
                '{:.5g}'.format(100*changes[j]) if j < len(changes) and changes[j] is not None else 'N/A'
                for j in range(3)]+[str(comparison.get('local_far_step_passed', 'unchecked'))]
            for j, value in enumerate(values):
                self.validation.setItem(i, j, QTableWidgetItem(value))
        self.setCurrentWidget(self.validation)

    def set_logarithmic(self, value):
        self.logarithmic = value
        self.continuous_view.logarithmic = value
        self.continuous_view.render()
        self.dirty = set(self.figures)
        self.render_current()

    def render_current(self, *_):
        if self.data is None or self.currentIndex() >= len(self.figures):
            return
        key = tuple(self.figures)[self.currentIndex()]
        if key not in self.dirty:
            return
        config, report, arrays = self.data
        figure = self.figures[key]
        if key == 'geometry':
            geometry_figure(figure, config)
        elif key in ('xy', 'xz', 'yz'):
            if arrays['near']:
                near_figure(figure, arrays['near'], key, self.logarithmic)
            else:
                figure.clear()
        elif key in ('far', 'far3d'):
            far_figure(figure, arrays['far'], config, self.logarithmic, key == 'far3d')
        else:
            image_figure(figure, arrays['images'], report['collection'], self.logarithmic)
        for ax in figure.axes:
            ax.tick_params(labelsize=8)
            ax.title.set_fontsize(10)
        self.canvases[key].draw_idle()
        self.dirty.remove(key)

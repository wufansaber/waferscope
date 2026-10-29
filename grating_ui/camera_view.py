"""Raw camera figure, numeric brightness table and exact pixel readout."""
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QAbstractItemView, QComboBox, QHeaderView, QLabel, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget)

from grating_scene.camera_display import brightness_rows, count_text, pixel_text
from grating_scene.camera_plot import camera_figure


class CameraView(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        self.canvas = FigureCanvasQTAgg(Figure(figsize=(12, 10)))
        self.canvas.setMinimumSize(1020, 900)
        self.readout = QLabel('Pixel: N/A')
        self.readout.setTextInteractionFlags(Qt.TextSelectableByMouse)
        toolbar = NavigationToolbar2QT(self.canvas, self)
        self.dn_range = QComboBox()
        self.dn_range.addItem('Observed DN (shared)', False)
        self.dn_range.addItem('Full ADC range', True)
        self.dn_range.setToolTip('One common raw DN color range for all Sample and Reference images')
        toolbar.addWidget(self.dn_range)
        layout.addWidget(toolbar)
        layout.addWidget(self.readout)
        layout.addWidget(self.canvas)
        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(('CH', 'Frame', 'Min (DN)', 'Max (DN)',
            'Mean (DN)', 'Sum (DN)', 'At ADC max', 'Full well'))
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.frames = {}
        self.report = self.arrays = None
        self.dn_range.currentIndexChanged.connect(self.redraw)
        self.canvas.mpl_connect('motion_notify_event', self.inspect_pixel)
        self.canvas.mpl_connect('button_press_event', self.inspect_pixel)

    def clear(self):
        self.frames.clear()
        self.report = self.arrays = None
        self.readout.setText('Pixel: N/A')
        self.table.setRowCount(0)
        self.canvas.figure.clear()
        self.canvas.draw_idle()

    def show_camera(self, report, arrays):
        self.report, self.arrays = report, arrays
        self.redraw()
        self.readout.setText('Pixel: N/A')
        rows = brightness_rows(arrays, report)
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            values = [str(row['channel']), row['title']]+[
                '{:.10g}'.format(row[key]) for key in ('min_dn', 'max_dn', 'mean_dn', 'sum_dn')]+[
                count_text(row[key], row['pixels']) for key in ('at_adc_max_pixels', 'full_well_pixels')]
            for j, value in enumerate(values):
                self.table.setItem(i, j, QTableWidgetItem(value))

    def redraw(self, *_):
        if self.arrays is not None:
            self.frames = camera_figure(self.canvas.figure, self.arrays, self.report,
                full_adc=self.dn_range.currentData())
            self.canvas.draw_idle()

    def inspect_pixel(self, event):
        frame = self.frames.get(event.inaxes)
        if frame is None:
            return
        row, data = frame
        value = pixel_text(data, event.xdata, event.ydata)
        if value:
            self.readout.setText('CH{} {} | {}'.format(row['channel'], row['title'], value))

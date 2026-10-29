"""Continuous-angle display and explicit frozen-current input selection."""
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox, QFileDialog,
    QFormLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QSpinBox, QVBoxLayout, QWidget)

from rcwa_prototype.angular_distribution import plot_distribution


class ContinuousView(QWidget):
    def __init__(self):
        super().__init__()
        self.arrays = None
        self.caption = QLabel('Load a completed continuous radiation folder or use Continuous angles to calculate one.')
        self.caption.setWordWrap(True)
        self.surface = QCheckBox('3D direction surface')
        self.surface.toggled.connect(self.render)
        self.logarithmic = True
        self.figure = Figure(figsize=(11, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumSize(1000, 480)
        layout = QVBoxLayout(self)
        for widget in (self.caption, self.surface, NavigationToolbar2QT(self.canvas, self), self.canvas):
            layout.addWidget(widget)

    def show_distribution(self, report, arrays):
        self.arrays = arrays
        s = report['scene']
        self.caption.setText('Continuous particle-current DSC | integral = {:.7g} um^2 | '
            'wavelength {:.4g} nm, sphere {:.4g} nm, gap {:.4g} nm | '
            'periodic source interactions retained; angular convergence and isolation unchecked'.format(
                report['hemisphere_cross_section_um2'], 1000*s['wavelength_um'],
                1000*s['diameter_um'], 1000*s['gap_um']))
        self.render()

    def render(self, *_):
        if self.arrays is not None:
            plot_distribution(self.figure, self.arrays, self.logarithmic, self.surface.isChecked())
            self.canvas.draw_idle()

    def clear(self):
        self.arrays = None
        self.caption.setText('No continuous radiation loaded for this scene.')
        self.figure.clear()
        self.canvas.draw_idle()


class ContinuousInputs(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Continuous angles from exported particle current')
        self.setMinimumWidth(650)
        layout = QVBoxLayout(self)
        caption = QLabel('Select matching grcwa particle-current and bare-reference folders (report.json + modes.npz). '
            'S4 scenes contain discrete modes only; they do not contain the required internal current. '
            'This calculation radiates one cell current; periodic-copy interactions remain.')
        caption.setWordWrap(True)
        layout.addWidget(caption)
        form = QFormLayout()
        layout.addLayout(form)
        self.source, self.reference = QLineEdit(), QLineEdit()
        for title, edit in (('Particle current folder', self.source), ('Bare reference folder', self.reference)):
            row = QWidget()
            line = QHBoxLayout(row)
            line.setContentsMargins(0, 0, 0, 0)
            line.addWidget(edit)
            button = QPushButton('Browse')
            button.clicked.connect(lambda checked=False, e=edit: self.browse(e))
            line.addWidget(button)
            form.addRow(title, row)
        self.polar, self.azimuthal, self.halfwidth = QSpinBox(), QSpinBox(), QSpinBox()
        for title, control, lower, upper, value in (
            ('Polar samples (0 < theta < 90 deg)', self.polar, 4, 256, 24),
            ('Azimuth samples (0 to 360 deg)', self.azimuthal, 8, 1024, 72),
            ('Reciprocal background halfwidth', self.halfwidth, 1, 128, 40)):
            control.setRange(lower, upper)
            control.setValue(value)
            form.addRow(title, control)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def browse(self, edit):
        folder = QFileDialog.getExistingDirectory(self, 'Select current/reference export', edit.text())
        if folder:
            edit.setText(folder)

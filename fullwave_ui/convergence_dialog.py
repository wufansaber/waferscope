"""Editable scan levels; empty rows omit an axis from the current study."""
from PyQt5.QtWidgets import (QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
                            QLineEdit, QMessageBox, QScrollArea, QSpinBox, QVBoxLayout,
                            QWidget, QFrame, QApplication)

from fullwave_prototype.haze_convergence import validate_plan


class ConvergenceDialog(QDialog):
    def __init__(self, plan, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Haze convergence')
        self.setMinimumWidth(520)
        self.plan = dict(plan)
        self.levels = {}
        layout = QVBoxLayout(self)
        body = QWidget()
        form = QFormLayout(body)
        self.samples = QSpinBox()
        self.samples.setRange(1, 128)
        self.samples.setValue(plan['samples'])
        form.addRow('Surfaces per case', self.samples)
        self.tolerance = QDoubleSpinBox()
        self.tolerance.setRange(.1, 50)
        self.tolerance.setValue(100*plan['tolerance'])
        self.tolerance.setSuffix(' %')
        form.addRow('Numerical tolerance', self.tolerance)
        self.confidence_target = QDoubleSpinBox()
        self.confidence_target.setRange(.1, 50)
        self.confidence_target.setValue(100*plan['target_relative_ci'])
        self.confidence_target.setSuffix(' %')
        form.addRow('95% interval half-width target', self.confidence_target)
        for name, levels in plan['scans'].items():
            field = QLineEdit(', '.join('{:g}'.format(level) for level in levels))
            field.setToolTip('Comma-separated levels in comparison order; empty omits this axis')
            self.levels[name] = field
            form.addRow(name.replace('_', ' '), field)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidget(body)
        layout.addWidget(scroll)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        screen = QApplication.primaryScreen()
        height = min(760, int(screen.availableGeometry().height()*.85)) if screen else 700
        self.resize(600, height)

    def accept(self):
        try:
            self.plan.update(samples=self.samples.value(), tolerance=self.tolerance.value()/100,
                target_relative_ci=self.confidence_target.value()/100,
                scans={name: [float(value.strip()) for value in field.text().split(',')]
                       for name, field in self.levels.items() if field.text().strip()})
            validate_plan(self.plan)
        except (ValueError, OSError, KeyError) as exc:
            QMessageBox.warning(self, 'Invalid convergence plan', str(exc))
            return
        super().accept()

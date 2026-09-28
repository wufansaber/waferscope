"""FEM inputs and conversion to the public pipeline command line."""
import os

from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QLabel,
                            QSpinBox, QWidget)

from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.imaging_options import IMAGING_FLAGS, imaging_settings, DETECTOR_OPTIONS


class SimulationControls(QWidget):
    def __init__(self):
        super().__init__()
        self.fields = {}
        self.form = QFormLayout(self)
        self.form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        for name, label, value, low, high, decimals in (
            ('wavelength-nm', 'Wavelength (nm)', 355, 100, 2000, 1),
            ('diameter-nm', 'Particle diameter (nm)', 100, 1, 1000, 1),
            ('particle-n', 'Particle n', 1.59, .01, 20, 4),
            ('particle-k', 'Particle k', 0, 0, 20, 4),
            ('substrate-n', 'Wafer n', 5.57, .01, 20, 4),
            ('substrate-k', 'Wafer k', .387, 0, 20, 4),
            ('incidence-deg', 'Incidence / normal (deg)', 75, 0, 89, 1),
            ('gap-nm', 'Sphere / wafer gap (nm)', 10, 0, 500, 1),
            ('rms-nm', 'Roughness RMS (nm)', 2, 0, 30, 2),
            ('domain-radius-um', 'Domain radius (um)', .22, .12, 2, 3),
            ('magnification', 'Magnification', 20, 1, 200, 1),
            ('objective-focal-mm', 'Objective focal length (mm)', 20, 1, 500, 1),
            ('irradiance-w-um2', 'Irradiance (W/um2)', 1e-6, 1e-12, 1, 12),
            ('defocus-um', 'Image defocus (um)', 0, -1000, 1000, 2),
            ('transmission', 'Power transmission', 1, 0, 1, 3),
        ):
            self.number(name, label, value, low, high, decimals)
        self.choice('polarization', 'Incident polarization', ('p', 's'))
        self.choice('analyzer', 'Analyzer', ('none', 'x', 'y', '45', '135'))
        seed = QSpinBox()
        seed.setRange(0, 1000000)
        seed.setValue(7)
        self.fields['roughness-seed'] = seed
        self.form.addRow('Roughness seed', seed)
        for i, tilt in enumerate((45, -45, 0), 1):
            self.number('na'+str(i), 'Channel {} NA'.format(i), .5, 0, .99, 3)
            self.number('tilt'+str(i), 'Channel {} tilt (deg)'.format(i), tilt, -89.9, 89.9, 1)
        self.choice('accuracy', 'FEM profile', ('balanced', 'preview', 'refined'))
        threads = QSpinBox()
        threads.setRange(1, max(1, os.cpu_count() or 1))
        threads.setValue(min(2, threads.maximum()))
        threads.setKeyboardTracking(False)
        threads.setToolTip('NGSolve worker threads. More threads may require more memory. '
                           'BLAS/OpenMP stay at one to avoid nested parallelism.')
        self.fields['threads'] = threads
        self.form.addRow('FEM threads', threads)
        self.choice('direct-solver', 'Direct solver', ('sparsecholesky', 'umfpack'))
        self.choice('pml-order', 'PML order (0=same)', ('0', '2'))
        self.number('pml-inner-fraction', 'PML full-order fraction', .5, 0, .95, 2)
        self.number('local-mesh-scale', 'Local mesh scale', 1, .4, 1, 2)
        self.choice('pupil-samples', 'Pupil samples', ('128', '64', '32', '96'))
        self.number('beam-waist-um', 'Beam waist (um; 0=plane)', 0, 0, 10000, 3)
        self.choice('beam-order', 'Beam quadrature order', ('5', '3', '7', '9'))
        for key in ('astigmatism', 'coma', 'spherical'):
            self.number(key+'-nm', key.capitalize()+' RMS (nm)', 0, -1000, 1000, 2)
        self.choice('archive-mode', 'Image archive', ('compact', 'full'))
        detector = QCheckBox()
        self.fields['detector'] = detector
        self.form.addRow('Camera response', detector)
        for key, value in DETECTOR_OPTIONS.items():
            if isinstance(value, int):
                control = QSpinBox()
                control.setRange(0 if key == 'seed' else 1, 1000000 if key == 'seed' else 512)
                control.setValue(value)
                self.fields['detector-'+key] = control
                self.form.addRow('Camera '+key.replace('-', ' '), control)
            else:
                self.number('detector-'+key, 'Camera '+key.replace('-', ' '), value,
                            0 if key in ('read-noise-e', 'dark-current-e-s', 'background-e-s', 'quantum-efficiency') else .001,
                            1 if key == 'quantum-efficiency' else 1000000, 3)
        detector.toggled.connect(self.set_camera_enabled)
        self.set_camera_enabled(False)
        self.form.addRow(QLabel('Finite rough patch / homogeneous wafer'))

    def set_camera_enabled(self, enabled):
        for name, field in self.fields.items():
            if name.startswith('detector-'):
                field.setEnabled(enabled)

    def number(self, name, label, value, low, high, decimals):
        field = QDoubleSpinBox()
        field.setDecimals(decimals)
        field.setRange(low, high)
        field.setValue(value)
        field.setKeyboardTracking(False)
        field.setMinimumWidth(125)
        self.fields[name] = field
        self.form.addRow(label, field)

    def choice(self, name, label, values):
        field = QComboBox()
        field.addItems(values)
        self.fields[name] = field
        self.form.addRow(label, field)

    def arguments(self, imaging_only=False):
        result = []
        for name, field in self.fields.items():
            if name.startswith(('na', 'tilt')):
                continue
            if imaging_only and name not in IMAGING_FLAGS:
                continue
            if isinstance(field, QCheckBox):
                if field.isChecked():
                    result.append('--'+name)
                continue
            value = field.currentText() if isinstance(field, QComboBox) else str(field.value())
            result.extend(('--'+name, value))
        for flag, prefix in (('--na', 'na'), ('--tilts', 'tilt')):
            result.append(flag)
            result.extend(str(self.fields[prefix+str(i)].value()) for i in (1, 2, 3))
        args = build_parser().parse_args(result)
        if not imaging_only:
            scattering_settings(args)
        imaging_settings(args)
        return result

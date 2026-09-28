"""FEM inputs and conversion to the public pipeline command line."""
from PyQt5.QtWidgets import (QComboBox, QDoubleSpinBox, QFormLayout, QLabel,
                            QSpinBox, QWidget)

from fullwave_prototype.imaging import ImagingSettings
from fullwave_prototype.imaging_run import build_parser, scattering_settings


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
        self.choice('direct-solver', 'Direct solver', ('sparsecholesky', 'umfpack'))
        self.choice('pml-order', 'PML order (0=same)', ('0', '2'))
        self.number('pml-inner-fraction', 'PML full-order fraction', .5, 0, .95, 2)
        self.number('local-mesh-scale', 'Local mesh scale', 1, .4, 1, 2)
        self.choice('pupil-samples', 'Pupil samples', ('128', '64', '32', '96'))
        self.form.addRow(QLabel('Finite rough patch / homogeneous wafer'))

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

    def arguments(self):
        result = []
        for name, field in self.fields.items():
            if name.startswith(('na', 'tilt')):
                continue
            value = field.currentText() if isinstance(field, QComboBox) else str(field.value())
            result.extend(('--'+name, value))
        for flag, prefix in (('--na', 'na'), ('--tilts', 'tilt')):
            result.append(flag)
            result.extend(str(self.fields[prefix+str(i)].value()) for i in (1, 2, 3))
        args = build_parser().parse_args(result)
        scattering_settings(args)
        ImagingSettings(numerical_apertures=tuple(args.na), tilt_degrees=tuple(args.tilts),
            magnification=args.magnification, focal_length_mm=args.objective_focal_mm,
            irradiance_w_um2=args.irradiance_w_um2, image_defocus_um=args.defocus_um,
            transmission=args.transmission, analyzer=args.analyzer,
            pupil_samples=args.pupil_samples).validate()
        return result

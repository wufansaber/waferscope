"""FEM inputs and conversion to the public pipeline command line."""
import os

from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QLabel,
                            QPushButton, QSpinBox, QWidget)

from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.imaging_options import IMAGING_FLAGS, imaging_settings, DETECTOR_OPTIONS


class SimulationControls(QWidget):
    def __init__(self):
        super().__init__()
        self.fields = {}
        self.form = QFormLayout(self)
        self.form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        self.preset = QComboBox()
        self.preset.addItems(('10 nm PSL preview', '100 nm standard', 'Custom'))
        self.preset_button = QPushButton('Apply preset')
        self.form.addRow('Parameter set', self.preset)
        self.form.addRow(self.preset_button)
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
            ('correlation-length-um', 'Rough correlation L (um)', .05, .001, 1, 4),
            ('rough-width-um', 'Rough patch width (um)', .2, .05, 2, 3),
            ('domain-radius-um', 'Domain radius (um)', .22, .12, 2, 3),
            ('magnification', 'Magnification', 20, 1, 200, 1),
            ('objective-focal-mm', 'Objective focal length (mm)', 20, 1, 500, 1),
            ('illumination-power-w', 'Laser power (W)', 1, 1e-15, 1000, 15),
            ('illumination-length-um', 'Wafer spot length (um)', 2, .01, 10000, 3),
            ('illumination-width-um', 'Wafer spot width (um)', 2, .01, 10000, 3),
            ('defocus-um', 'Image defocus (um)', 0, -1000, 1000, 2),
            ('transmission', 'Power transmission', 1, 0, 1, 3),
        ):
            self.number(name, label, value, low, high, decimals)
        self.fields['illumination-power-w'].setToolTip('Total incident power over the specified wafer-surface footprint.')
        for key in ('illumination-length-um', 'illumination-width-um'):
            self.fields[key].setToolTip('Effective rectangular footprint on the wafer. The FEM uses a '
                                        'locally uniform plane wave, not a finite square top-hat beam.')
        self.irradiance_label = QLabel()
        self.form.addRow('Beam-normal irradiance', self.irradiance_label)
        for key in ('illumination-power-w', 'illumination-length-um',
                    'illumination-width-um', 'incidence-deg'):
            self.fields[key].valueChanged.connect(self.update_irradiance_label)
        self.update_irradiance_label()
        self.choice('polarization', 'Incident polarization', ('p', 's'))
        self.choice('roughness-model', 'Rough surface model', ('legacy', 'gaussian'))
        self.fields['roughness-model'].setToolTip('Gaussian uses stationary RMS and lateral correlation length before '
                                                  'the patch edge taper. Selecting it also sets first-order FEM preview.')
        self.fields['rms-nm'].setToolTip('For Gaussian roughness this is the stationary surface RMS before taper; '
                                       'the FEM patch RMS is reported separately.')
        self.fields['correlation-length-um'].setToolTip('L in C(r)=RMS^2 exp[-(r/L)^2] for Gaussian roughness.')
        self.choice('analyzer', 'Analyzer', ('none', 'x', 'y', '45', '135'))
        rough_grid = QSpinBox()
        rough_grid.setRange(4, 17)
        rough_grid.setValue(5)
        rough_grid.setToolTip('For Gaussian correlation, at least two grid intervals per L are required. '
                              'More nodes can substantially increase FEM memory and runtime.')
        self.fields['rough-grid-points'] = rough_grid
        self.form.addRow('Rough surface nodes / axis', rough_grid)
        seed = QSpinBox()
        seed.setRange(0, 1000000)
        seed.setValue(7)
        self.fields['roughness-seed'] = seed
        self.form.addRow('Roughness seed', seed)
        for i, tilt in enumerate((45, -45, 0), 1):
            self.number('na'+str(i), 'Channel {} NA (0=off)'.format(i), .5, 0, .9999, 4)
            self.number('tilt'+str(i), 'Channel {} tilt (deg)'.format(i), tilt, -89.9, 89.9, 1)
        self.choice('accuracy', 'FEM profile', ('balanced', 'preview', 'refined'))
        self.choice('order', 'FEM order (0=profile)', ('0', '1', '2', '3'))
        self.fields['order'].setToolTip('Gaussian mode starts with first-order FEM to fit this computer. '
                                        'Higher order needs a separate memory/convergence check.')
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
        self.number('mesh-scale', 'Global mesh scale', 1, .4, 4, 2)
        self.number('local-mesh-scale', 'Local mesh scale', 1, .4, 1, 2)
        self.number('particle-mesh-um', 'Particle mesh (um; 0=auto)', 0, 0, .1, 5)
        self.fields['particle-mesh-um'].setSingleStep(.0005)
        self.choice('quadrature-order', 'Source quadrature order', ('auto', '5', '7'))
        self.choice('pupil-samples', 'Pupil samples', ('128', '64', '32', '96'))
        self.choice('padding', 'Image padding', ('8', '4', '2'))
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
        self.fields['roughness-model'].setCurrentText('gaussian')
        self.adjust_rough_grid('gaussian')
        standard = {name: (field.isChecked() if isinstance(field, QCheckBox) else
                           field.currentText() if isinstance(field, QComboBox) else field.value())
                    for name, field in self.fields.items()}
        small_particle = dict(standard, **{'diameter-nm': 10, 'accuracy': 'preview',
            'threads': 2, 'pml-inner-fraction': 0, 'mesh-scale': 1.5,
            'particle-mesh-um': .0015, 'quadrature-order': '7'})
        self._presets = (small_particle, standard)
        self._last_preset_index = 0
        self._applying_preset = False
        self.preset.currentIndexChanged.connect(self.apply_preset)
        self.preset_button.clicked.connect(self.restore_preset)
        self.fields['roughness-model'].currentTextChanged.connect(self.adjust_rough_grid)
        self.apply_preset()
        for field in self.fields.values():
            if isinstance(field, QCheckBox):
                field.toggled.connect(self.mark_custom)
            elif isinstance(field, QComboBox):
                field.currentIndexChanged.connect(self.mark_custom)
            else:
                field.valueChanged.connect(self.mark_custom)

    def apply_preset(self, *_):
        index = self.preset.currentIndex()
        if index >= len(self._presets):
            return
        self._last_preset_index = index
        self._applying_preset = True
        try:
            for name, value in self._presets[index].items():
                field = self.fields[name]
                if isinstance(field, QCheckBox):
                    field.setChecked(value)
                elif isinstance(field, QComboBox):
                    field.setCurrentText(value)
                else:
                    field.setValue(value)
        finally:
            self._applying_preset = False

    def restore_preset(self):
        if self.preset.currentIndex() >= len(self._presets):
            self.preset.setCurrentIndex(self._last_preset_index)
        else:
            self.apply_preset()

    def mark_custom(self, *_):
        if not self._applying_preset and self.preset.currentIndex() != len(self._presets):
            self.preset.setCurrentIndex(len(self._presets))

    def set_camera_enabled(self, enabled):
        for name, field in self.fields.items():
            if name.startswith('detector-'):
                field.setEnabled(enabled)

    def adjust_rough_grid(self, model):
        if model == 'gaussian':
            if self.fields['rough-grid-points'].value() < 9:
                self.fields['rough-grid-points'].setValue(9)
            self.fields['accuracy'].setCurrentText('preview')
            self.fields['order'].setCurrentText('1')
            self.fields['direct-solver'].setCurrentText('umfpack')
            self.fields['pml-inner-fraction'].setValue(0)

    def update_irradiance_label(self, *_):
        from math import cos, radians
        power = self.fields['illumination-power-w'].value()
        length = self.fields['illumination-length-um'].value()
        width = self.fields['illumination-width-um'].value()
        incidence = self.fields['incidence-deg'].value()
        self.irradiance_label.setText('{:.6g} W/um2'.format(
            power/(length*width*cos(radians(incidence)))))

    def load_imaging_settings(self, saved, incidence_deg=None):
        from math import cos, radians
        incidence = (incidence_deg if incidence_deg is not None else
                     self.fields['incidence-deg'].value())
        length = saved.get('illumination_length_um') or 2.
        width = saved.get('illumination_width_um') or 2.
        power = saved.get('illumination_power_w')
        if power is None:
            power = saved['irradiance_w_um2']*length*width*cos(radians(incidence))
        values = {
            'incidence-deg': incidence,
            'magnification': saved['magnification'],
            'objective-focal-mm': saved['focal_length_mm'],
            'illumination-power-w': power,
            'illumination-length-um': length,
            'illumination-width-um': width,
            'defocus-um': saved['image_defocus_um'],
            'transmission': saved['transmission'],
            'analyzer': saved['analyzer'],
            'pupil-samples': str(saved['pupil_samples']),
            'padding': str(saved['padding']),
        }
        for index in (1, 2, 3):
            values['na'+str(index)] = saved['numerical_apertures'][index-1]
            values['tilt'+str(index)] = saved['tilt_degrees'][index-1]
        for name, value in zip(('astigmatism-nm', 'coma-nm', 'spherical-nm'), saved['aberrations_nm']):
            values[name] = value
        detector = saved['detector']
        values['detector'] = detector['enabled']
        values.update({'detector-'+name.replace('_', '-'): value
                       for name, value in detector.items() if name != 'enabled'})
        self._applying_preset = True
        try:
            for name, value in values.items():
                field = self.fields[name]
                if isinstance(field, QCheckBox):
                    field.setChecked(value)
                elif isinstance(field, QComboBox):
                    field.setCurrentText(value)
                else:
                    field.setValue(value)
            self.preset.setCurrentIndex(len(self._presets))
        finally:
            self._applying_preset = False

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

    def arguments(self, imaging_only=False, incidence_deg=None):
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
            if name == 'order' and value == '0':
                continue
            if name == 'quadrature-order' and value == 'auto':
                continue
            result.extend(('--'+name, value))
        for flag, prefix in (('--na', 'na'), ('--tilts', 'tilt')):
            result.append(flag)
            result.extend(str(self.fields[prefix+str(i)].value()) for i in (1, 2, 3))
        args = build_parser().parse_args(result)
        if not imaging_only:
            scattering_settings(args)
        imaging_settings(args, incidence_deg=(incidence_deg if incidence_deg is not None else
                                             self.fields['incidence-deg'].value()))
        return result

"""FEM inputs and conversion to the public pipeline command line."""
import os
from PyQt5.QtCore import QSignalBlocker

from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
                            QFormLayout, QHBoxLayout, QLabel, QLineEdit, QMessageBox,
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
        self.preset.addItems(('10 nm PSL preview', '100 nm standard', '266 nm low-haze camera',
                             '266 nm line X preview', '266 nm line Y preview', '266 nm cross preview', 'Custom'))
        self.preset_button = QPushButton('Apply preset')
        self.form.addRow('Parameter set', self.preset)
        self.form.addRow(self.preset_button)
        for name, label, value, low, high, decimals in (
            ('wavelength-nm', 'Wavelength (nm)', 355, 100, 2000, 1),
            ('diameter-nm', 'Sphere diameter (nm)', 100, 1, 1000, 1),
            ('particle-n', 'Defect n', 1.59, .01, 20, 4),
            ('particle-k', 'Defect k', 0, 0, 20, 4),
            ('substrate-n', 'Wafer n', 5.57, .01, 20, 4),
            ('substrate-k', 'Wafer k', .387, 0, 20, 4),
            ('incidence-deg', 'Incidence / normal (deg)', 75, 0, 89, 1),
            ('gap-nm', 'Defect bottom / wafer gap (nm)', 10, 0, 500, 1),
            ('rms-nm', 'Roughness RMS (nm)', 2, 0, 30, 2),
            ('correlation-length-um', 'PSD length B / L (um)', .05, .001, 1, 4),
            ('psd-exponent', 'ABC exponent C', 3, 2.01, 10, 3),
            ('rough-width-um', 'Rough patch width (um)', .2, .05, 2, 3),
            ('domain-radius-um', 'Domain radius (um)', .22, .12, 10, 3),
            ('magnification', 'Magnification', 20, 1, 200, 1),
            ('objective-focal-mm', 'Objective focal length (mm)', 20, 1, 500, 1),
            ('illumination-power-w', 'Laser power (W)', 1, 1e-15, 1000, 15),
            ('illumination-length-um', 'Wafer spot length (um)', 2, .01, 10000, 3),
            ('illumination-width-um', 'Wafer spot width (um)', 2, .01, 10000, 3),
            ('defocus-um', 'Image defocus (um)', 0, -1000, 1000, 2),
            ('transmission', 'Power transmission', 1, 0, 1, 8),
        ):
            self.number(name, label, value, low, high, decimals)
        from .defect_controls import DefectControls
        self.defect_controls = DefectControls(self)
        self.choice('beam-profile', 'Illumination', ('gaussian', 'plane'))
        self.choice('beam-definition', 'Finite beam definition', ('wafer_gaussian', 'beam_cross_section',
                    'wafer_rectangle', 'measured_wafer'))
        self.choice('beam-phase-assumption', 'Missing envelope phase', ('flat', 'reject'))
        self.beam_measurement_path = QLineEdit()
        beam_browse = QPushButton('...')
        beam_browse.setToolTip('Select wafer intensity/phase NPZ')
        beam_browse.clicked.connect(self.choose_beam_measurement)
        beam_row = QHBoxLayout()
        beam_row.addWidget(self.beam_measurement_path, 1)
        beam_row.addWidget(beam_browse)
        self.form.addRow('Measured beam file', beam_row)
        self.beam_measurement_browse = beam_browse
        self.fields['illumination-power-w'].setToolTip('Gaussian: total beam power. Plane: power over effective rectangular footprint.')
        self.fields['illumination-power-w'].setMinimumWidth(170)
        for key in ('illumination-length-um', 'illumination-width-um'):
            self.fields[key].setToolTip('Effective rectangular footprint on the wafer. The FEM uses a '
                                        'locally uniform plane wave, not a finite square top-hat beam.')
        self.irradiance_label = QLabel()
        self.form.addRow('Reference irradiance', self.irradiance_label)
        for key in ('illumination-power-w', 'illumination-length-um',
                    'illumination-width-um', 'incidence-deg'):
            self.fields[key].valueChanged.connect(self.update_irradiance_label)
        self.update_irradiance_label()
        self.choice('polarization', 'Incident polarization', ('p', 's'))
        self.choice('roughness-model', 'Rough surface model', ('legacy', 'gaussian', 'abc',
                    'measured_height', 'measured_psd'))
        self.measurement_path = QLineEdit()
        self.measurement_browse = QPushButton('...')
        self.measurement_browse.setToolTip('Select measured height or PSD file')
        self.measurement_browse.clicked.connect(self.choose_measurement)
        measurement_row = QHBoxLayout()
        measurement_row.addWidget(self.measurement_path, 1)
        measurement_row.addWidget(self.measurement_browse)
        self.form.addRow('Measurement file', measurement_row)
        self.choice('detrend', 'Height detrend', ('plane', 'mean', 'none'))
        self.choice('missing-policy', 'Missing heights', ('reject', 'nearest'))
        self.fields['roughness-model'].setToolTip('Gaussian and ABC use stationary RMS before the patch edge taper. '
                                                  'ABC has a power-law PSD tail and needs a finer height grid.')
        self.fields['rms-nm'].setToolTip('For Gaussian roughness this is the stationary surface RMS before taper; '
                                       'the FEM patch RMS is reported separately.')
        self.fields['correlation-length-um'].setToolTip('Gaussian: covariance length L. '
                                                       'ABC: NIST PSD roll-off length B, not the same L.')
        self.fields['psd-exponent'].setToolTip('ABC high-frequency exponent C > 2; unused for Gaussian.')
        self.choice('analyzer', 'Analyzer', ('none', 'x', 'y', '45', '135'))
        rough_grid = QSpinBox()
        rough_grid.setRange(4, 41)
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
        self.choice('beam-order', 'Beam quadrature order', ('31', '41', '61', '81', '101', '151', '21', '15', '9', '7', '5', '3'))
        self.fields['beam-profile'].currentTextChanged.connect(self.update_beam_profile)
        self.fields['beam-definition'].currentTextChanged.connect(self.update_beam_profile)
        self.beam_measurement_path.textChanged.connect(self.update_irradiance_label)
        self.fields['beam-order'].currentTextChanged.connect(self.update_irradiance_label)
        self.fields['polarization'].currentTextChanged.connect(self.update_irradiance_label)
        self.fields['beam-phase-assumption'].currentTextChanged.connect(self.update_irradiance_label)
        self.fields['wavelength-nm'].valueChanged.connect(self.update_irradiance_label)
        self.update_beam_profile()
        for key in ('astigmatism', 'coma', 'spherical'):
            self.number(key+'-nm', key.capitalize()+' RMS (nm)', 0, -1000, 1000, 2)
        self.choice('archive-mode', 'Image archive', ('compact', 'full'))
        detector = QCheckBox()
        self.fields['detector'] = detector
        self.form.addRow('Camera response', detector)
        for key, value in DETECTOR_OPTIONS.items():
            if isinstance(value, int):
                control = QSpinBox()
                limits = ((8, 16) if key == 'bit-depth' else (1, 128) if key == 'reference-frames'
                          else (0, 1000000) if key.endswith('seed') else (1, 512))
                control.setRange(*limits)
                control.setValue(value)
                self.fields['detector-'+key] = control
                self.form.addRow('Camera '+key.replace('-', ' '), control)
            else:
                self.number('detector-'+key, 'Camera '+key.replace('-', ' '), value,
                            0 if key in ('read-noise-e', 'dark-current-e-s', 'background-e-s', 'quantum-efficiency',
                                        'offset-dn', 'prnu-fraction', 'dsnu-e', 'dark-current-std-e-s',
                                        'haze-mismatch-fraction') else .001,
                            .5 if key == 'prnu-fraction' else 10 if key == 'haze-mismatch-fraction'
                            else 1 if key in ('quantum-efficiency', 'fill-factor') else 1000000,
                            6 if key in ('quantum-efficiency', 'prnu-fraction', 'fill-factor', 'haze-mismatch-fraction') else 4)
        detector.toggled.connect(self.set_camera_enabled)
        self.set_camera_enabled(False)
        for key, tip in {
            'read-noise-e': 'EMVA 1288 temporal dark noise in electrons, excluding ADC quantization; '
                           'use the minimum-exposure value. Exposure-dependent dark shot noise is added separately.',
            'gain-e-dn': 'Electron/DN conversion gain; EMVA system gain K in DN/e- is its reciprocal.',
            'prnu-fraction': 'EMVA relative RMS response variation; 0.01 means 1%. Fixed between frames.',
            'dsnu-e': 'EMVA DSNU1288: fixed dark-offset spatial RMS in electrons, not e-/s.',
            'dark-current-std-e-s': 'Optional DCNU rate variation, not part of EMVA 1288 Linear.',
            'fill-factor': 'Central-square sensitive area. Keep at 1 if effective QE already includes fill factor.',
            'reference-frames': 'Independent noise frames of the same no-defect rough surface. '
                                  'Detection includes reference haze photon/read noise divided by this count.',
            'haze-mismatch-fraction': 'Residual haze RMS / mean after reference subtraction. '
                                      '0 means a stable, identical surface. SNR uncertainty only; '
                                      'a user assumption, not PSD ensemble statistics.'
        }.items():
            self.fields['detector-'+key].setToolTip(tip)
        self.form.addRow(QLabel('Finite rough patch / homogeneous wafer'))
        self.fields['roughness-model'].setCurrentText('gaussian')
        self.adjust_rough_grid('gaussian')
        from .field_coverage import FieldCoverageControls
        self.field_coverage = FieldCoverageControls(self)
        from .camera_controls import CameraControls
        self.camera_controls = CameraControls(self)
        self.set_camera_enabled(False)
        standard = {name: (field.isChecked() if isinstance(field, QCheckBox) else
                           field.currentText() if isinstance(field, QComboBox) else field.value())
                    for name, field in self.fields.items()}
        small_particle = dict(standard, **{'diameter-nm': 10, 'accuracy': 'preview',
            'threads': 2, 'pml-inner-fraction': 0, 'mesh-scale': 1.5,
            'particle-mesh-um': .0015, 'quadrature-order': '7'})
        from camera_profiles import profile_values
        low_haze = dict(small_particle, **{'diameter-nm': 100, 'wavelength-nm': 266,
            'rms-nm': .2, 'correlation-length-um': .05, 'rough-grid-points': 9,
            'particle-mesh-um': .0225, 'order': '1', 'direct-solver': 'umfpack',
            'pml-order': '0', 'transmission': .0001, 'illumination-fit': 'camera',
            'camera-profile': 'gemi_area', 'gray-snr-pairs': 5000})
        low_haze.update({'detector' if name == 'enabled' else 'detector-'+name.replace('_', '-'): value
                         for name, value in profile_values('gemi_area').items()})
        pattern = dict(low_haze, **{'defect-length-nm': 200, 'defect-width-nm': 50,
            'defect-height-nm': 50, 'defect-rotation-deg': 0, 'particle-mesh-um': .015,
            'rms-nm': 0, 'illumination-fit': 'manual', 'beam-profile': 'plane',
            'illumination-length-um': 2, 'illumination-width-um': 2, 'transmission': .000001,
            'pupil-samples': '32', 'padding': '4', 'gray-snr-pairs': 100})
        self._presets = (small_particle, standard, low_haze)+tuple(
            dict(pattern, **{'defect-shape': shape}) for shape in ('line_x', 'line_y', 'cross'))
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
        self.beam_measurement_path.textChanged.connect(self.mark_custom)

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
        self.field_coverage.update()

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
            if name.startswith(('detector-', 'gray-snr-')):
                field.setEnabled(enabled)

    def adjust_rough_grid(self, model):
        self.fields['psd-exponent'].setEnabled(model == 'abc')
        measured = model.startswith('measured_')
        self.measurement_path.setEnabled(measured)
        self.measurement_browse.setEnabled(measured)
        self.fields['rms-nm'].setEnabled(not measured)
        self.fields['correlation-length-um'].setEnabled(model in ('gaussian', 'abc'))
        self.fields['detrend'].setEnabled(model == 'measured_height')
        self.fields['missing-policy'].setEnabled(model == 'measured_height')
        if model in ('gaussian', 'abc', 'measured_height', 'measured_psd'):
            if self.fields['rough-grid-points'].value() < 9:
                self.fields['rough-grid-points'].setValue(9)
            self.fields['accuracy'].setCurrentText('preview')
            self.fields['order'].setCurrentText('1')
            self.fields['direct-solver'].setCurrentText('umfpack')
            self.fields['pml-inner-fraction'].setValue(0)

    def choose_measurement(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Experimental height / PSD',
            self.measurement_path.text(), 'Measurement (*.npz *.csv)')
        if not path:
            return
        try:
            from fullwave_prototype.experimental_surface import load_height, load_psd
            model = self.fields['roughness-model'].currentText()
            if model == 'measured_height':
                _, report = load_height(path, self.fields['detrend'].currentText(),
                                       self.fields['missing-policy'].currentText())
                self.fields['rough-width-um'].setValue(min(self.fields['rough-width-um'].value(),
                                                          min(report['footprint_xy_um'])))
            else:
                _, _, _, report = load_psd(path)
            self.measurement_path.setText(path)
            self.fields['rms-nm'].setMaximum(max(30., report['rms_nm']))
            self.fields['rms-nm'].setValue(report['rms_nm'])
            self.measurement_path.setToolTip('RMS {:.5g} nm | SHA256 {}'.format(
                report['rms_nm'], report['sha256']))
            self.mark_custom()
        except (ValueError, OSError, KeyError) as exc:
            QMessageBox.warning(self, 'Invalid measurement', str(exc))

    def update_irradiance_label(self, *_):
        from math import cos, radians
        power = self.fields['illumination-power-w'].value()
        length = self.fields['illumination-length-um'].value()
        width = self.fields['illumination-width-um'].value()
        incidence = self.fields['incidence-deg'].value()
        if self.fields['beam-profile'].currentText() == 'gaussian':
            from fullwave_prototype.parameters import Parameters
            from fullwave_prototype.finite_illumination import power_area
            try:
                p = Parameters(beam_profile='gaussian', beam_length_um=length, beam_width_um=width,
                    beam_definition=self.fields['beam-definition'].currentText(),
                    incidence_deg=incidence, wavelength_nm=self.fields['wavelength-nm'].value(),
                    polarization=self.fields['polarization'].currentText() if 'polarization' in self.fields else 'p',
                    beam_order=int(self.fields['beam-order'].currentText()) if 'beam-order' in self.fields else 31)
                if p.beam_definition == 'measured_wafer':
                    from dataclasses import replace
                    from fullwave_prototype.beam_inputs import fingerprint, measured_envelope
                    path = self.beam_measurement_path.text()
                    p = replace(p, beam_measurement_path=path, beam_measurement_sha256=fingerprint(path),
                                beam_phase_assumption=self.fields['beam-phase-assumption'].currentText())
                    _, _, _, info = measured_envelope(path, p.beam_measurement_sha256, p.beam_phase_assumption)
                    for key, extent in (('illumination-length-um', info['input_extent_y_um']),
                                        ('illumination-width-um', info['input_extent_x_um'])):
                        field = self.fields[key]
                        with QSignalBlocker(field):
                            field.setValue(extent)
                irradiance = power/power_area(p)
            except (ValueError, OSError, KeyError):
                self.irradiance_label.setText('Invalid finite beam')
                return
        else:
            irradiance = power/(length*width*cos(radians(incidence)))
        self.irradiance_label.setText('{:.6g} W/um2'.format(irradiance))

    def update_beam_profile(self, *_):
        finite = self.fields['beam-profile'].currentText() == 'gaussian'
        definition = self.fields['beam-definition'].currentText()
        measured = finite and definition == 'measured_wafer'
        self.fields['beam-definition'].setEnabled(finite)
        self.fields['beam-phase-assumption'].setEnabled(measured)
        self.beam_measurement_path.setEnabled(measured)
        self.beam_measurement_browse.setEnabled(measured)
        self.fields['beam-waist-um'].setEnabled(not finite)
        if finite:
            self.fields['beam-waist-um'].setValue(0)
            if definition in ('wafer_rectangle', 'measured_wafer') and int(self.fields['beam-order'].currentText()) < 61:
                self.fields['beam-order'].setCurrentText('61')
        for key, axis, name in (('illumination-length-um', 'Y', 'length'),
                                ('illumination-width-um', 'X', 'width')):
            self.form.labelForField(self.fields[key]).setText(
                (('Gaussian 1/e^2 ' if definition == 'wafer_gaussian' else
                  'Beam cross-section ' if definition == 'beam_cross_section' else
                  'Wafer rectangle ' if definition == 'wafer_rectangle' else 'Measured extent ')+name+' '+axis+' (um)')
                if finite else 'Wafer spot '+name+' (um)')
            self.fields[key].setEnabled(not measured)
            self.fields[key].setToolTip(({
                'wafer_gaussian': 'Target intensity 1/e^2 diameter on the wafer; propagating reconstruction may deform it.',
                'beam_cross_section': 'Beam-normal Gaussian 1/e^2 diameter; Y projects by 1/cos(incidence).',
                'wafer_rectangle': 'Coherent uniform rectangle on the wafer before propagating spectral cutoff.',
                'measured_wafer': 'Pixel support extent from the measured file, not a fitted Gaussian diameter.'
            }[definition]) if finite else 'Effective rectangular footprint for plane-wave power scaling only.')
        self.update_irradiance_label()

        if hasattr(self, 'field_coverage'):
            self.field_coverage.update()

    def choose_beam_measurement(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Measured coherent wafer beam', '', 'Beam NPZ (*.npz)')
        if path:
            self.beam_measurement_path.setText(path)
            self.mark_custom()

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
            'beam-profile': saved.get('illumination_model', 'plane'),
            'beam-definition': saved.get('illumination_definition', 'wafer_gaussian'),
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
            'illumination-fit': saved.get('illumination_fit', 'manual'),
            'minimum-field-illumination': saved.get('minimum_field_illumination', .5),
            'camera-profile': saved.get('camera_profile', 'custom'),
            'gray-snr-pairs': saved.get('gray_snr_pairs', 0),
            'gray-snr-seed': saved.get('gray_snr_seed', 20261016),
        }
        for index in (1, 2, 3):
            values['na'+str(index)] = saved['numerical_apertures'][index-1]
            values['tilt'+str(index)] = saved['tilt_degrees'][index-1]
        for name, value in zip(('astigmatism-nm', 'coma-nm', 'spherical-nm'), saved['aberrations_nm']):
            values[name] = value
        detector = dict({key.replace('-', '_'): value for key, value in DETECTOR_OPTIONS.items()},
                        **saved['detector'])
        values['detector'] = detector['enabled']
        values.update({'detector-'+name.replace('_', '-'): value
                       for name, value in detector.items() if name != 'enabled'})
        self._applying_preset = True
        self._loading_field_coverage = True
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
            self._loading_field_coverage = False
        self.field_coverage.update(resize=False)

    def load_beam_parameters(self, saved):
        self.beam_measurement_path.setText(saved.get('beam_measurement_path', ''))
        for key, source in (('wavelength-nm', 'wavelength_nm'), ('polarization', 'polarization'),
                            ('beam-order', 'beam_order'), ('beam-phase-assumption', 'beam_phase_assumption')):
            if source not in saved:
                continue
            field = self.fields[key]
            if isinstance(field, QComboBox):
                field.setCurrentText(str(saved[source]))
            else:
                field.setValue(saved[source])

    def load_physical_parameters(self, saved, roughness=None):
        self.defect_controls.load(saved)
        for name in ('diameter-nm', 'particle-n', 'particle-k', 'substrate-n', 'substrate-k',
                     'gap-nm', 'domain-radius-um', 'particle-mesh-um', 'threads'):
            if name.replace('-', '_') in saved:
                self.fields[name].setValue(saved[name.replace('-', '_')])
        if roughness and roughness.get('model') in ('gaussian', 'abc'):
            self.fields['roughness-model'].setCurrentText(roughness['model'])
            for name, key in (('rms-nm', 'rms_nm'), ('correlation-length-um', 'correlation_length_um'),
                ('psd-exponent', 'psd_exponent'), ('rough-width-um', 'width_um'),
                ('rough-grid-points', 'grid_points'), ('roughness-seed', 'seed')):
                if key in roughness:
                    self.fields[name].setValue(roughness[key])
        for name in ('order', 'direct-solver', 'pml-order'):
            if name.replace('-', '_') in saved:
                self.fields[name].setCurrentText(str(saved[name.replace('-', '_')]))
        if 'pml_inner_fraction' in saved:
            self.fields['pml-inner-fraction'].setValue(saved['pml_inner_fraction'])

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

    def arguments(self, imaging_only=False, incidence_deg=None, particle_enabled=True, parameters=None):
        if not imaging_only:
            self.defect_controls.update()
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
        if not imaging_only and self.fields['roughness-model'].currentText().startswith('measured_'):
            result.extend(('--measurement-path', self.measurement_path.text()))
        if not imaging_only and self.fields['beam-profile'].currentText() == 'gaussian' and self.fields['beam-definition'].currentText() == 'measured_wafer':
            result.extend(('--beam-measurement-path', self.beam_measurement_path.text()))
        args = build_parser().parse_args(result)
        if imaging_only and parameters is not None and args.illumination_fit == 'camera':
            for flag, value in (('--illumination-length-um', parameters.beam_length_um),
                                ('--illumination-width-um', parameters.beam_width_um)):
                result[result.index(flag)+1] = str(value)
            args = build_parser().parse_args(result)
        if not imaging_only:
            parameters, roughness = scattering_settings(args, particle_enabled=particle_enabled)
            if args.illumination_fit == 'camera':
                for flag, value in (('--illumination-length-um', args.illumination_length_um),
                                    ('--illumination-width-um', args.illumination_width_um)):
                    result[result.index(flag)+1] = str(value)
            if roughness.model.startswith('measured_'):
                self.fields['rms-nm'].setMaximum(max(30., roughness.rms_nm))
                self.fields['rms-nm'].setValue(roughness.rms_nm)
        settings = imaging_settings(args, incidence_deg=(incidence_deg if incidence_deg is not None else
                                             self.fields['incidence-deg'].value()), parameters=parameters)
        self.irradiance_label.setText('{:.6g} W/um2'.format(settings.irradiance_w_um2))
        return result

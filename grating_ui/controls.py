"""Compact, scrollable scene controls with physical units in every dimension."""
from dataclasses import replace

from PyQt5.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout,
    QGridLayout, QGroupBox, QLabel, QSpinBox, QVBoxLayout, QWidget)

from grating_scene.config import SceneConfig
from rcwa_prototype.geometry import Scene


def number(value, low, high, suffix='', integer=False, decimals=4):
    widget = QSpinBox() if integer else QDoubleSpinBox()
    if not integer:
        widget.setDecimals(decimals)
    widget.setRange(low, high)
    widget.setValue(value)
    widget.setSuffix(suffix)
    widget.setMinimumWidth(105)
    return widget


class SceneControls(QWidget):
    def __init__(self):
        super().__init__()
        self.scene_fields, self.fields = {}, {}
        self.base_scene = Scene()
        self.layout = QVBoxLayout(self)
        form = self.section('Structure and illumination')
        self.kind = QComboBox()
        self.kind.addItem('Periodic line grating', 'grating')
        self.kind.addItem('Flat wafer (control)', 'flat')
        form.addRow('Wafer', self.kind)
        specs = [('wavelength_um', 'Wavelength', 266., 100., 2000., ' nm', 1000.),
            ('incidence_deg', 'Incidence from normal', 75., 0., 89.9, ' deg', 1.),
            ('pitch_um', 'Grating pitch', .5, .01, 50., ' um', 1.),
            ('ridge_width_um', 'Ridge width', .25, .001, 50., ' um', 1.),
            ('depth_um', 'Etch depth', .1, .001, 10., ' um', 1.),
            ('lx_um', 'X cell length', .5, .02, 50., ' um', 1.),
            ('diameter_um', 'Sphere diameter', 100., 1., 10000., ' nm', 1000.),
            ('gap_um', 'Sphere bottom gap', 0., 0., 10000., ' nm', 1000.),
            ('substrate_n', 'Wafer n', 5.57, .01, 20., '', 1.),
            ('substrate_k', 'Wafer k', .387, 0., 20., '', 1.),
            ('particle_n', 'Sphere n', 1.59, .01, 20., '', 1.),
            ('particle_k', 'Sphere k', 0., 0., 20., '', 1.)]
        for key, title, value, low, high, unit, scale in specs:
            control = number(value, low, high, unit, decimals=7 if key.endswith(('_n', '_k')) else 4)
            self.scene_fields[key] = (control, scale)
            form.addRow(title, control)
        self.periods = number(2, 1, 100, integer=True)
        form.addRow('Y grating periods / cell', self.periods)
        self.particle = QCheckBox()
        self.particle.setChecked(True)
        form.addRow('Include sphere', self.particle)
        self.polarization = QComboBox()
        for label in ('p', 's'):
            self.polarization.addItem(label.upper(), label)
        form.addRow('Polarization', self.polarization)
        form = self.section('Independent collection channels')
        channels, grid = QWidget(), QGridLayout()
        channels.setLayout(grid)
        grid.addWidget(QLabel('Channel'), 0, 0)
        grid.addWidget(QLabel('NA'), 0, 1)
        grid.addWidget(QLabel('Tilt about Y'), 0, 2)
        self.nas, self.tilts = [], []
        for i, tilt in enumerate((45., -45., 0.), 1):
            na, angle = number(.6, 0., .999, decimals=3), number(tilt, -89.9, 89.9, ' deg', decimals=2)
            na.setToolTip('Air NA; zero disables this channel')
            self.nas.append(na)
            self.tilts.append(angle)
            for column, widget in enumerate((QLabel('CH'+str(i)), na, angle)):
                grid.addWidget(widget, i, column)
        form.addRow(channels)
        form = self.section('Ideal wafer-conjugate relay')
        specs = [('magnification', 'Magnification', 20., 1., 200., ' x'),
            ('transmission', 'Power transmission', 1., 0., 1., ''),
            ('defocus_um', 'Image defocus', 0., -1000., 1000., ' um'),
            ('illumination_power_w', 'Reference-area incident power', 1., 1e-9, 100., ' W'),
            ('normalization_length_um', 'Wafer normalization length', 2., .01, 10000., ' um'),
            ('normalization_width_um', 'Wafer normalization width', 2., .01, 10000., ' um')]
        for key, title, value, low, high, unit in specs:
            control = number(value, low, high, unit, decimals=9 if key in ('illumination_power_w', 'transmission') else 4)
            self.fields[key] = control
            form.addRow(title, control)
        self.analyzer = QComboBox()
        for value in ('none', 'x', 'y', '45', '135'):
            self.analyzer.addItem(value, value)
        form.addRow('Linear analyzer', self.analyzer)
        form = self.section('Discretization and observation')
        self.preset = QComboBox()
        self.preset.addItem('Preview', (3, 16, 4, 25, 128))
        self.preset.addItem('Higher cutoff', (3, 40, 8, 33, 256))
        self.preset.addItem('Custom', None)
        form.addRow('Sampling preset', self.preset)
        self.mx, self.my = number(3, 0, 32, integer=True), number(16, 0, 128, integer=True)
        self.slices = number(4, 1, 64, integer=True)
        self.near_samples = number(25, 5, 129, integer=True)
        self.image_samples = number(128, 16, 1024, integer=True)
        for title, control in (('X halfwidth Mx', self.mx), ('Y halfwidth My', self.my),
            ('Sphere slices', self.slices), ('Near samples / axis', self.near_samples),
            ('Image samples / axis', self.image_samples)):
            form.addRow(title, control)
        self.fields['near_xy_z_um'] = number(.02, -10., 10., ' um')
        self.fields['near_air_height_um'] = number(.15, .001, 100., ' um')
        self.near_enabled = QCheckBox()
        self.near_enabled.setChecked(True)
        form.addRow('Sample native near planes', self.near_enabled)
        form.addRow('XY observation z', self.fields['near_xy_z_um'])
        form.addRow('Air observation height', self.fields['near_air_height_um'])
        self.preset.currentIndexChanged.connect(self.apply_preset)
        self.mx.valueChanged.connect(self.ensure_image_sampling)
        self.my.valueChanged.connect(self.ensure_image_sampling)
        self.layout.addStretch()

    def section(self, title):
        group = QGroupBox(title)
        form = QFormLayout(group)
        form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.layout.addWidget(group)
        return form

    def apply_preset(self):
        values = self.preset.currentData()
        if values:
            for control, value in zip((self.mx, self.my, self.slices, self.near_samples, self.image_samples), values):
                control.setValue(value)

    def ensure_image_sampling(self):
        minimum = 4*max(self.mx.value(), self.my.value())+1
        if self.image_samples.value() < minimum:
            self.image_samples.setValue(min(1024, 2**(minimum-1).bit_length()))

    def config(self):
        scene = replace(self.base_scene, **{key: control.value()/scale for key, (control, scale) in self.scene_fields.items()},
            kind=self.kind.currentData(), periods_y=self.periods.value(), slices=self.slices.value(),
            harmonics_x=self.mx.value(), harmonics_y=self.my.value(), grating_method='analytic_li')
        result = SceneConfig(scene=scene, particle=self.particle.isChecked(),
            polarization=self.polarization.currentData(), numerical_apertures=tuple(c.value() for c in self.nas),
            tilts_deg=tuple(c.value() for c in self.tilts), analyzer=self.analyzer.currentData(),
            near_samples=self.near_samples.value(), image_samples=self.image_samples.value(),
            near_enabled=self.near_enabled.isChecked(),
            **{key: control.value() for key, control in self.fields.items()})
        result.validate()
        return result

    def load_config(self, config):
        self.base_scene = config.scene
        self.preset.setCurrentIndex(2)
        self.kind.setCurrentIndex(self.kind.findData(config.scene.kind))
        for key, (control, scale) in self.scene_fields.items():
            control.setValue(getattr(config.scene, key)*scale)
        self.periods.setValue(config.scene.periods_y)
        self.mx.setValue(config.scene.harmonics_x)
        self.my.setValue(config.scene.harmonics_y)
        self.slices.setValue(config.scene.slices)
        self.particle.setChecked(config.particle)
        self.polarization.setCurrentIndex(self.polarization.findData(config.polarization))
        self.analyzer.setCurrentIndex(self.analyzer.findData(config.analyzer))
        for key, control in self.fields.items():
            control.setValue(getattr(config, key))
        for controls, values in ((self.nas, config.numerical_apertures), (self.tilts, config.tilts_deg)):
            for control, value in zip(controls, values):
                control.setValue(value)
        self.near_samples.setValue(config.near_samples)
        self.near_enabled.setChecked(config.near_enabled)
        self.image_samples.setValue(config.image_samples)

"""Reuse the application's detector controls for the periodic relay camera."""
from dataclasses import asdict, replace

from PyQt5.QtWidgets import QFormLayout, QGroupBox, QPushButton, QSpinBox, QVBoxLayout

from camera_profiles import GEMI_NOTE, profile_values
from scatmech_ui.detector_controls import DetectorGroup

PARAMETER_FIELDS = {'pixel_pitch_um': 'pitch', 'pixels_x': 'nx', 'pixels_y': 'ny', 'exposure_ms': 'exposure',
    'quantum_efficiency': 'qe', 'read_noise_e': 'read_noise', 'dark_current_e_s': 'dark',
    'background_e_s': 'background', 'full_well_e': 'full_well', 'seed': 'seed',
    'gain_e_dn': 'gain', 'offset_dn': 'offset', 'bit_depth': 'bits', 'fill_factor': 'fill',
    'prnu_fraction': 'prnu', 'dsnu_e': 'dsnu', 'dark_current_std_e_s': 'dcnu',
    'fixed_pattern_seed': 'fixed_seed'}


class PeriodicCameraControls(QGroupBox):
    def __init__(self):
        super().__init__('Area camera | periodic relay')
        layout = QVBoxLayout(self)
        self.detector = DetectorGroup()
        self.detector.prnu.setDecimals(7)
        layout.addWidget(self.detector)
        form = QFormLayout()
        self.references, self.pairs, self.seed = QSpinBox(), QSpinBox(), QSpinBox()
        self.references.setRange(1, 128)
        self.pairs.setRange(0, 20000)
        self.pairs.setValue(1000)
        self.pairs.setSingleStep(100)
        self.seed.setRange(0, 2**31-1)
        self.seed.setValue(20261027)
        for title, control in (('Reference frames / exposure', self.references),
                               ('Raw SNR pairs (0=off)', self.pairs), ('Acquisition seed', self.seed)):
            form.addRow(title, control)
        layout.addLayout(form)
        apply = QPushButton('Apply Gemi area parameters')
        apply.setToolTip(GEMI_NOTE)
        apply.clicked.connect(self.apply_gemi)
        layout.addWidget(apply)

    def apply_gemi(self):
        self.load_settings(dict(detector=profile_values('gemi_area'), gray_pairs=self.pairs.value(), gray_seed=self.seed.value()))

    def load_settings(self, settings):
        from detector_model import DetectorParameters
        parameters = DetectorParameters(**settings['detector'])
        parameters.validate()
        values = asdict(parameters)
        self.detector.setChecked(values['enabled'])
        for key, name in PARAMETER_FIELDS.items():
            getattr(self.detector, name).setValue(values[key])
        self.references.setValue(values['reference_frames'])
        self.pairs.setValue(settings['gray_pairs'])
        self.seed.setValue(settings['gray_seed'])

    def settings(self):
        p = replace(self.detector.parameters(), reference_frames=self.references.value())
        p.validate()
        return dict(detector=asdict(p), gray_pairs=self.pairs.value(), gray_seed=self.seed.value())

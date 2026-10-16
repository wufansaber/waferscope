"""Camera presets and raw-gray acquisition controls, separated from FEM geometry."""
from PyQt5.QtWidgets import QComboBox, QLabel, QPushButton, QSpinBox

from camera_profiles import GEMI_NOTE, profile_values


class CameraControls:
    def __init__(self, controls):
        self.controls, self.applying = controls, False
        controls.choice('camera-profile', 'Camera parameter source', ('custom', 'gemi_area'))
        self.apply = QPushButton('Apply Gemi area parameters')
        self.apply.setToolTip(GEMI_NOTE)
        self.apply.clicked.connect(self.apply_profile)
        controls.form.addRow(self.apply)
        self.note = QLabel('')
        self.note.setWordWrap(True)
        controls.form.addRow(self.note)
        for key, label, value, maximum in (('gray-snr-pairs', 'Gray frame pairs (0=off)', 0, 20000),
                                         ('gray-snr-seed', 'Gray acquisition seed', 20261016, 2**31-1)):
            spin = QSpinBox()
            spin.setRange(0, maximum)
            spin.setValue(value)
            spin.setKeyboardTracking(False)
            controls.fields[key] = spin
            controls.form.addRow(label, spin)
        controls.fields['gray-snr-pairs'].setSingleStep(100)
        controls.fields['gray-snr-pairs'].setToolTip('0 disables acquisition. Otherwise 100..20000 independent sample/reference pairs; '
            'SNR from raw DN ROI sums. Same static haze; no TDI accumulation or spatial mismatch sampling.')
        controls.fields['camera-profile'].currentTextChanged.connect(self.update)
        for name, field in controls.fields.items():
            if name.startswith('detector-'):
                field.valueChanged.connect(self.mark_custom)
        self.update()

    def apply_profile(self):
        self.applying = True
        try:
            for name, value in profile_values('gemi_area').items():
                field = self.controls.fields['detector' if name == 'enabled' else 'detector-'+name.replace('_', '-')]
                field.setChecked(value) if name == 'enabled' else field.setValue(value)
            self.controls.fields['camera-profile'].setCurrentText('gemi_area')
            self.controls.fields['gray-snr-pairs'].setValue(5000)
        finally:
            self.applying = False
        self.update()

    def mark_custom(self, *_):
        if not self.applying and not getattr(self.controls, '_applying_preset', False):
            self.controls.fields['camera-profile'].setCurrentText('custom')

    def update(self, *_):
        gemi = self.controls.fields['camera-profile'].currentText() == 'gemi_area'
        self.note.setText('Area approximation | no TDI | QE: 266 nm only | gain 16 e-/DN: estimated' if gemi else '')
        self.note.setToolTip(GEMI_NOTE if gemi else '')

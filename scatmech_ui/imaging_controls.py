"""Controls for the independent objective/tube-lens image planes."""

from PyQt5.QtWidgets import QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox, QLabel, QPushButton, QSpinBox

from .imaging import ImagingParameters


class ImagingGroup(QGroupBox):
    def __init__(self, objective_focal_length):
        super().__init__("Independent image planes")
        form = QFormLayout(self)
        form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        self.magnification = self._spin(form, "Magnification", 1, 200, 20, " x")
        self.tube = QLabel()
        form.addRow("Tube focal length", self.tube)
        def update_tube():
            self.tube.setText("{:.3f} mm".format(
                objective_focal_length.value() * self.magnification.value()))
        self.magnification.valueChanged.connect(update_tube)
        objective_focal_length.valueChanged.connect(update_tube)
        update_tube()
        self.irradiance = self._spin(form, "Incident irradiance", 0.000001, 1e9, 0.001,
                                     " mW/um^2", 6)
        self.irradiance.setToolTip("Incident power per area normal to the beam; uniform plane wave.")
        self.defocus = self._spin(form, "Image defocus", -10000, 10000, 0, " um")
        self.transmission = self._spin(form, "Power transmission", 0, 1, 1, "")
        self.analyzer = QComboBox()
        for label, value in (("None", "none"), ("Local X", "x"), ("Local Y", "y"),
                             ("45 deg", "45"), ("135 deg", "135")):
            self.analyzer.addItem(label, value)
        form.addRow("Pupil analyzer", self.analyzer)
        self.samples = QComboBox()
        for n in (32, 64, 96, 128):
            self.samples.addItem(str(n), n)
        self.samples.setCurrentIndex(1)
        form.addRow("Pupil samples / axis", self.samples)
        self.multipoles = QSpinBox()
        self.multipoles.setRange(0, 60)
        self.multipoles.setSpecialValueText("Automatic")
        form.addRow("Multipole order", self.multipoles)
        self.run_button = QPushButton("Calculate three image planes")
        form.addRow(self.run_button)

    @staticmethod
    def _spin(form, label, minimum, maximum, value, suffix, decimals=3):
        box = QDoubleSpinBox()
        box.setDecimals(decimals)
        box.setRange(minimum, maximum)
        box.setValue(value)
        box.setSuffix(suffix)
        form.addRow(label, box)
        return box

    def parameters(self):
        return ImagingParameters(magnification=self.magnification.value(),
                                 irradiance_mw_um2=self.irradiance.value(),
                                 image_defocus_um=self.defocus.value(),
                                 transmission=self.transmission.value(),
                                 analyzer=self.analyzer.currentData(),
                                 pupil_samples=self.samples.currentData(),
                                 multipole_order=self.multipoles.value())

"""Camera parameters shared by three independently sampled detectors."""

from PyQt5.QtWidgets import QFormLayout, QGroupBox, QSpinBox

from .detector import DetectorParameters
from .imaging_controls import ImagingGroup


class DetectorGroup(QGroupBox):
    def __init__(self):
        super().__init__("Detectors (three independent ROIs)")
        self.setCheckable(True)
        self.setChecked(True)
        form = QFormLayout(self)
        form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        spin = ImagingGroup._spin
        self.pitch = spin(form, "Pixel pitch", 0.01, 100, 6.5, " um")
        self.nx = self._integer(form, "ROI columns", 1, 512, 32)
        self.ny = self._integer(form, "ROI rows", 1, 512, 32)
        self.exposure = spin(form, "Exposure", 0.000001, 1000000, 1, " ms", 6)
        self.qe = spin(form, "Quantum efficiency", 0, 1, 0.6, "")
        self.qe.setToolTip("Set from the detector data at the laser wavelength; default is illustrative.")
        self.read_noise = spin(form, "Read noise (RMS)", 0, 10000, 2, " e-")
        self.dark = spin(form, "Dark current / pixel", 0, 1e9, 0.1, " e-/s")
        self.background = spin(form, "Background / pixel", 0, 1e12, 0, " e-/s")
        self.background.setToolTip("Uniform independent electron rate, not a wafer roughness or stray-light solution.")
        self.full_well = spin(form, "Full well / pixel", 1, 1e9, 30000, " e-", 0)
        self.seed = self._integer(form, "Noise seed", 0, 2147483647, 0)

    @staticmethod
    def _integer(form, label, minimum, maximum, value):
        box = QSpinBox()
        box.setRange(minimum, maximum)
        box.setValue(value)
        form.addRow(label, box)
        return box

    def parameters(self):
        return DetectorParameters(enabled=self.isChecked(), pixel_pitch_um=self.pitch.value(),
                                  pixels_x=self.nx.value(), pixels_y=self.ny.value(),
                                  exposure_ms=self.exposure.value(), quantum_efficiency=self.qe.value(),
                                  read_noise_e=self.read_noise.value(), dark_current_e_s=self.dark.value(),
                                  background_e_s=self.background.value(), full_well_e=self.full_well.value(),
                                  seed=self.seed.value())

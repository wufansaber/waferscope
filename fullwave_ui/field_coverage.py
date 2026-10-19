"""Camera-linked illumination inputs, separate from the scattering controls."""
from PyQt5.QtCore import QSignalBlocker
from PyQt5.QtWidgets import QLabel

from detector_model import DetectorParameters
from fullwave_prototype.field_coverage import camera_footprints, nominal_dimensions


class FieldCoverageControls:
    def __init__(self, controls):
        self.controls = controls
        controls.choice('illumination-fit', 'Illumination sizing', ('camera', 'manual'))
        controls.number('minimum-field-illumination', 'Min. field / center illumination', .5, .05, .95, 3)
        controls.fields['minimum-field-illumination'].setSingleStep(.05)
        controls.fields['minimum-field-illumination'].setToolTip(
            'Minimum sampled incident wafer flux / center flux over all active camera footprints. '
            'Does not enforce uniform haze or suppress coherent speckle.')
        controls.fields['illumination-fit'].setToolTip(
            'Camera: fit a Gaussian wafer beam to the independent sensor footprints; '
            'Calculate uses whole-surface SPM + local FEM. Manual preserves the entered beam.')
        self.extent = QLabel()
        self.extent.setWordWrap(True)
        controls.form.addRow('Wafer camera field X x Y', self.extent)
        for name in ('illumination-fit', 'beam-profile', 'beam-definition', 'roughness-model'):
            controls.fields[name].currentTextChanged.connect(self.update)
        for name in ('magnification', 'detector-pixel-pitch-um', 'detector-pixels-x',
                     'detector-pixels-y', 'minimum-field-illumination', 'incidence-deg',
                     'na1', 'na2', 'na3', 'tilt1', 'tilt2', 'tilt3'):
            controls.fields[name].valueChanged.connect(self.update)
        self.update()

    def update(self, *_, resize=True):
        import numpy as np
        f = self.controls.fields
        auto = f['illumination-fit'].currentText() == 'camera'
        supported = (f['beam-profile'].currentText() == 'gaussian' and
                     f['beam-definition'].currentText() in ('wafer_gaussian', 'beam_cross_section') and
                     f['roughness-model'].currentText() in ('gaussian', 'abc'))
        if auto and not supported:
            f['illumination-fit'].setCurrentText('manual')
            return
        f['minimum-field-illumination'].setEnabled(auto)
        for name in ('illumination-length-um', 'illumination-width-um'):
            f[name].setEnabled(not auto and f['beam-definition'].currentText() != 'measured_wafer')
        try:
            detector = DetectorParameters(pixel_pitch_um=f['detector-pixel-pitch-um'].value(),
                pixels_x=f['detector-pixels-x'].value(), pixels_y=f['detector-pixels-y'].value())
            footprints = camera_footprints(detector, f['magnification'].value(),
                [f['na'+str(i)].value() for i in (1, 2, 3)],
                [f['tilt'+str(i)].value() for i in (1, 2, 3)])
            self.extent.setText(' | '.join('Ch{} {:.3f} x {:.3f} um'.format(
                q['channel'], q['width_x_um'], q['height_y_um']) for q in footprints) or 'No active objective')
            if not auto or not resize or getattr(self.controls, '_loading_field_coverage', False):
                return
            if not supported:
                self.extent.setText('Camera sizing requires a Gaussian beam')
                return
            length, width = nominal_dimensions(footprints, f['minimum-field-illumination'].value())
            if f['beam-definition'].currentText() == 'beam_cross_section':
                length *= np.cos(np.deg2rad(f['incidence-deg'].value()))
            for name, value in (('illumination-length-um', length), ('illumination-width-um', width)):
                with QSignalBlocker(f[name]):
                    f[name].setValue(value)
            self.controls.update_irradiance_label()
        except ValueError as exc:
            self.extent.setText(str(exc))

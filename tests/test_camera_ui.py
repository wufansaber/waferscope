"""Both UI input paths preserve the EMVA camera units and new ADC controls."""
import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from PyQt5.QtWidgets import QApplication
from fullwave_ui.controls import SimulationControls
from fullwave_prototype.imaging_run import build_parser
from fullwave_prototype.imaging_options import imaging_settings
from scatmech_ui.detector_controls import DetectorGroup


def test_scatmech_camera_parameters():
    app = QApplication.instance() or QApplication([])
    group = DetectorGroup()
    group.bits.setValue(12)
    group.gain.setValue(3.25)
    group.dsnu.setValue(1.2)
    group.prnu.setValue(.015)
    p = group.parameters()
    p.validate()
    assert p.gain_e_dn == 3.25 and p.bit_depth == 12
    assert p.dsnu_e == 1.2 and p.prnu_fraction == .015
    group.close()
    app.processEvents()


def test_fem_camera_cli_controls():
    app = QApplication.instance() or QApplication([])
    controls = SimulationControls()
    controls.fields['detector'].setChecked(True)
    controls.fields['detector-gain-e-dn'].setValue(3.25)
    controls.fields['detector-dsnu-e'].setValue(1.2)
    controls.fields['detector-bit-depth'].setValue(12)
    controls.fields['detector-reference-frames'].setValue(4)
    controls.fields['detector-haze-mismatch-fraction'].setValue(.12)
    args = build_parser().parse_args(controls.arguments())
    p = imaging_settings(args).detector
    assert p.enabled and p.gain_e_dn == 3.25
    assert p.dsnu_e == 1.2 and p.bit_depth == 12
    assert p.reference_frames == 4 and p.haze_mismatch_fraction == .12
    controls.close()
    app.processEvents()

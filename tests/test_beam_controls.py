"""Measured-beam UI dimensions, phase policy and source restoration."""
import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import numpy as np
import pytest
pytest.importorskip('PyQt5')
from PyQt5.QtWidgets import QApplication
from fullwave_ui.controls import SimulationControls
from fullwave_prototype.imaging_run import build_parser, scattering_settings


def test_measured_beam_controls_export_and_restore(tmp_path):
    app = QApplication.instance() or QApplication([])
    widget = SimulationControls()
    axis = (np.arange(16)-8)*.15
    x, y = np.meshgrid(axis, axis)
    path = tmp_path/'beam.npz'
    np.savez(path, x_um=axis, y_um=axis, intensity=np.exp(-2*(x*x+y*y)))
    widget.fields['beam-definition'].setCurrentText('measured_wafer')
    widget.beam_measurement_path.setText(str(path))
    assert not widget.fields['illumination-length-um'].isEnabled()
    assert widget.fields['beam-phase-assumption'].isEnabled()
    assert widget.fields['illumination-length-um'].value() == pytest.approx(2.4)
    args = build_parser().parse_args(widget.arguments())
    p, _ = scattering_settings(args)
    assert p.beam_definition == 'measured_wafer'
    assert p.beam_length_um == pytest.approx(2.4)
    assert p.beam_width_um == pytest.approx(2.4)
    widget.fields['beam-phase-assumption'].setCurrentText('reject')
    with pytest.raises(ValueError, match='does not determine phase'):
        widget.arguments()
    np.savez(path, x_um=axis, y_um=axis, intensity=np.ones_like(x), phase_envelope_rad=.1*x)
    p, _ = scattering_settings(build_parser().parse_args(widget.arguments()))
    widget.load_beam_parameters(dict(beam_measurement_path=str(path), wavelength_nm=266,
                                    beam_phase_assumption='reject', beam_order=101, polarization='s'))
    assert widget.fields['beam-order'].currentText() == '101'
    assert widget.fields['wavelength-nm'].value() == 266
    widget.fields['beam-definition'].setCurrentText('wafer_rectangle')
    assert widget.fields['illumination-length-um'].isEnabled()
    assert not widget.beam_measurement_path.isEnabled()
    widget.close()

"""UI geometry-domain correction without changing physical or accuracy inputs."""
from dataclasses import asdict
import os

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from PyQt5.QtCore import QSignalBlocker
from PyQt5.QtWidgets import QApplication
from fullwave_prototype.defects import required_radius_um
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.parameters import Parameters
from fullwave_ui.controls import SimulationControls


@pytest.fixture
def controls():
    app = QApplication.instance() or QApplication([])
    widget = SimulationControls()
    widget.preset.setCurrentText('266 nm cross preview')
    yield widget
    widget.close()
    app.processEvents()


def test_domain_grows_without_changing_accuracy_or_defect(controls):
    f = controls.fields
    original = {name: f[name].currentText() if hasattr(f[name], 'currentText') else f[name].value()
        for name in ('mesh-scale', 'local-mesh-scale', 'particle-mesh-um', 'accuracy',
                     'order', 'pml-order', 'threads', 'wavelength-nm', 'na1', 'na2', 'na3')}
    f['defect-length-nm'].setValue(1000)
    f['defect-height-nm'].setValue(200)
    args = build_parser().parse_args(controls.arguments())
    p, _ = scattering_settings(args)
    assert p.domain_radius_um > required_radius_um(p)
    assert (p.defect_length_nm, p.defect_width_nm, p.defect_height_nm) == (1000, 50, 200)
    assert 'Invalid' not in controls.defect_controls.info.text()
    for name, value in original.items():
        current = f[name].currentText() if hasattr(f[name], 'currentText') else f[name].value()
        assert current == value


def test_calculate_preflight_handles_blocked_edit_signals(controls):
    f = controls.fields
    with QSignalBlocker(f['defect-length-nm']):
        f['defect-length-nm'].setValue(1000)
    assert f['domain-radius-um'].value() == .22
    args = build_parser().parse_args(controls.arguments())
    p, _ = scattering_settings(args)
    assert p.domain_radius_um >= 1.2*required_radius_um(p)
    assert p.domain_radius_um < 1.2*required_radius_um(p)+.00101


def test_reported_screenshot_limit_is_rounded_up(controls):
    controls.defect_controls.fit_domain(.539931)
    assert controls.fields['domain-radius-um'].value() == .648


def test_already_sufficient_domain_is_preserved(controls):
    f = controls.fields
    f['domain-radius-um'].setValue(1)
    f['defect-length-nm'].setValue(1000)
    f['defect-length-nm'].setValue(200)
    assert f['domain-radius-um'].value() == 1


def test_small_manual_domain_is_corrected_and_presets_can_reset(controls):
    f = controls.fields
    f['defect-length-nm'].setValue(1000)
    f['domain-radius-um'].setValue(.22)
    assert f['domain-radius-um'].value() > .5
    controls.preset.setCurrentText('266 nm cross preview')
    assert f['domain-radius-um'].value() == .22


def test_saved_physics_is_not_resized_during_load(controls):
    saved = Parameters(defect_shape='cross', defect_length_nm=1000, gap_nm=10, domain_radius_um=.51)
    saved.validate()
    controls.load_physical_parameters(asdict(saved))
    assert controls.fields['domain-radius-um'].value() == .51
    assert controls.fields['defect-length-nm'].value() == 1000


def test_backend_still_rejects_invalid_domain():
    with pytest.raises(ValueError, match='entirely inside'):
        Parameters(defect_shape='cross', defect_length_nm=1000, domain_radius_um=.22).validate()

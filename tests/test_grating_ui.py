from dataclasses import replace
import importlib.util
import os

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
pytestmark = pytest.mark.skipif(importlib.util.find_spec('PyQt5') is None, reason='PyQt5 environment required')


@pytest.fixture
def app():
    from PyQt5.QtWidgets import QApplication
    instance = QApplication.instance() or QApplication([])
    yield instance


def test_controls_roundtrip_and_independent_apertures(app):
    from grating_scene.config import SceneConfig
    from grating_ui.controls import SceneControls
    controls = SceneControls()
    config = replace(SceneConfig(), numerical_apertures=(.3, 0., .96), tilts_deg=(30., -45., 0.), magnification=40.)
    controls.load_config(config)
    assert controls.config() == config
    controls.my.setValue(48)
    assert controls.image_samples.value() >= 193
    controls.close()


def test_new_page_starts_without_native_engine_or_heavy_calculation(app):
    from grating_ui.widget import GratingWidget
    widget = GratingWidget()
    assert not widget.cancel.isEnabled()
    assert not widget.reimage.isEnabled()
    assert widget.results.count() == 14
    assert widget.results.table.rowCount() == 4
    assert not widget.check_button.isEnabled()
    assert not widget.camera_button.isEnabled()
    assert widget.wait_for_shutdown()
    widget.close()


def test_continuous_result_load_toggle_and_self_contained_zip(app, tmp_path, monkeypatch):
    from tests.test_continuous_angles import fixtures, mock_gate
    from rcwa_prototype.angular_distribution import run, load_distribution
    from grating_ui.widget import GratingWidget
    from PyQt5.QtWidgets import QFileDialog
    import zipfile
    source, reference = fixtures(tmp_path)
    mock_gate(monkeypatch)
    output = tmp_path/'continuous'
    run(source, reference, output, 4, 8, 2)
    widget = GratingWidget()
    assert widget.open_result(output)
    assert widget.results.currentWidget() is widget.results.continuous_view
    assert widget.loaded_result is None
    assert widget.export.isEnabled()
    widget.results.continuous_view.surface.setChecked(True)
    widget.results.set_logarithmic(False)
    assert widget.results.continuous_view.figure.axes[0].name == '3d'
    assert not widget.results.continuous_view.logarithmic
    target = tmp_path/'export.zip'
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', lambda *a: (str(target), 'ZIP (*.zip)'))
    widget.export_result()
    with zipfile.ZipFile(target) as bundle:
        bundle.extractall(tmp_path/'relocated')
    report, _ = load_distribution(tmp_path/'relocated')
    assert report['hemisphere_cross_section_um2'] > 0
    assert widget.wait_for_shutdown()
    widget.close()


def test_validation_table_does_not_treat_blocked_case_as_passed(app):
    from grating_ui.results import SceneResults
    results = SceneResults()
    results.show_validation(dict(cases=dict(y_cutoff=dict(status='blocked_memory_preflight'))))
    assert results.validation.item(0, 1).text() == 'blocked_memory_preflight'
    assert results.validation.item(0, 5).text() == 'unchecked'
    results.close()


def test_diagnostics_table_uses_explicit_na_for_zero_difference_field(app):
    from grating_ui.results import SceneResults
    results = SceneResults()
    fields = dict(reference=dict(reflectance=.2, passive_flux_check_passed=True),
        sample=dict(reflectance=.2, passive_flux_check_passed=True), delta={})
    results.show_diagnostics(dict(fields=fields,
        side_channel_power_asymmetry=dict(reference=0., sample=0., delta=None)))
    assert results.diagnostics.item(0, 1).text() == '0.2'
    assert results.diagnostics.item(0, 3).text() == 'N/A'
    assert results.diagnostics.item(3, 2).text() == 'True'
    results.close()


def test_periodic_camera_uses_existing_gemi_area_inputs_without_tdi(app):
    from camera_profiles import gemi_area_parameters
    from grating_ui.camera_controls import PeriodicCameraControls
    controls = PeriodicCameraControls()
    controls.apply_gemi()
    from detector_model import DetectorParameters
    assert DetectorParameters(**controls.settings()['detector']) == gemi_area_parameters()
    assert controls.settings()['gray_pairs'] == 1000
    controls.close()


def test_grating_is_attached_to_the_actual_standalone_entry(app):
    from scatmech_ui.app import ScatmechStandaloneWindow
    window = ScatmechStandaloneWindow()
    assert window.tabs.count() == 3
    assert window.tabs.widget(2) is window.grating_widget
    assert window.tabs.tabText(2) == 'S4 Grating (periodic)'
    window.close()


@pytest.mark.parametrize('index', [0, 1, 2])
def test_contact_grating_examples_apply_without_starting_a_job(app, index):
    from grating_ui.widget import GratingWidget
    from grating_scene.examples import example_folder
    from PyQt5.QtCore import QProcess
    widget = GratingWidget()
    widget.example_choice.setCurrentIndex(index)
    if not example_folder(widget.root, widget.example_choice.currentData()).is_dir():
        pytest.skip('Original demonstration fields are not available')
    widget.apply_example()
    config = widget.controls.config()
    assert config.scene.gap_um == 0 and config.particle and config.scene.kind == 'grating'
    assert widget.process.state() == QProcess.NotRunning
    if index == 2:
        assert config.magnification == 40 and config.transmission == 1e-6
        assert widget.camera_controls.settings()['gray_pairs'] == 5000
    widget.set_busy(True)
    assert not widget.open_example_button.isEnabled()
    widget.set_busy(False)
    widget.close()


def test_loading_far_only_camera_while_a_near_tab_is_active(app):
    from grating_ui.results import SceneResults
    from grating_scene.examples import EXAMPLES, example_folder
    from grating_scene.pipeline import load_scene
    from fullwave_ui.runtime_paths import application_root
    a = example_folder(application_root(), EXAMPLES[1]['key'])
    b = example_folder(application_root(), EXAMPLES[2]['key'])/'optical_scene'
    if not a.is_dir() or not b.is_dir():
        pytest.skip('Original demonstration fields are not available')
    results = SceneResults()
    results.set_scene(*load_scene(a))
    results.setCurrentIndex(2)
    results.set_scene(*load_scene(b))
    assert all(not results.isTabEnabled(i) for i in (1, 2, 3))
    assert results.currentIndex() not in (1, 2, 3)
    results.close()


def test_camera_numeric_table_and_pixel_readout_for_legacy_report(app):
    from types import SimpleNamespace
    import numpy as np
    from grating_ui.camera_view import CameraView
    view = CameraView()
    arrays = {}
    for i in (1, 2, 3):
        arrays['ch{}_digital_dn'.format(i)] = np.full((2, 3), 1023, dtype=np.uint16)
        arrays['ch{}_reference_mean_dn'.format(i)] = np.full((2, 3), 64.5)
        arrays['ch{}_difference_dn'.format(i)] = np.full((2, 3), 958.5)
        arrays['ch{}_saturated'.format(i)] = np.ones((2, 3), dtype=bool)
    view.show_camera(dict(settings=dict(detector=dict(bit_depth=10))), arrays)
    assert view.table.rowCount() == 9
    assert view.table.item(0, 3).text() == '1023'
    assert view.table.item(0, 6).text() == '6/6 (100%)'
    assert view.table.item(1, 7).text() == 'N/A'
    assert next(iter(view.frames)).images[0].get_clim()[0] == 64.5
    view.dn_range.setCurrentIndex(1)
    assert next(iter(view.frames)).images[0].get_clim() == (0, 1023)
    ax = list(view.frames)[2]
    view.inspect_pixel(SimpleNamespace(inaxes=ax, xdata=2., ydata=1.))
    assert view.readout.text() == 'CH1 Sample - reference | x=2  y=1  DN=958.5'
    view.clear()
    assert view.table.rowCount() == 0 and view.readout.text() == 'Pixel: N/A'
    assert not view.frames
    view.table.close()
    view.close()

"""Check examples and actual GUI -> frozen guard -> S4 -> camera execution."""
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace
from unittest.mock import patch
import zipfile


def wait(app, widget, timeout=240):
    from PyQt5.QtCore import QProcess
    start = time.monotonic()
    while widget.process.state() != QProcess.NotRunning:
        app.processEvents()
        if time.monotonic()-start > timeout:
            widget.cancel_job()
            widget.wait_for_shutdown()
            raise TimeoutError('Packaged grating worker exceeded smoke timeout')
        time.sleep(.03)
    app.processEvents()
    if widget.process.exitCode() or widget.process.exitStatus() != QProcess.NormalExit:
        raise RuntimeError(widget.log.toPlainText())


def run(output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    os.environ['FARFIELD_RESULTS_DIR'] = str(output/'results')
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtGui import QFont, QFontDatabase
    from scatmech_ui.app import ScatmechStandaloneWindow
    from grating_scene.camera import load_camera
    from grating_scene.pipeline import digest, load_scene
    app = QApplication.instance() or QApplication([])
    font = Path(os.environ.get('WINDIR', 'C:/Windows'))/'Fonts/segoeui.ttf'
    if font.is_file():
        QFontDatabase.addApplicationFont(str(font))
        app.setFont(QFont('Segoe UI', 9))
    window = ScatmechStandaloneWindow()
    window.resize(1500, 960)
    window.tabs.setCurrentIndex(2)
    widget = window.grating_widget
    window.show()
    try:
        if not (Path(widget.engine.text()).is_file() and Path(widget.executable.text()).is_file()):
            raise RuntimeError('Portable worker or S4 executable missing')
        widget.example_choice.setCurrentIndex(1)
        widget.open_example()
        if widget.loaded_result is None:
            raise RuntimeError('Packaged original near-field example did not load')
        example = widget.loaded_result
        raw_example = digest(example/'sample/modes.npz')
        widget.results.setCurrentIndex(2)
        app.processEvents()
        window.grab().save(str(output/'grating_near_ui.png'))
        widget.example_choice.setCurrentIndex(2)
        widget.apply_example()
        if widget.controls.config().magnification != 40.:
            raise RuntimeError('Camera example parameters were not applied')
        widget.open_example()
        if widget.camera_result is None:
            raise RuntimeError('Packaged camera example did not load')
        load_camera(widget.camera_result)
        view = widget.results.camera_view
        if view.table.rowCount() != 9:
            raise RuntimeError('Packaged camera brightness table is missing')
        view.dn_range.setCurrentIndex(1)
        if next(iter(view.frames)).images[0].get_clim() != (0, 1023):
            raise RuntimeError('Packaged ADC display range failed')
        view.dn_range.setCurrentIndex(0)
        ax = next(iter(view.frames))
        view.inspect_pixel(SimpleNamespace(inaxes=ax, xdata=0., ydata=0.))
        if 'DN=' not in view.readout.text():
            raise RuntimeError('Packaged exact camera pixel readout failed')
        app.processEvents()
        window.grab().save(str(output/'grating_camera_ui.png'))
        widget.results.setCurrentWidget(view.table)
        app.processEvents()
        window.grab().save(str(output/'grating_camera_brightness_ui.png'))
        widget.example_choice.setCurrentIndex(0)
        widget.apply_example()
        widget.start_job()
        wait(app, widget)
        calculated = widget.loaded_result
        if calculated is None or calculated == example:
            raise RuntimeError('Actual packaged S4 calculation did not load')
        config, native, arrays = load_scene(calculated)
        if config.scene.gap_um != 0 or not config.particle or config.scene.kind != 'grating' or len(arrays['near']) != 21:
            raise RuntimeError('Packaged calculation is not the contact sphere + grating near-field scene')
        if native['native']['sample']['backend'] != 'S4':
            raise RuntimeError('Packaged result is not native S4')
        raw = digest(calculated/'sample/modes.npz')
        widget.controls.fields['magnification'].setValue(40.)
        widget.controls.fields['transmission'].setValue(.000001)
        widget.controls.nas[0].setValue(0.)
        widget.start_job(reuse=True)
        wait(app, widget)
        updated = widget.loaded_result
        if updated == calculated or digest(updated/'sample/modes.npz') != raw:
            raise RuntimeError('Packaged optical replay failed or changed native fields')
        widget.camera_controls.apply_gemi()
        widget.camera_controls.detector.nx.setValue(4)
        widget.camera_controls.detector.ny.setValue(8)
        widget.camera_controls.pairs.setValue(1000)
        widget.start_camera()
        wait(app, widget)
        if widget.camera_result is None:
            raise RuntimeError('Actual packaged area camera failed')
        camera, camera_arrays = load_camera(widget.camera_result)
        if len(camera.get('brightness', [])) != 9:
            raise RuntimeError('Packaged camera did not save raw brightness statistics')
        if camera['channels'][0]['signal_e'] != 0:
            raise RuntimeError('NA=0 was ignored by the packaged camera')
        archive = output/'original_grating_bundle.zip'
        with patch('grating_ui.widget.QFileDialog.getSaveFileName', return_value=(str(archive), 'ZIP (*.zip)')):
            widget.export_result()
        with zipfile.ZipFile(archive) as bundle:
            if bundle.read('optical_scene/sample/modes.npz') != (updated/'sample/modes.npz').read_bytes():
                raise RuntimeError('Packaged export lost original S4 data')
        if digest(example/'sample/modes.npz') != raw_example:
            raise RuntimeError('Loading/replaying modified the bundled original example')
        app.processEvents()
        window.grab().save(str(output/'grating_calculated_camera_ui.png'))
        report = dict(status='passed', examples=3, worker=widget.engine.text(), s4=widget.executable.text(),
            native_scene=str(calculated), replay_scene=str(updated), camera=str(widget.camera_result),
            original_native_sha256=raw, native_backend='S4', particle=True, gap_um=0.,
            near_planes=['XY', 'XZ', 'YZ'], disabled_camera_signal_e=camera['channels'][0]['signal_e'],
            raw_zip=str(archive), gray_pairs=1000, quantitative_hybrid_certified=False)
        report.update(camera_brightness_rows=view.table.rowCount(), exact_pixel_readout=True,
            shared_dn_display_ranges=True)
        (output/'grating_package_smoke.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    finally:
        widget.wait_for_shutdown()
        window.close()

"""Exercise Calculate, Update optics and cancellation in the real Qt process bridge."""
import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import time
from unittest.mock import patch
import zipfile

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def wait(app, widget, timeout=180):
    from PyQt5.QtCore import QProcess
    start = time.monotonic()
    while widget.process.state() != QProcess.NotRunning:
        app.processEvents()
        if time.monotonic()-start > timeout:
            widget.cancel_job()
            widget.wait_for_shutdown()
            raise TimeoutError('Qt worker did not complete in time')
        time.sleep(.02)
    app.processEvents()


def main():
    from PyQt5.QtWidgets import QApplication
    from grating_ui.widget import GratingWidget
    from grating_scene.config import SceneConfig
    from grating_scene.pipeline import digest, load_scene
    from rcwa_prototype.geometry import Scene
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ['FARFIELD_RESULTS_DIR'] = str(args.output.resolve())
    app = QApplication.instance() or QApplication([])
    widget = GratingWidget()
    config = SceneConfig(scene=Scene(harmonics_x=1, harmonics_y=3, slices=2,
        grating_method='analytic_li', gap_um=0.), near_samples=5, image_samples=32)
    widget.controls.load_config(config)
    widget.start_job()
    wait(app, widget)
    if widget.loaded_result is None:
        raise RuntimeError('Calculate failed: '+widget.log.toPlainText())
    first = widget.loaded_result
    load_scene(first)
    raw = digest(first/'sample/modes.npz')
    widget.controls.nas[0].setValue(0.)
    widget.controls.fields['magnification'].setValue(40.)
    widget.controls.fields['transmission'].setValue(.000001)
    widget.start_job(reuse=True)
    wait(app, widget)
    second = widget.loaded_result
    if second == first:
        raise RuntimeError('Update optics failed: '+widget.log.toPlainText())
    _, _, arrays = load_scene(second)
    if digest(second/'sample/modes.npz') != raw or arrays['images']['ch1_sample_w_um2'].any():
        raise RuntimeError('Reused optics changed the native data or ignored the disabled channel')
    archive = args.output/'original_scene_bundle.zip'
    with patch('grating_ui.widget.QFileDialog.getSaveFileName', return_value=(str(archive), 'ZIP (*.zip)')):
        widget.export_result()
    with zipfile.ZipFile(archive) as bundle:
        if hashlib.sha256(bundle.read('sample/modes.npz')).hexdigest() != raw:
            raise RuntimeError('Export did not retain the original native data')
    widget.camera_controls.detector.nx.setValue(4)
    widget.camera_controls.detector.ny.setValue(4)
    widget.camera_controls.detector.pitch.setValue(2.5)
    widget.camera_controls.detector.exposure.setValue(.001)
    widget.camera_controls.pairs.setValue(1000)
    widget.start_camera()
    wait(app, widget)
    from grating_scene.camera import load_camera
    if widget.camera_result is None:
        raise RuntimeError('Camera job failed: '+widget.log.toPlainText())
    camera = widget.camera_result
    load_camera(camera)
    camera_archive = args.output/'camera_raw_bundle.zip'
    with patch('grating_ui.widget.QFileDialog.getSaveFileName', return_value=(str(camera_archive), 'ZIP (*.zip)')):
        widget.export_result()
    with zipfile.ZipFile(camera_archive) as bundle:
        if hashlib.sha256(bundle.read('optical_scene/sample/modes.npz')).hexdigest() != raw:
            raise RuntimeError('Camera export lost original optical data')
    cancellation_config = replace(config, scene=replace(config.scene, harmonics_x=3, harmonics_y=16, slices=4),
                                  near_samples=25, image_samples=128)
    widget.controls.load_config(cancellation_config)
    widget.start_job()
    deadline = time.monotonic()+2.
    while time.monotonic() < deadline:
        app.processEvents()
        time.sleep(.02)
    cancelled = widget.output
    widget.cancel_job()
    wait(app, widget)
    if (cancelled/'scene_report.json').is_file():
        raise RuntimeError('Cancellation unexpectedly published a completed result')
    if widget.loaded_result != second:
        raise RuntimeError('A cancelled calculation replaced the previous result')
    resource = json.loads(cancelled.with_suffix('.resources.json').read_text())
    if resource['status'] != 'stopped_cancelled' or resource['surviving_owned_pids']:
        raise RuntimeError('Native child cancellation did not complete cleanly')
    report = dict(status='passed', calculated=str(first), optics_updated=str(second),
        camera=str(camera), camera_archive=str(camera_archive), cancelled=str(cancelled),
        original_native_sha256=raw, raw_archive=str(archive), cancellation=resource)
    (args.output/'qprocess_smoke.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    widget.close()
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

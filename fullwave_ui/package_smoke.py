"""Exercise the frozen GUI -> guard -> FEM -> result-view path offscreen."""
import json
import os
from pathlib import Path


def run(output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    os.environ['FARFIELD_RESULTS_DIR'] = str(output/'results')
    from PyQt5.QtCore import QTimer, QProcess
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication
    from scatmech_ui.app import ScatmechStandaloneWindow
    app = QApplication.instance() or QApplication([])
    font = Path(os.environ.get('WINDIR', 'C:/Windows'))/'Fonts/segoeui.ttf'
    if font.is_file():
        QFontDatabase.addApplicationFont(str(font))
        app.setFont(QFont('Segoe UI', 9))
    window = ScatmechStandaloneWindow()
    window.tabs.setCurrentIndex(1)
    widget = window.fullwave_widget
    if widget.controls.fields['gap-nm'].value() != 0 or widget.controls.fields['rms-nm'].value() != 0:
        raise RuntimeError('New simulations must default to a contact sphere on a flat wafer')
    # Keep this bounded solver/cache smoke independent of whole-beam defaults.
    widget.controls.fields['beam-profile'].setCurrentText('plane')
    widget.controls.fields['illumination-fit'].setCurrentText('manual')
    widget.controls.fields['accuracy'].setCurrentText('preview')
    widget.controls.fields['pupil-samples'].setCurrentText('32')
    widget.controls.fields['rms-nm'].setValue(0)
    widget.controls.fields['threads'].setValue(4)
    widget.controls.fields['beam-waist-um'].setValue(0)
    widget.controls.fields['beam-order'].setCurrentText('31')
    widget.controls.fields['coma-nm'].setValue(15)
    widget.controls.fields['detector'].setChecked(True)
    widget.controls.arguments()
    widget.timeout.setValue(3)
    widget.process.finished.connect(app.quit)
    widget.process.errorOccurred.connect(app.quit)
    window.show()
    QTimer.singleShot(0, widget.start_job)
    QTimer.singleShot(210000, app.quit)
    app.exec_()
    if not widget.wait_for_shutdown():
        raise RuntimeError('Packaged worker did not stop')
    if widget.output is None:
        raise RuntimeError('Packaged worker was not started')
    resources = json.loads((widget.output/'resources.json').read_text(encoding='utf-8'))
    if resources['status'] != 'completed' or resources['surviving_owned_pids']:
        raise RuntimeError('Packaged calculation failed: '+str(resources))
    if widget.process.exitStatus() != QProcess.NormalExit or widget.process.exitCode() != 0:
        raise RuntimeError('Packaged worker exited abnormally')
    if 'Displayed result:' not in widget.loaded.text():
        raise RuntimeError('Packaged GUI did not load result')
    original = widget.output
    original_scattering = json.loads((original/'scattering/report.json').read_text(encoding='utf-8'))
    background = original_scattering['metrics']['without']
    if not background['zero_rhs_skipped_factorization'] or background['factor_seconds'] != 0:
        raise RuntimeError('Packaged flat reference did not take the exact zero-RHS path')
    widget.controls.fields['coma-nm'].setValue(0)
    QTimer.singleShot(0, lambda: widget.start_job(source=original/'scattering'))
    app.exec_()
    if not widget.wait_for_shutdown():
        raise RuntimeError('Packaged re-imaging did not stop')
    second = json.loads((widget.output/'resources.json').read_text(encoding='utf-8'))
    report = json.loads((widget.output/'imaging_report.json').read_text(encoding='utf-8'))
    physical = report['scattering_parameters']
    if physical['gap_nm'] != 0 or physical['defect_shape'] != 'sphere':
        raise RuntimeError('Packaged FEM lost the contact-sphere geometry')
    if (second['status'] != 'completed' or second['surviving_owned_pids']
            or widget.process.exitCode() != 0 or widget.process.exitStatus() != QProcess.NormalExit
            or not all(c['pupil_cache_hit'] for c in report['channels'])):
        raise RuntimeError('Packaged cache re-imaging failed')
    if not (widget.output/'detector.png').is_file() or report['archive_mode'] != 'full':
        raise RuntimeError('Packaged detector or full archive missing')
    widget.results.setCurrentIndex(2)
    app.processEvents()
    if not window.grab().save(str(output/'packaged_gui.png')):
        raise RuntimeError('Packaged GUI screenshot failed')
    (output/'package_smoke.json').write_text(json.dumps(dict(status='passed',
        engine=widget.interpreter.text(), original_result=str(original), result=str(widget.output),
        pupil_cache_hits=[c['pupil_cache_hit'] for c in report['channels']],
        gap_nm=physical['gap_nm'], center_height_nm=physical['diameter_nm']/2,
        zero_rhs_skipped_factorization=True,
        displayed=widget.loaded.text()), indent=2), encoding='utf-8')
    window.close()

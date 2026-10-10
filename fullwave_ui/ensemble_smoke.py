"""Exercise the UI -> NIST reference -> independent FEM -> result-view path."""
import json
import os
from pathlib import Path


def run(output, model='gaussian', measurement_path=None):
    if model not in ('gaussian', 'abc', 'measured_psd'):
        raise ValueError('Smoke PSD model must be gaussian, abc or measured_psd')
    if model == 'measured_psd' and measurement_path is None:
        raise ValueError('Measured PSD smoke requires a measurement file')
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    os.environ['FARFIELD_RESULTS_DIR'] = str(output/'results')
    from PyQt5.QtCore import QProcess, QTimer
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
    widget.controls.fields['wavelength-nm'].setValue(266)
    widget.controls.fields['roughness-model'].setCurrentText(model)
    if measurement_path is not None:
        widget.controls.measurement_path.setText(str(Path(measurement_path).resolve()))
        widget.controls.fields['direct-solver'].setCurrentText('sparsecholesky')
    widget.controls.fields['rms-nm'].setValue(2)
    widget.controls.fields['correlation-length-um'].setValue(.2 if model == 'abc' else .05)
    widget.controls.fields['psd-exponent'].setValue(4.5 if model == 'abc' else 3.)
    widget.controls.fields['rough-width-um'].setValue(.2)
    widget.controls.fields['rough-grid-points'].setValue(9)
    widget.controls.fields['beam-waist-um'].setValue(0)
    widget.controls.fields['beam-profile'].setCurrentText('plane')
    widget.controls.arguments(particle_enabled=False)
    widget.ensemble_count.setValue(2)
    widget.ensemble_taper.setValue(.025)
    widget.ensemble_area.setValue(.04)
    widget.timeout.setValue(3)
    widget.process.finished.connect(app.quit)
    widget.process.errorOccurred.connect(app.quit)
    window.show()
    QTimer.singleShot(0, widget.start_ensemble)
    QTimer.singleShot(210000, app.quit)
    app.exec_()
    if not widget.wait_for_shutdown():
        raise RuntimeError('Haze ensemble worker did not stop')
    if widget.output is None or widget.process.exitStatus() != QProcess.NormalExit or widget.process.exitCode() != 0:
        raise RuntimeError('Haze ensemble worker failed: '+widget.log.toPlainText()[-2000:])
    resources = json.loads((widget.output/'resources.json').read_text(encoding='utf-8'))
    report = json.loads((widget.output/'haze_ensemble_report.json').read_text(encoding='utf-8'))
    if (resources['status'] != 'completed' or resources['surviving_owned_pids'] or
            report['status'] != 'complete' or report['metrics']['samples'] != 2 or
            not report['nist_reference'] or not report.get('finite_patch_reference') or
            report['roughness']['model'] != model or
            not (widget.output/'haze_ensemble_finite.png').is_file() or
            'FEM / finite patch' not in widget.results.ensemble_metrics.toPlainText() or
            not widget.results.isTabEnabled(widget.results.indexOf(widget.results.ensemble_tabs))):
        raise RuntimeError('Haze ensemble result failed validation')
    if not window.grab().save(str(output/'haze_ensemble_ui.png')):
        raise RuntimeError('Haze ensemble UI screenshot failed')
    (output/'haze_ensemble_smoke.json').write_text(json.dumps(dict(status='passed',
        result=str(widget.output), metrics=report['metrics']), indent=2), encoding='utf-8')
    window.close()


if __name__ == '__main__':
    import sys
    if len(sys.argv) not in (2, 3, 4):
        raise SystemExit('Usage: python -m fullwave_ui.ensemble_smoke OUTPUT [gaussian|abc|measured_psd] [MEASUREMENT]')
    run(sys.argv[1], sys.argv[2] if len(sys.argv) >= 3 else 'gaussian',
        sys.argv[3] if len(sys.argv) == 4 else None)

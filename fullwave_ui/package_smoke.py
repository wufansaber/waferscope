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
    widget.controls.fields['accuracy'].setCurrentText('preview')
    widget.controls.fields['pupil-samples'].setCurrentText('32')
    widget.controls.fields['rms-nm'].setValue(0)
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
    widget.results.setCurrentIndex(2)
    app.processEvents()
    if not window.grab().save(str(output/'packaged_gui.png')):
        raise RuntimeError('Packaged GUI screenshot failed')
    (output/'package_smoke.json').write_text(json.dumps(dict(status='passed',
        engine=widget.interpreter.text(), result=str(widget.output),
        displayed=widget.loaded.text()), indent=2), encoding='utf-8')
    window.close()

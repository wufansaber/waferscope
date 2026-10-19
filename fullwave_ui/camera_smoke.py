"""Source/frozen GUI -> engine -> raw ADC statistics -> result-table validation."""
import json
import os
from pathlib import Path


def run(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    os.environ['FARFIELD_RESULTS_DIR'] = str(output/'results')
    from PyQt5.QtCore import QTimer, QProcess
    from PyQt5.QtGui import QFont, QFontDatabase
    from PyQt5.QtWidgets import QApplication
    import numpy as np
    from gray_statistics import measure_roi_series
    from scatmech_ui.app import ScatmechStandaloneWindow
    from .results import ResultsTabs
    app = QApplication.instance() or QApplication([])
    font = Path(os.environ.get('WINDIR', 'C:/Windows'))/'Fonts/segoeui.ttf'
    if font.is_file():
        QFontDatabase.addApplicationFont(str(font))
        app.setFont(QFont('Segoe UI', 9))
    window = ScatmechStandaloneWindow()
    window.tabs.setCurrentIndex(1)
    widget = window.fullwave_widget
    if not widget.load_hybrid_result(source):
        raise RuntimeError('Cannot load the low-RMS optical source')
    widget.controls.camera_controls.apply.click()
    widget.controls.fields['gray-snr-pairs'].setValue(5000)
    widget.process.finished.connect(app.quit)
    widget.process.errorOccurred.connect(app.quit)
    watchdog = QTimer()
    watchdog.setSingleShot(True)
    watchdog.timeout.connect(app.quit)
    watchdog.start(180000)
    window.show()
    QTimer.singleShot(0, widget.start_camera_job)
    app.exec_()
    watchdog.stop()
    if not widget.wait_for_shutdown() or widget.output is None:
        raise RuntimeError('Camera worker did not complete')
    resources = json.loads((widget.output/'resources.json').read_text(encoding='utf-8'))
    if (resources['status'] != 'completed' or resources['surviving_owned_pids']
            or widget.process.exitCode() or widget.process.exitStatus() != QProcess.NormalExit):
        raise RuntimeError('Camera worker failed: '+str(resources))
    report = json.loads((widget.output/'imaging_report.json').read_text(encoding='utf-8'))
    measured = []
    with np.load(widget.output/'gray_raw_stacks.npz') as saved:
        for i, channel in enumerate(report['channels'], 1):
            raw = saved['ch{}_gray_sample_raw_dn'.format(i)]
            ref = saved['ch{}_gray_reference_raw_dn'.format(i)]
            if raw.dtype != np.uint16 or raw.max() >= 1024 or raw.shape != (5000, 20, 20):
                raise RuntimeError('Incorrect raw ADC shape/range')
            result = measure_roi_series(raw.sum(axis=(1, 2), dtype=float),
                ref.mean(axis=1).sum(axis=(1, 2)), 400)
            if not channel['detector']['gray_valid'] or not np.isclose(
                    result['snr'], channel['detector']['gray_snr'], rtol=1e-13):
                raise RuntimeError('Gray SNR cannot be independently reproduced from raw DN')
            measured.append(result['snr'])
    tabs = ResultsTabs()
    tabs.load(widget.output)
    tabs.setCurrentWidget(tabs.camera_tabs)
    tabs.camera_tabs.setCurrentWidget(tabs.camera_metrics)
    for width in (1200, 640):
        tabs.resize(width, 850)
        tabs.show()
        app.processEvents()
        if not tabs.grab().save(str(output/'camera_metrics_{}.png'.format(width))):
            raise RuntimeError('Camera metric screenshot failed')
    tabs.camera_tabs.setCurrentWidget(tabs.gray_acquisition)
    app.processEvents()
    tabs.grab().save(str(output/'gray_acquisition_ui.png'))
    window.resize(1680, 980)
    app.processEvents()
    window.grab().save(str(output/'camera_gui.png'))
    result = dict(status='passed', optical_source=str(source), result=str(widget.output),
        engine=widget.interpreter.text(), profile=report['settings']['camera_profile'],
        gray_snr_from_reloaded_raw=measured, frames=5000, tdi_accumulation=False)
    (output/'camera_ui_smoke.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    tabs.close()
    window.close()
    app.processEvents()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    run(args.source, args.output)

"""Re-acquire a camera from saved optical images, without another FEM solve."""
import codecs
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path

from PyQt5.QtCore import QProcessEnvironment
from detector_model import DetectorParameters
from gray_acquisition import validate_acquisition
from fullwave_prototype.imaging_options import DETECTOR_OPTIONS
from .runtime_paths import application_root, results_root


def start_camera(widget):
    fields = widget.controls.fields
    if not fields['detector'].isChecked():
        raise ValueError('Enable Camera response before acquiring gray frames')
    source = widget.camera_optical_result
    if source is None:
        raise ValueError('Load a completed optical image result first')
    parameters = DetectorParameters(enabled=True, **{key.replace('-', '_'):
        fields['detector-'+key].value() for key in DETECTOR_OPTIONS})
    parameters.validate()
    pairs, seed = (fields[name].value() for name in ('gray-snr-pairs', 'gray-snr-seed'))
    validate_acquisition(parameters, pairs, seed)
    interpreter = Path(widget.interpreter.text())
    if not interpreter.is_file():
        raise ValueError('Select the packaged FEM engine or native Python environment')
    root = application_root()
    widget.output = results_root(root)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    widget.output.mkdir(parents=True, exist_ok=False)
    config = widget.output/'camera_settings.json'
    config.write_text(json.dumps(dict(detector=asdict(parameters), gray_snr_pairs=pairs,
        gray_snr_seed=seed, camera_profile=fields['camera-profile'].currentText()), indent=2), encoding='utf-8')
    widget.active_report = widget.output/'resources.json'
    widget.cancel_file = widget.output/'cancel.request'
    command = ['-m', 'fullwave_prototype.guarded_run', '--report', str(widget.active_report),
        '--max-memory-gib', str(widget.memory.value()), '--timeout-seconds', str(widget.timeout.value()*60),
        '--cancel-file', str(widget.cancel_file), 'detection_replay', '--source', str(source),
        '--output', str(widget.output), '--camera-settings-json', str(config)]
    widget.job_kind, widget.imaging_only_job = 'camera', True
    widget.process.setWorkingDirectory(str(root))
    env = QProcessEnvironment.systemEnvironment()
    for key, value in (('PYTHONIOENCODING', 'utf-8'), ('PYTHONUNBUFFERED', '1'), ('PYINSTALLER_RESET_ENVIRONMENT', '1')):
        env.insert(key, value)
    widget.process.setProcessEnvironment(env)
    widget.decoder = codecs.getincrementaldecoder('utf-8')('replace')
    widget.log.clear()
    widget.set_running(True)
    widget.status.setText('Camera-only acquisition | saved optics unchanged | '+str(source))
    widget.process.start(str(interpreter), command)

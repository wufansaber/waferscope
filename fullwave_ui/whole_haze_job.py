"""Launch the explicitly approximate whole-illumination haze calculation."""
import codecs
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path

from PyQt5.QtCore import QProcessEnvironment
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.whole_surface import plan_surface
from .runtime_paths import application_root, results_root


def start_whole_haze(widget):
    args = build_parser().parse_args(widget.controls.arguments(particle_enabled=False))
    p, roughness = scattering_settings(args, particle_enabled=False)
    plan_surface(p, roughness, widget.whole_coverage.value(), widget.whole_variance.value(),
                 widget.ensemble_taper.value())
    settings = imaging_settings(args, parameters=p)
    interpreter = Path(widget.interpreter.text())
    if not interpreter.is_file():
        raise ValueError('Select the isolated solver Python environment')
    root = application_root()
    widget.output = results_root(root)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    widget.output.mkdir(parents=True, exist_ok=False)
    inputs = dict(parameters=asdict(p), roughness=asdict(roughness), imaging_settings=asdict(settings),
        samples=widget.ensemble_count.value(), coverage=widget.whole_coverage.value(),
        variance_fraction=widget.whole_variance.value(), taper_width_um=widget.ensemble_taper.value())
    path = widget.output/'whole_haze_input.json'
    path.write_text(json.dumps(inputs, indent=2), encoding='utf-8')
    widget.cancel_file = widget.output/'cancel.request'
    widget.active_report = widget.output/'resources.json'
    widget.job_kind = 'whole_haze'
    command = ['-m', 'fullwave_prototype.guarded_run', '--report', str(widget.active_report),
        '--max-memory-gib', str(widget.memory.value()), '--timeout-seconds', str(widget.timeout.value()*60),
        '--cancel-file', str(widget.cancel_file), 'whole_haze', '--input', str(path), '--output', str(widget.output)]
    (widget.output/'request.json').write_text(json.dumps(dict(interpreter=str(interpreter), arguments=command),
                                            indent=2), encoding='utf-8')
    widget.process.setWorkingDirectory(str(root))
    env = QProcessEnvironment.systemEnvironment()
    for key, value in (('PYTHONIOENCODING', 'utf-8'), ('PYTHONUNBUFFERED', '1'), ('PYINSTALLER_RESET_ENVIRONMENT', '1')):
        env.insert(key, value)
    widget.process.setProcessEnvironment(env)
    widget.decoder = codecs.getincrementaldecoder('utf-8')('replace')
    widget.log.clear()
    widget.set_running(True)
    widget.status.setText('Whole-beam haze | first-order surface approximation | '+str(widget.output))
    widget.process.start(str(interpreter), command)

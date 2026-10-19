"""Resource-guarded new/re-imaged hybrid jobs, separate from local FEM sources."""
import codecs
from dataclasses import asdict, replace
from datetime import datetime
import json
from pathlib import Path

from PyQt5.QtCore import QProcessEnvironment
from fullwave_prototype.hybrid_run import load_sources
from fullwave_prototype.imaging_run import build_parser, scattering_settings
from fullwave_prototype.imaging_options import imaging_settings
from fullwave_prototype.parameters import Parameters
from .runtime_paths import application_root, results_root


def start_hybrid(widget, source=None):
    if source is None:
        args = build_parser().parse_args(widget.controls.arguments())
        p, rough = scattering_settings(args)
        taper = max(widget.ensemble_taper.value(), rough.width_um/(rough.grid_points-1))
        rough = replace(rough, taper_fraction=taper/rough.width_um)
        rough.validate(p)
        settings = imaging_settings(args, parameters=p)
        inputs = dict(parameters=asdict(p), roughness=asdict(rough), imaging_settings=asdict(settings),
            mother_range_factor=int(widget.mother_factor.currentText()),
            coverage=widget.whole_coverage.value(), variance_fraction=widget.whole_variance.value(),
            taper_width_um=widget.ensemble_taper.value(), quadrature_order=args.quadrature_order or (7 if p.order >= 3 else 5))
    else:
        _, metadata = load_sources(source)
        p = Parameters(**metadata['parameters'])
        args = build_parser().parse_args(widget.controls.arguments(imaging_only=True, parameters=p,
                                                                   incidence_deg=p.incidence_deg))
        settings = imaging_settings(args, parameters=p)
    interpreter = Path(widget.interpreter.text())
    if not interpreter.is_file():
        raise ValueError('Select the isolated FEM Python environment')
    root = application_root()
    widget.output = results_root(root)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    widget.output.mkdir(parents=True, exist_ok=False)
    command_tail = []
    if source is None:
        path = widget.output/'hybrid_request_input.json'
        path.write_text(json.dumps(inputs, indent=2), encoding='utf-8')
        command_tail = ['--input', str(path)]
    else:
        path = widget.output/'hybrid_imaging_settings.json'
        path.write_text(json.dumps(asdict(settings), indent=2), encoding='utf-8')
        command_tail = ['--source', str(source), '--imaging-settings-json', str(path)]
    widget.cancel_file = widget.output/'cancel.request'
    widget.active_report = widget.output/'resources.json'
    widget.job_kind = 'hybrid'
    command = ['-m', 'fullwave_prototype.guarded_run', '--report', str(widget.active_report),
        '--max-memory-gib', str(widget.memory.value()), '--timeout-seconds', str(widget.timeout.value()*60),
        '--cancel-file', str(widget.cancel_file), 'hybrid_run', '--output', str(widget.output)]+command_tail
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
    widget.status.setText('Hybrid whole SPM + shared local FEM | '+('re-imaging' if source else 'new solve'))
    widget.process.start(str(interpreter), command)

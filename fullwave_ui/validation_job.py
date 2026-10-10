"""Launch source-locked model validation using the same guarded UI process."""
import codecs
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path

from PyQt5.QtCore import QProcessEnvironment
from fullwave_prototype.hybrid_run import load_sources
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.imaging_run import build_parser
from fullwave_prototype.imaging_options import imaging_settings
from .runtime_paths import application_root, results_root


def start_validation(widget, mode):
    if mode not in ('range', 'mesh', 'power', 'local', 'local_p2', 'ensemble128'):
        raise ValueError('Unknown model validation')
    if mode == 'ensemble128':
        if widget.loaded_whole is None:
            raise ValueError('Load a complete whole-beam haze result before extending statistics')
        report = json.loads((widget.loaded_whole/'whole_haze_report.json').read_text(encoding='utf-8'))
        if report['metrics']['samples'] >= 128:
            raise ValueError('This ensemble already has 128 surfaces')
    elif widget.loaded_hybrid is None:
        raise ValueError('Load a complete hybrid source before validating it')
    else:
        _, metadata = load_sources(widget.loaded_hybrid)
    interpreter = Path(widget.interpreter.text())
    if not interpreter.is_file():
        raise ValueError('Select the isolated FEM interpreter')
    root = application_root()
    widget.output = results_root(root)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    widget.output.mkdir(parents=True, exist_ok=False)
    if mode == 'ensemble128':
        module = 'whole_haze'
        values = dict(parameters=report['parameters'], roughness=report['roughness'],
            imaging_settings=report['imaging_settings'], samples=128,
            coverage=report['plan']['target_power_coverage'], target_relative_ci=.05)
        path = widget.output/'ensemble_input.json'
        path.write_text(json.dumps(values, indent=2), encoding='utf-8')
        tail = ['--input', str(path), '--output', str(widget.output), '--resume', str(widget.loaded_whole)]
    else:
        tail = ['--source', str(widget.loaded_hybrid), '--output', str(widget.output)]
    if mode == 'mesh':
        module = 'surface_mesh_audit'
    elif mode == 'range':
        module = 'hybrid_range_validation'
        p = Parameters(**metadata['parameters'])
        args = build_parser().parse_args(widget.controls.arguments(imaging_only=True,
            incidence_deg=p.incidence_deg))
        args.illumination_length_um = p.beam_length_um
        args.illumination_width_um = p.beam_width_um
        settings = imaging_settings(args, parameters=p)
        path = widget.output/'validation_imaging_settings.json'
        path.write_text(json.dumps(asdict(settings), indent=2), encoding='utf-8')
        tail += ['--imaging-settings-json', str(path)]
    elif mode != 'ensemble128':
        module = 'local_spm_validation'
        tail += ['--orders', '2' if mode == 'local_p2' else '1']
        if mode == 'local_p2':
            tail += ['--direct-solver', 'sparsecholesky']
        if mode == 'power':
            tail += ['--height-scales', '.1', '--power-checks']
    widget.cancel_file = widget.output/'cancel.request'
    widget.active_report = widget.output/'resources.json'
    widget.job_kind = 'whole_haze' if mode == 'ensemble128' else 'model_validation'
    command = ['-m', 'fullwave_prototype.guarded_run', '--report', str(widget.active_report),
        '--max-memory-gib', str(widget.memory.value()), '--timeout-seconds', str(widget.timeout.value()*60),
        '--cancel-file', str(widget.cancel_file), module]+tail
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
    widget.status.setText('Source-locked model check | '+mode)
    widget.process.start(str(interpreter), command)

"""Source-locked validation launch and completed-report UI behavior."""
from dataclasses import asdict
import json
import os
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import pytest
pytest.importorskip('PyQt5')
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication

from fullwave_prototype.parameters import Parameters
from fullwave_ui.widget import FullwaveWidget
from fullwave_ui.results import ResultsTabs


@pytest.fixture
def app():
    return QApplication.instance() or QApplication([])


@pytest.mark.parametrize('mode', ['range', 'mesh', 'power', 'local', 'local_p2'])
def test_validation_guard_and_saved_source_lock(app, tmp_path, mode):
    interpreter = tmp_path/'python.exe'
    interpreter.touch()
    source = tmp_path/'source'
    source.mkdir()
    p = Parameters(beam_profile='gaussian', beam_length_um=3, beam_width_um=4,
                   wavelength_nm=266, incidence_deg=60)
    widget = FullwaveWidget()
    widget.controls.fields['illumination-fit'].setCurrentText('manual')
    widget.interpreter.setText(str(interpreter))
    widget.loaded_hybrid = source
    widget.controls.fields['diameter-nm'].setValue(200)
    widget.controls.fields['wavelength-nm'].setValue(355)
    widget.controls.fields['na2'].setValue(0)
    with patch('fullwave_ui.validation_job.application_root', return_value=tmp_path), \
         patch('fullwave_ui.validation_job.results_root', return_value=tmp_path), \
         patch('fullwave_ui.validation_job.load_sources', return_value=({}, {'parameters': asdict(p)})), \
         patch.object(widget.process, 'start') as start, \
         patch('fullwave_ui.widget.QMessageBox.warning') as warning:
        widget.start_model_validation(mode)
    warning.assert_not_called()
    arguments = start.call_args.args[1]
    assert 'fullwave_prototype.guarded_run' in arguments
    assert '--cancel-file' in arguments and '--max-memory-gib' in arguments
    assert arguments[arguments.index('--source')+1] == str(source)
    assert '--diameter-nm' not in arguments and '--wavelength-nm' not in arguments
    assert widget.job_kind == 'model_validation' and not widget.validation_button.isEnabled()
    if mode == 'range':
        settings = json.loads(Path(arguments[arguments.index('--imaging-settings-json')+1]).read_text())
        assert settings['numerical_apertures'][1] == 0
        assert settings['illumination_length_um'] == 3
        assert settings['illumination_width_um'] == 4
        assert settings['illumination_incidence_deg'] == 60
    elif mode == 'local_p2':
        assert arguments[arguments.index('--direct-solver')+1] == 'sparsecholesky'
    elif mode == 'power':
        assert '--power-checks' in arguments
        assert arguments[arguments.index('--height-scales')+1] == '.1'
    elif mode == 'mesh':
        assert 'surface_mesh_audit' in arguments
    widget.set_running(False)
    assert widget.validation_button.isEnabled()
    widget.close()


def test_fixed_ensemble_extension_is_locked_to_loaded_source(app, tmp_path):
    interpreter = tmp_path/'python.exe'
    interpreter.touch()
    source = tmp_path/'source'
    source.mkdir()
    report = dict(parameters={'diameter_nm': 10}, roughness={'seed': 7},
        imaging_settings={'numerical_apertures': [.5, .5, .5]},
        plan={'target_power_coverage': .99}, metrics={'samples': 16})
    (source/'whole_haze_report.json').write_text(json.dumps(report))
    widget = FullwaveWidget()
    widget.interpreter.setText(str(interpreter))
    widget.loaded_whole = source
    with patch('fullwave_ui.validation_job.application_root', return_value=tmp_path), \
         patch('fullwave_ui.validation_job.results_root', return_value=tmp_path), \
         patch.object(widget.process, 'start') as start:
        widget.start_model_validation('ensemble128')
    arguments = start.call_args.args[1]
    values = json.loads(Path(arguments[arguments.index('--input')+1]).read_text())
    assert values['samples'] == 128 and values['parameters'] == report['parameters']
    assert arguments[arguments.index('--resume')+1] == str(source)
    assert widget.job_kind == 'whole_haze'
    widget.set_running(False)
    report['metrics']['samples'] = 128
    (source/'whole_haze_report.json').write_text(json.dumps(report))
    with patch.object(widget.process, 'start') as start, \
         patch('fullwave_ui.widget.QMessageBox.warning') as warning:
        widget.start_model_validation('ensemble128')
    start.assert_not_called()
    warning.assert_called_once()
    widget.close()


def test_partial_or_nonfinite_validation_preserves_current_result(app, tmp_path):
    report = dict(model='recovered_parent_outer_spm_range_validation',
                  status='complete_with_physical_limits', last_adjacent_complex_target_met=False,
                  cases=[dict(input_power_coverage=.999, collected_spm_power_w=[1, 2, 3])])
    path = tmp_path/'range_report.json'
    path.write_text(json.dumps(report), encoding='utf-8')
    pixmap = QPixmap(64, 64)
    pixmap.fill()
    assert pixmap.save(str(tmp_path/'range_validation.png'))
    widget = FullwaveWidget()
    widget.loaded_hybrid = tmp_path/'original_source'
    assert widget.load_model_validation(tmp_path)
    assert 'NOT met' in widget.status.text()
    assert widget.loaded_hybrid == tmp_path/'original_source'
    old = widget.results.validation_details.toPlainText()
    for status, power in [('running', 1), ('complete_with_physical_limits', float('nan'))]:
        report['status'] = status
        report['cases'][0]['collected_spm_power_w'][0] = power
        path.write_text(json.dumps(report), encoding='utf-8')
        with pytest.raises(ValueError):
            widget.results.load_validation(tmp_path)
        assert widget.results.validation_details.toPlainText() == old
    widget.close()

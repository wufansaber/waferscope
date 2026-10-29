import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from grating_scene.config import SceneConfig
from grating_scene.examples import EXAMPLES, example_folder, example_inputs, example_record
from grating_scene.runtime import default_engine, default_s4, distribution_root
from tools.grating_hybrid_crosscheck import coincident_xy


def test_portable_runtime_resolves_both_gui_and_worker(tmp_path):
    root = tmp_path.resolve()
    native = root/'s4/S4.exe'
    native.parent.mkdir()
    native.touch()
    worker = root/'fem/WaferScope_FEM.exe'
    worker.parent.mkdir()
    worker.touch()
    assert default_engine(root) == worker
    assert default_s4(root) == native
    for executable in (root/'WaferScope.exe', worker):
        with patch('sys.frozen', True, create=True), patch('sys.executable', str(executable)):
            assert distribution_root() == root
            assert default_s4() == native


def test_source_runtime_fallback_does_not_require_an_installed_bundle(tmp_path):
    assert default_engine(tmp_path) == tmp_path/'build/rcwa_env/Scripts/python.exe'
    assert default_s4(tmp_path) == tmp_path/'build/s4_runtime/bin_stripes/S4.exe'


def test_bundled_example_parameters_are_portable(tmp_path):
    key = EXAMPLES[0]['key']
    folder = tmp_path/'examples/grating'/key
    folder.mkdir(parents=True)
    config = SceneConfig()
    (folder/'config.json').write_text(json.dumps(config.to_dict()), encoding='utf-8')
    loaded, camera = example_inputs(tmp_path, key)
    assert loaded == config and camera is None
    assert example_folder(tmp_path, key) == folder
    with pytest.raises(ValueError):
        example_record('../escape')


def test_hybrid_near_comparison_uses_coincident_points_without_phase_fit():
    x, y = np.meshgrid([-.5, 0., .5], [-.5, 0., .5])
    points = np.stack((x, y, np.full_like(x, .02)), axis=-1)
    native = dict(xy_points_um=points, xy_sample_electric=np.ones((3, 3, 3), complex))
    raw = dict(near_x_um=np.arange(4)/4-.5, near_y_um=np.arange(4)/4-.5,
        near_z_um=np.array(.02), near_electric=np.ones((4, 4, 3), complex))
    assert coincident_xy(native, raw)['points'] == 4
    assert coincident_xy(native, raw)['relative_complex_e_l2'] == 0
    raw['near_electric'] *= np.exp(.2j)
    assert coincident_xy(native, raw)['relative_complex_e_l2'] == pytest.approx(2*np.sin(.1))
    raw['near_z_um'] = .03
    with pytest.raises(ValueError, match='Near planes differ'):
        coincident_xy(native, raw)
